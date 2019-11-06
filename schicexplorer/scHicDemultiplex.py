import argparse
import os
import gzip
import shutil
import errno
from multiprocessing import Process, Queue
import time
import logging
log = logging.getLogger(__name__)

from schicexplorer.utilities import opener


def parse_arguments(args=None):

    parser = argparse.ArgumentParser(
        add_help=False
    )

    parserRequired = parser.add_argument_group('Required arguments')

    # define the arguments
    parserRequired.add_argument('--fastq', '-f',
                                help='The fastq files to demultiplex of Nagano 2017 Cell cycle dynamics of chromosomal organization at single-cell resolution'
                                'GEO: GSE94489. Files need to have four FASTQ lines per read:'
                                '/1 forward; /2 barcode forward; /3 bardcode reverse; /4 reverse read',
                                metavar='list of fastq files to demultiplex',
                                nargs='+',
                                required=True)
    parserRequired.add_argument('--barcodeFile', '-bf',
                                help='The fastq files to demultiplex',
                                metavar='list of fastq files to demultiplex. Use GSE94489_README.txt file.',
                                required=True)
    parserRequired.add_argument('--srrToSampleFile', '-s',
                                help='The mappings from SRR number to sample id as given in the barcode file.',
                                required=True)
    parserRequired.add_argument('--outputFolder', '-o',
                                help='Path of folder to save the demultiplexed files',
                                metavar='FOLDER',
                                required=False,
                                default='demultiplexed')
    parserRequired.add_argument('--threads',
                                help='Number of threads. Using the python multiprocessing module.',
                                required=False,
                                default=4,
                                type=int)
    parserRequired.add_argument('--bufferSize', '-bs',
                                help='Number of lines to buffer in memory, if full, write the data to disk.',
                                required=False,
                                default=20e6,
                                type=int)
    return parser


def writeFile(pFileName, pReadArray):

    with gzip.open(pFileName, 'ab') as f:
        for read in pReadArray:
            f.write(read)


def readSrrToSampleFile(pFileName):
    sample_to_srr = {}
    with open(pFileName, 'rb') as file:
        line = file.readline()
        while line:
            try:
                line = file.readline()
                line_ = line.strip().decode("utf-8").split('\t')
                if line_[0] in sample_to_srr:
                    continue
                else:
                    sample_to_srr[line_[0]] = line_[-1]
            except Exception:
                line = False
    return sample_to_srr


def readBarcodeFile(pFileName):
    barcode_sample = {}
    sample_run_to_individual_run = {}
    with open(pFileName, 'rb') as file:
        line = file.readline()
        while line:
            try:
                line = file.readline()
                line_ = line.strip().decode("utf-8").split('\t')
                if line_[-1] in barcode_sample:
                    barcode_sample[line_[-1]].append(line_[2] + line_[3])
                else:
                    barcode_sample[line_[-1]] = [line_[2] + line_[3]]
                if line_[0] in sample_run_to_individual_run:
                    if line_[-1] not in sample_run_to_individual_run[line_[0]]:
                        sample_run_to_individual_run[line_[0]].append(line_[-1])
                else:
                    sample_run_to_individual_run[line_[0]] = [line_[-1]]
            except Exception:
                line = False
    return barcode_sample, sample_run_to_individual_run


def reverseComplement(pString):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    pString = pString[::-1]
    new_string = ''
    for i in pString:
        new_string += complement[i]
    return new_string


def compressFiles(pFileNameList, pQueue):
    for cells in pFileNameList:
        for cell_ in cells:
            with open(cell_, 'rb') as f_in:
                with gzip.open(cell_ + '.gz', 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            os.remove(cell_)
    pQueue.put('Done')


def writeSampleFiles(pFileNameList, pBuffer, pQueue):
    for i, cells in enumerate(pFileNameList):
        for j, cell_ in enumerate(cells):
            writeFile(cell_, pBuffer[i][j])
    pQueue.put('Done')


def handleCompressingMulticore(pFileNameList, pBuffer, pThreads):
    filesPerThread = len(pFileNameList) // pThreads
    if filesPerThread == 0:
        writeSampleFiles(pFileNameList, pBuffer, None)
    else:
        queue = [None] * pThreads
        process = [None] * pThreads
        file_list_sample = [None] * pThreads
        buffer_sample = [None] * pThreads

        all_data_collected = False

        for i in range(pThreads):

            if i < pThreads - 1:
                file_list_sample = pFileNameList[i * filesPerThread:(i + 1) * filesPerThread]
                buffer_sample = pBuffer[i * filesPerThread:(i + 1) * filesPerThread]
            else:
                file_list_sample = pFileNameList[i * filesPerThread:]
                buffer_sample = pBuffer[i * filesPerThread:]

            queue[i] = Queue()
            process[i] = Process(target=writeSampleFiles, kwargs=dict(
                pFileNameList=file_list_sample,
                pBuffer=buffer_sample,
                pQueue=queue[i]
            )
            )

            process[i].start()

        while not all_data_collected:
            for i in range(pThreads):
                if queue[i] is not None and not queue[i].empty():
                    _ = queue[i].get()
                    process[i].join()
                    process[i].terminate()
                    process[i] = None

            all_data_collected = True

            for i in range(pThreads):
                if process[i] is not None:
                    all_data_collected = False
            time.sleep(1)


def splitFastq(pFastqFile, pOutputFolder, pBarcodeSampleDict, pSampleToIndividualSampleDict, pSrrToSampleDict, pThreads, pBufferSize):
    file_writer = []
    cell_index = {}
    cell_counter = 0
    cell_index_write = 0
    line_count = 0
    lines_out_buffer = []
    log.debug('len(pFastqFile) {}'.format(len(pFastqFile)))
    for fastq in pFastqFile:
        if fastq.split('/')[-1].split(".")[0] in pSrrToSampleDict:
            sample = pSrrToSampleDict[fastq.split('/')[-1].split(".")[0]]
        else:
            log.warning('No sample known with SRR: {}'.format(fastq.split('/')[-1].split(".")[0]))
            continue
        log.debug('sample {}'.format(sample))
        if sample not in pSampleToIndividualSampleDict:
            log.warning('No sample known with SRR: {}'.format(fastq.split('/')[-1].split(".")[0]))
            continue
        log.debug('pSampleToIndividualSampleDict[sample]: {}'.format(pSampleToIndividualSampleDict[sample]))
        for sample_ in pSampleToIndividualSampleDict[sample]:
            if sample_ in pBarcodeSampleDict:
                log.debug(' sample {} pBarcodeSampleDict[sample_] {}'.format(sample_, pBarcodeSampleDict[sample_]))
        fh = opener(fastq)

        line = True
        while line:
            fastq_read = []
            for i in range(0, 16):
                line_to_append = fh.readline()
                if not line_to_append:
                    line = False
                    log.debug('end of file! line: {}'.format(line_count))
                    break
                fastq_read.append(line_to_append)  # 5 and 9
            if not line:
                break
            line_count += 16
            forward_barcode = fastq_read[5].strip().decode("utf-8")
            reverse_barcode = fastq_read[9].strip().decode("utf-8")
            sample_name = None
            for sample_ in pSampleToIndividualSampleDict[sample]:
                if sample_ in pBarcodeSampleDict:
                    if forward_barcode + reverse_barcode in pBarcodeSampleDict[sample_]:
                        sample_name = sample_
                        break

            # use reverse complement of second barcode, reason is an error of NextSeq
            # https://bitbucket.org/tanaylab/schic2/src/68d7972f64ac2fd32b7c31c5041b39a7176bf14d/map3c/split_barcodes?at=default&fileviewer=file-view-default#split_barcodes-40
            if sample_name is None:
                reverse_barcode = reverseComplement(reverse_barcode)
                for sample_ in pSampleToIndividualSampleDict[sample]:
                    if sample_ in pBarcodeSampleDict:
                        if forward_barcode + reverse_barcode in pBarcodeSampleDict[sample_]:
                            sample_name = sample_
                            break

            if sample_name is None:
                continue

            if forward_barcode + reverse_barcode + '_' + sample_name in cell_index:
                cell_index_write = cell_index[forward_barcode + reverse_barcode + '_' + sample_name]
                lines_out_buffer[cell_index_write][0].extend(fastq_read[:4])
                lines_out_buffer[cell_index_write][1].extend(fastq_read[12:])
            else:
                cell_index_write = cell_counter
                cell_index[forward_barcode + reverse_barcode + '_' + sample_name] = cell_counter
                cell_counter += 1
                forward_read_cell = pOutputFolder + '/' + sample_name + '_' + forward_barcode + '_' + reverse_barcode + '_R1.fastq.gz'
                reverse_read_cell = pOutputFolder + '/' + sample_name + '_' + forward_barcode + '_' + reverse_barcode + '_R2.fastq.gz'

                file_writer.append([forward_read_cell, reverse_read_cell])
                lines_out_buffer.append([fastq_read[:4], fastq_read[12:]])

            if line_count % pBufferSize == 0 or line is False:
                buffered_elements = 0
                for lines_buffered in lines_out_buffer:
                    buffered_elements += len(lines_buffered[0])
                    buffered_elements += len(lines_buffered[1])

                if buffered_elements > pBufferSize or line is False:
                    handleCompressingMulticore(file_writer, lines_out_buffer, pThreads)

                    lines_out_buffer = None
                    lines_out_buffer = []
                    cell_index = None
                    cell_index = {}
                    file_writer = None
                    file_writer = []
                    cell_counter = 0
        fh.close()


def main(args=None):

    args = parse_arguments().parse_args(args)
    if not os.path.exists(args.outputFolder):
        try:
            os.makedirs(args.outputFolder)
        except OSError as exc:  # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise
    barcode_sample_dict, sample_to_individual_sample_dict = readBarcodeFile(args.barcodeFile)

    srr_to_sample_dict = readSrrToSampleFile(args.srrToSampleFile)

    splitFastq(args.fastq, args.outputFolder, barcode_sample_dict, sample_to_individual_sample_dict, srr_to_sample_dict, args.threads, args.bufferSize)
