import argparse
import os
import gzip
import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)

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
                                help='The fastq files to demultiplex of Nagano 2017 "Cell cycle dynamics of chromosomal organisation at single-cell resolution"'
                                       'GEO: GSE94489. Files need to have four FASTQ lines per read:'
                                       '/1 forward; /2 barcode forward; /3 bardcode reverse; /4 reverse read',
                                metavar='list of fastq files to demultiplex',
                                nargs='+',
                                required=True)
    parserRequired.add_argument('--barcodeFile', '-b',
                                help='The fastq files to demultiplex',
                                metavar='list of fastq files to demultiplex. Use GSE94489_README.txt file.',
                                required=True)
    parserRequired.add_argument('--srrToSampleFile', '-s',
                                help='The mappings from SRR number to sample id as given in the barcode file.',
                                metavar='',
                                required=True)
    parserRequired.add_argument('--outputFolder', '-o',
                                help='Path of folder to save the demultiplexed files',
                                metavar='FOLDER',
                                required=False,
                                default='demultiplexed')
    return parser

def writeFile(pFileName, pReadArray):
    
    with open(pFileName, 'ab') as f:
        # f.write(content)
        for read in pReadArray:
            f.write(read)

def readSrrToSampleFile(pFileName):
    sample_to_srr = {}
    count = 0
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
            except:
                line = False
    return sample_to_srr

def readBarcodeFile(pFileName):
    barcode_sample = {}
    sample_run_to_individual_run = {}
    count = 0
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
            except:
                line = False
    return barcode_sample, sample_run_to_individual_run

def splitFastq(pFastqFile, pOutputFolder, pBarcodeSampleDict, pSampleToIndividualSampleDict, pSrrToSampleDict):
    # pass
    file_writer = []
    cell_index = {}
    cell_counter = 0
    cell_index_write = 0
    line_count = 0

    
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


        fh = opener(fastq)

        # line = fh.readline()
        line = True
        while line:
            # if line_count > 2000000:
            #     break
            # print(line)
            try:
                line = True
                fastq_read = []
                for i in range(0, 16):
                    fastq_read.append(fh.readline()) # 5 and 9
                line_count += 16
            except:
                line = False    
                # cell_index = 0
                # print(fastq_read)
            forward_barcode = fastq_read[5].strip().decode("utf-8")
            reverse_barcode = fastq_read[9].strip().decode("utf-8")
            # print(forward_barcode + reverse_barcode)
            sample_name = None
            for sample_ in pSampleToIndividualSampleDict[sample]:
                if sample_ in pBarcodeSampleDict:
                    if forward_barcode + reverse_barcode in pBarcodeSampleDict[sample_]:
                        sample_name = sample_
            # if forward_barcode + reverse_barcode in 
            if sample_name is None:
                continue

            if forward_barcode + reverse_barcode + '_' + sample_name in cell_index:
                cell_index_write = cell_index[forward_barcode + reverse_barcode + '_' + sample_name]

            else:
                cell_index_write = cell_counter
                cell_index[forward_barcode + reverse_barcode + '_' + sample_name] = cell_counter
                cell_counter += 1
                forward_read_cell = pOutputFolder + '/' +  sample_name + '_' + forward_barcode + '_' + reverse_barcode + '_R1.fastq'
                reverse_read_cell = pOutputFolder + '/' + sample_name + '_' + forward_barcode + '_' + reverse_barcode + '_R2.fastq'

                file_writer.append([forward_read_cell, reverse_read_cell])
                

            writeFile(file_writer[cell_index_write][0], fastq_read[:4])
            writeFile(file_writer[cell_index_write][1], fastq_read[12:])

                
            

            # print(fh.readlines())

    # for cells in file_writer:
    #     cells[0].close()
    #     cells[1].close()
    
    print('line count {}'.format(line_count))

def main(args=None):

    args = parse_arguments().parse_args(args)
    # print('hello world from scHiCExplorer')
    if not os.path.exists(args.outputFolder):
        try:
            os.makedirs(args.outputFolder)
        except OSError as exc: # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise
    barcode_sample_dict, sample_to_individual_sample_dict = readBarcodeFile(args.barcodeFile)

    log.debug('barcode_sample_dict {}'.format(barcode_sample_dict))
    log.debug('sample_to_individual_sample_dict {}'.format(sample_to_individual_sample_dict))

    srr_to_sample_dict = readSrrToSampleFile(args.srrToSampleFile)
    log.debug('srr_to_sample_dict {}'.format(srr_to_sample_dict))

    splitFastq(args.fastq, args.outputFolder, barcode_sample_dict, sample_to_individual_sample_dict, srr_to_sample_dict)
