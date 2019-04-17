import argparse
import os

from schicexplorer.utilities import opener
import gzip

def parse_arguments(args=None):

    parser = argparse.ArgumentParser(
        add_help=False
    )

    parserRequired = parser.add_argument_group('Required arguments')

    # define the arguments
    parserRequired.add_argument('--fastq', '-f',
                                help='The fastq files to demultiplex',
                                metavar='list of fastq files to demultiplex',
                                nargs='+',
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


def splitFastq(pFastqFile, pOutputFolder):
    # pass
    file_writer = []
    cell_index = {}
    cell_counter = 0
    cell_index_write = 0
    line_count = 0
    for fastq in pFastqFile:
        fh = opener(fastq)

        # line = fh.readline()
        line = True
        while line:
            # print(line)
            try:
                line = True
                fastq_read = []
                for i in range(0, 16):
                    fastq_read.append(fh.readline()) # 5 and 9
                    line_count += 1
                
                # cell_index = 0
                # print(fastq_read)
                forward_barcode = fastq_read[5].strip().decode("utf-8")
                reverse_barcode = fastq_read[9].strip().decode("utf-8")
                # print(forward_barcode + reverse_barcode)

                if forward_barcode + reverse_barcode in cell_index:
                    cell_index_write = cell_index[forward_barcode + reverse_barcode]

                else:
                    cell_index_write = cell_counter
                    cell_index[forward_barcode + reverse_barcode] = cell_counter
                    cell_counter += 1
                    forward_read_cell = pOutputFolder + '/' + forward_barcode + '_' + reverse_barcode + '_R1.fastq'
                    reverse_read_cell = pOutputFolder + '/' + forward_barcode + '_' + reverse_barcode + '_R2.fastq'

                    file_writer.append([forward_read_cell, reverse_read_cell])
                    

                writeFile(file_writer[cell_index_write][0], fastq_read[:4])
                writeFile(file_writer[cell_index_write][1], fastq_read[12:])

                
            except:
                line = False

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
    splitFastq(args.fastq, args.outputFolder)
