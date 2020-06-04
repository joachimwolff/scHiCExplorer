import gzip
import cooler
import re

def opener(filename):
    """
    Determines if a file is compressed or not
    """

    f = open(filename, 'rb')
    # print("gzip or not?", f.read(2))

    if f.read(2) == b'\x1f\x8b':
        f.seek(0)
        return gzip.GzipFile(fileobj=f)
    else:
        f.seek(0)
        return f

def cell_name_list(pScoolUri):
    try:
        matrices_list = cooler.fileops.list_coolers(pScoolUri)

        # old and non-standard scool format stored all cells in root
        # no '/' in matrices_list and no '/cell/*'

        if not '/' in matrices_list:
            return matrices_list
        # new standard scool format, all cells are stored under '/cell/'
        r = re.compile('/cell/*')
        if '/' in matrices_list and any(r.match(line) for line in matrices_list):
            matrices_list.remove('/')
            return matrices_list
        raise Exception('Wrong data format. Please use a scool file.')
        # ['/', '/cells/cell1', '/cells/cell2', '/cells/cell3']

    except Exception:
        raise Exception('Wrong data format. Please use a scool file.')
        exit(1)
