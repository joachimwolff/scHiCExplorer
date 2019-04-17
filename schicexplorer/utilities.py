import gzip

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