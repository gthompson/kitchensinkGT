# mymode.py
def countlines(name):
    numlines = -1
    fin = openfile(name)
    if fin:
        numlines = len(fin.readlines())
        fin.close()
    return numlines

def countChars(name):
    numChars = -1
    fin = openfile(name)
    if fin:
        numChars = len(fin.read())
        fin.close()
    return numChars

def test(name):
    nl = countlines(name)
    nc = countChars(name)
    print 'Number of lines is %d' % nl
    print 'Number of characters is %d' % nc

def openfile(name):
    try:
        fin = open(name, 'r')
    except IOError:
        print 'Cannot open %s ',name

    return fin
    

    


