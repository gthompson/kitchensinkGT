import glob, re
pyfiles = glob.glob('*.py')
for pyfile in pyfiles:
    print(' ')
    print(pyfile)
    f = open(pyfile, 'r')
    for line in f:
        if re.search('^class\s', line):
            print('- ',line, end='')
        if re.search('def\s', line):
            print('- ',line, end='')
        if re.search('return\s', line):
            print('   ',line, end='')
