import os, sys, glob
HALDIR = '/media/sda1/KSC'
os.chdir(HALDIR)
NEWTONDIR = '/raid/data/KennedySpaceCenter/duringPASSCAL/REFTEK_DATA'
net = '1R'
if not os.path.exists('DAYS'):
    os.system('ln -s %s/DAYS DAYS' % (NEWTONDIR))
for yyyy in range(2018, 2022):
    for jjj in range(1, 367):

        antelopedb = "antelope/db_BUD_%s_%4d_%03d" % (net, yyyy, jjj)
        if os.path.exists(antelopedb):
            next
        globstr = 'DAYS/%s/*/*.%s.EH?.%4d.%03d' % (net, net, yyyy, jjj)
        print(globstr, '\n')
        mseedfiles = glob.glob('DAYS/%s/*/*.%s..EH?.%4d.%03d' % (net, net, yyyy, jjj))
        print("%d %d %d \n" % (yyyy, jjj, len(mseedfiles) ) )
        if mseedfiles:
            for mseedfile in mseedfiles:
                os.system('miniseed2db %s %s' % (mseedfile, antelopedb) )
