#!/usr/bin/env python
import os
import datetime, pytz
import glob
utc = pytz.timezone('UTC')
start_dt_utc = utc.localize(datetime.datetime(2022, 4, 1))
end_dt_utc = utc.localize(datetime.datetime(2022, 12, 3))
delta = datetime.timedelta(days=1)
run1 = {'dboutroot':'db/dbgood', 'SDSdir':'SDSgood', 'chuckfiledir':'chuckfiles/good', 'inputdir':'miniseed/good'}
run2 = {'dboutroot':'db/dball', 'SDSdir':'SDSall', 'chuckfiledir':'chuckfiles/all', 'inputdir':'miniseed/*'}
runs = [run1, run2]
runs = [run2]
for run in runs:
    outputdir = run['SDSdir']
    if not os.path.isdir(outputdir):
        os.makedirs(outputdir)
    else:
        os.system(f"rm -rf {outputdir}/*")
    os.system(f"rm -rf {run['dboutroot']}*")
    if not os.path.isdir(run['chuckfiledir']):
        os.makedirs(run['chuckfiledir'])
    else:
        os.system(f"rm -rf {run['chuckfiledir']}/*")

    while (start_dt_utc <= end_dt_utc):
        print(start_dt_utc, end="\n")
        ymd = start_dt_utc.strftime("%Y%m%d")
        dbout = f"{run['dboutroot']}{ymd}"
        chuckfile = f"{run['chuckfiledir']}/badblocks{ymd}.msd"
        startepoch = int(start_dt_utc.timestamp())
        endepoch = startepoch + 86400
        allfiles = sorted(glob.glob(f"{run['inputdir']}/*{ymd}*"))
        if len(allfiles)>0:
            allfilesstr = " ".join(allfiles)
            os.system(f"miniseed2days -U -w '%Y/%{{net}}/%{{sta}}/%{{chan}}.D/%{{net}}.%{{sta}}.%{{loc}}.%{{chan}}.D.%Y.%j' -S {outputdir} -C {chuckfile} -s {startepoch} -e {endepoch} {allfilesstr}")
            file_stats = os.stat(chuckfile)
            if file_stats.st_size == 0:
                os.remove(chuckfile)
        start_dt_utc += delta
