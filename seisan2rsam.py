import os, datetime, pprint
import numpy as np
SEISAN_DATA_ROOT = os.path.join("/raid", "data", "seisan")
DB = "DSNC_"

def main(startdate, enddate):
    event_list=[]
    reapath = os.path.join(SEISAN_DATA_ROOT, 'REA', DB)
    years=range(startdate.year,enddate.year+1)
    for year in years:
        #print year
        if year==enddate.year and year==startdate.year:
            months=range(startdate.month,enddate.month+1)
        elif year==startdate.year:
            months=range(startdate.month,13)
        elif year==enddate.year:
            months=range(1,enddate.month+1)
        else:
            months=range(1,13)
        for month in months:
            #print month
            dir=os.path.join(reapath, "%04d" % year, "%02d" % month)
            #print dir
            flist=glob(os.path.join(dir,"%04d*" % year))
            for file in flist:
                % load file
                st = read(file)
                % for each channel
                for tr in st:

                    % detrend data
                    tr.detrend('linear')
                    % compute abs
                    a = abs(tr.data)
                    m = list()
                    % for each minute
                    for i=0:60:length(a)
                        % compute mean
                        % append to array for channel
                        m.append( mean(a[i:i+59]) )

                    % save to channel file
