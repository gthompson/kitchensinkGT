#!/bin/bash
conda activate iceweb
SRC='/home/thompsong/work/PROJECTS/KSC_EROSION/src'
#python $SRC/14_csv_to_mseed.py
#python $SRC/16_miniseed2days_to_SDS.py
#python $SRC/17_make_daily_antelope_dbs.py
#python $SRC/18_sds_availability_plot.py /home/thompsong/work/PROJECTS/KSC_EROSION/SDSgood
#python $SRC/18_sds_availability_plot.py /home/thompsong/work/PROJECTS/KSC_EROSION/SDSall
python $SRC/25_SDS_to_EventFiles.py
python $SRC/30_EventFiles_to_website.py
#python $SRC/19_KSC_sds_to_rsam.py


