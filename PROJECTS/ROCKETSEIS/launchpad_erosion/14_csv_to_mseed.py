#!/usr/bin/env python
# coding: utf-8

from IPython import get_ipython

get_ipython().run_line_magic('run', 'header.ipynb')
#print(paths)

# Parse lookuptable and convert good CSV files to SDS
lookuptableDF = LLE.removed_unnamed_columns(pd.read_csv(paths['lookuptable']))
lookuptableDF.to_csv('lookuptable_backup.csv')
lookuptableDF = lookuptableDF.sort_values(by=['starttime'])

transducersDF = LLE.removed_unnamed_columns(pd.read_csv(paths['transducersCSVfile']))

for index, row in lookuptableDF.iterrows():
    print(f"SCAFFOLD: {index}, {row['sourcefile']}, {row['passed']})      
    df2 = pd.read_csv(os.path.join(paths['CORRECTED'],row['outputfile']))
    if row['passed']: # only convert to SDS the rows that passed the quality metrics
        LLE.convert2mseed(df2, MSEED_DIR, transducersDF)
    else: # SCAFFOLD: what to do with rows that failed to pass
        LLE.convert2mseed_badfile(df2, MSEED_DIR, transducersDF)    


