#!/usr/bin/env python
import pandas as pd
import os
def fix_df(df):
    good_columns = []
    for thiscol in df.columns:
        if not 'ntitle' in thiscol:
            if not 'Unname' in thiscol:
                good_columns.append(thiscol)
    df2 = df[good_columns] # subset to correct columns
    return df2

def df2csv_without_index(df, csvfile):
    df = df.reset_index()  
    df.drop(df.filter(regex="Unname"),axis=1, inplace=True)
    df.to_csv(csvfile, index=False)

    
#full_unchecked_csv = os.path.join(SEISAN_DATA, 'MachineLearning', SEISAN_DB, 'labelling', '%scatalog.csv' % SEISAN_DB)
#full_unchecked_csv = os.path.join('/Volumes/shareddrive/thompsong', 'MachineLearning', \
#                                  SEISAN_DB, 'labelling', '%s11_master_catalog.csv' % SEISAN_DB)

full_unchecked_csv = '/home/thompsong/DATA/MVO/MachineLearning/MVOE_/labelling/MVOE_11_master_catalog.csv'
full_unchecked_cat = pd.read_csv(full_unchecked_csv)  
unchecked = fix_df(full_unchecked_cat)
print('unchecked has %d rows' % len(unchecked) )  
unchecked2 = unchecked.copy()

old_checked_csv = "../CSVfiles/catalog_unique.csv"
old_checked_csv = "11_merged_catalog.csv"
old_checked_cat = pd.read_csv(old_checked_csv)
checked = fix_df(old_checked_cat)
print('checked has %d rows' % len(checked) )  
checked2 = checked[['path', 'quality','subclass','D', 'R',
       'r', 'e', 'l', 'h', 't', 'new_subclass', 'weight', 'checked', 'split',
       'delete', 'ignore']]

# MERGE
unchecked2.update(checked2)
num_checked = sum(unchecked2['checked']==True)
print('After merging there are %d checked rows out of %d total' % (num_checked,len(unchecked2.index) ) )

# WRITE
if len(unchecked)==len(unchecked2):    
    df2csv_without_index(unchecked2, '11_remerged_catalog.csv')
    # manually copy this to ~/DATA/MVO/MachineLearning/SEISAN_DB/original/merged_catalog.csv
print('Done')
