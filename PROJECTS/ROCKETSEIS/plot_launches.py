# Call with an optional date like 2022-03-30 to start at a particular date
import os, sys
import pandas as pd
import matplotlib.pyplot as plt
HOME = os.getenv('HOME')

# Import file
xlsdir = os.path.join(HOME, 'Dropbox', 'PROFESSIONAL', 'RESEARCH', \
                        '3_Project_Documents', 'NASAprojects', \
                        '201602 Rocket Seismology', 'tracking_launches_recorded')
xlsfile = os.path.join(xlsdir, 'All_KSC_Rocket_Launches.xlsx')
df = pd.read_excel(xlsfile, sheet_name='launches')

# Process dataframe
dfall = df[['Date', 'Rocket_Payload','SLC']] # subset to critical columns
dfall['SLC'] = dfall['SLC'].astype(str) # turn numeric SLC like '40' to string
print(sys.argv)
if len(sys.argv)>0:
    startdate = sys.argv[1]
dfall = dfall[dfall['Date']>startdate]
dfall['Date']=pd.to_datetime(dfall['Date']) # convert Date string to datetime

# Convert datetime to decimal year
import datetime
def year_fraction(date):
    start = datetime.date(date.year, 1, 1).toordinal()
    year_length = datetime.date(date.year+1, 1, 1).toordinal() - start
    return date.year + float(date.toordinal() - start) / year_length
decyear=[]
for i in range(len(dfall.index)):
    decyear.append(year_fraction(dfall.iloc[i]['Date']))
dfall['decyear']=decyear

# Count occurrences of each SLC
occur = dfall.groupby(['SLC']).size()
print(occur)

# Add column of count for each SLC - easier to plot - probably some better pandas solution
for thisSLC in sorted(occur.keys()):
    count = 0
    count_list = []
    for i in range(len(dfall.index)):
        if dfall.iloc[i]['SLC']==thisSLC:
            count+=1
        count_list.append(count)
    dfall[thisSLC]=count_list

# Plot cumulative launches per SLC against decimal year
plt.rcParams["figure.figsize"] = [19.2, 10.8]
plt.rcParams["figure.autolayout"] = True
plt.rc('font', size=24)
plt.style.use('dark_background')
COLOR = 'white'
plt.rcParams['text.color'] = COLOR
plt.rcParams['axes.labelcolor'] = COLOR
plt.rcParams['xtick.color'] = COLOR
plt.rcParams['ytick.color'] = COLOR
fig = plt.figure(figsize=(19.2,10.8))
fig.set_dpi(100)
ax = dfall.plot.area(x='decyear',y=occur.keys())
ax.legend(loc='upper left')
ax.set_title('Cumulative number of launches by launch complex')
ax.set_xlabel('Year')
ax.set_ylabel('Number of launches')
plt.savefig(os.path.join(xlsdir,'cum_launches_by_SLC_vs_date.png'), transparent=False, dpi=100)

# Count occurrences of each rocket type
rockets = ['Atlas V', 'Delta IV', 'Falcon 9', 'Falcon Heavy', 'SLS']
for thisRocket in rockets:
    count = 0
    count_list = []
    for i in range(len(dfall.index)):
        try:
            if thisRocket in dfall.iloc[i]['Rocket_Payload']:
                count+=1
            
        except:
            pass
        count_list.append(count)
    dfall[thisRocket]=count_list
    
fig2 = plt.figure(figsize=(19.2,10.8))
fig2.set_dpi(100)
ax2 = dfall.plot.area(x='decyear',y=rockets)
ax2.legend(loc='upper left')
ax2.set_title('Cumulative number of launches by rocket')
ax2.set_xlabel('Year')
ax2.set_ylabel('Number of launches')
plt.savefig(os.path.join(xlsdir,'cum_launches_by_rocket_vs_date.png'), transparent=False, dpi=100)    
