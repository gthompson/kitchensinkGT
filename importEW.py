from obspy.clients.earthworm import Client
import numpy as np
import matplotlib.pyplot as plt
from obspy.core import stream
client = Client("pubavo1.wr.usgs.gov", 16022)
scnl = list()
scnl.append({'station': 'ACH', 'network': 'AV', 'channel': 'EHE', 'location': ''})
num_scnls = len(scnl)
ind = np.arange(num_scnls)
width = 0.30
fig, ax = plt.subplots()
st = stream.Stream(traces=[])
prev_endtime = -1
while 1:
	response = client.get_availability(scnl[0]['network'], scnl[0]['station'], channel=scnl[0]['channel'])
	#print(response)  
	endtime = response[0][5]
	if prev_endtime == -1:
		prev_endtime = endtime - 600
	if prev_endtime >= endtime:
		continue
	scnl_num = 0
	av600 = list()
	av60 = list()
	avnow = list()
	for this in scnl:
		stnew = stream.Stream()
		try:
			stnew = client.get_waveforms(this['network'], this['station'], this['location'], this['channel'], prev_endtime, endtime)
			stnew.merge(method=1)
			if len(stnew)==1:
				trnew = stnew[0]
				st.append(trnew)
				st.merge(method=1)
				tr = st[scnl_num].copy()
				tr.detrend('linear')
				av600.append(np.mean(np.abs(tr.data))) 
				tr.trim(endtime - 60, endtime)
				av60.append(np.mean(np.abs(tr.data))) 
				tr.trim(endtime - 2.56, endtime)
				avnow.append(np.mean(np.abs(tr.data)))
				print('600s = %f,  60s = %f, 2.56s = %f' % (av600[scnl_num], av60[scnl_num], avnow[scnl_num]))
			else:
				#raise
				continue
		except:
			print("Error getting data")
			av600.append(0)
			av60.append(0)
			avnow.append(0)
		scnl_num+=1
	rects1 = ax.bar(ind, av600, width, color='g')
	rects2 = ax.bar(ind + width, av60, width, color='y')
	rects3 = ax.bar(ind + width * 2, avnow, width, color='r')
		
	ax.set_ylabel('RSAM')
	ax.set_title('RSAM over different time intervals')
	ax.set_xticks(ind + width * 2)
	#ax.set_xticklabels(('AV.ACH.EHE'))
	
	plt.show()
	prev_endtime = endtime

