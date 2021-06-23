#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  9 09:44:20 2021

@author: thompsong
"""

def tellme(s):
    print(s)
    plt.title(s, fontsize=16)
    plt.draw()

        #a = input('any ley')
        plt.close('all')
        # select start and end points
        N = len(st)
        fig = plt.figure()
        #st.detrend()
        st.trim(starttime=min_stime, endtime=max_etime, pad=True, fill_value=0)
        #fig = st.plot()
        ax = []
        
        for c in range(N):
            tr = st[c]
            if c==0:
                ax.append(fig.add_subplot(N, 1, N-c))
            else:
                ax.append(fig.add_subplot(N, 1, N-c, sharex=ax[0]))
            ax[c].plot(tr.times(), tr.data, "b-")
            #ax.xaxis_date()
            ax[c].set_ylabel(tr.id, rotation=0)
            #ax.set_xlim(datenum(min_stime), datenum(max_etime))
        #fig.autofmt_xdate()
        fig.show()
       
        #st.trim(starttime = x[0], endtime = x[1])

        tellme('Now do a manual zoom if needed. Click to proceed')
        fig.waitforbuttonpress()

        got_two_points = False
        while not got_two_points:
            tellme('Select two corners of zoom, middle mouse button to finish')
            pts = fig.ginput(2, timeout=-1)
            if len(pts) < 2:
                print("user did not select two points")
            else:
                (x0, y0), (x1, y1) = pts
                xmin, xmax = sorted([x0, x1])
                ymin, ymax = sorted([y0, y1])
                #ax.set_xlim(xmin, xmax)
                #ax.set_ylim(ymin, ymax)
                got_two_points = True
            

        tellme('All Done!')
        plt.close('all')