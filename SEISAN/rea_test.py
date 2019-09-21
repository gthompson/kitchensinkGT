from rea import REA
r = REA('/raid/data/Montserrat/MASTERING/seisan/REA/MVOE_/2002/01/31-2356-00L.S200201') # create the object
r.read_file() # print lines from the file to screen
print((vars(r))) # print attributes
print(r) # fails because there needs to be an origin. Default behaviour needed.
r.parse_reaheader() # populates attributes
#r.plot()
st = r.to_obspy() # not working!
print(st)
#st.plot()
