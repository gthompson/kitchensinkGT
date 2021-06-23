#! /usr/bin/env python
""" location class, Glenn Thompson 2010/04/29 """
    
class location(object):
    
    def __init__(self, lon=0, lat=0, z=0, zdir='up'):
        """ initialise a Location object 
            zdir is down for depth (origin), up for height (station) """
        self.lon = lon
        self.lat = lat
        self.z = z
        self.zdir = zdir

    def __str__(self):
        """ print out the attributes of a Location object """
        print('longitude = %.4f' % self.lon,', latitude = %.4f' % self.lat,', z = %.1f' % self.z, ', z-direction is %s' % self.zdir)
            

if __name__ == "__main__":
    import random
    l = [location(), location(random.random() * 360 - 180, random.random() * 180 - 90, random.random() * 30, 'down')]
    for i in range(len(l)):
        print(l[i].__str__())
    



    


    

