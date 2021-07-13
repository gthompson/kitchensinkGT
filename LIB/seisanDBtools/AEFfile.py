import os
from libMVO import correct_nslc
class AEFfile:

    def __init__(self, path=None):
        self.aefrows = list() # each element is an AEFrow dict
        self.trigger_window =  None
        self.average_window = None
        self.path = path.strip()
        if not os.path.exists(path):
            print("%s does not exist" % path)
            return
        try:
            file = open(path,'r') 
            lines = file.readlines()
        except IOError:
            print("Error: The specified file does not exist - %s" % (path))
            raise e

        for line in lines:

            if len(line) < 80:
                continue

            if line[79] == '3' or  (line[79] == ' ' and line[1:5]=="VOLC"):
                #print line[1:10]
                if line[1:5]=="VOLC":
                    if line[1:10]=="VOLC MAIN":
                        pass
                    else:
                        aefrow = AEFfile.parse_aefline(line)
                        if aefrow:
                            self.aefrows.append(aefrow)
                _i = line.find("trigger window")
                if _i > -1:
                    _str1 = line[_i:_i+24]
                    _i2 = _str1.find('=') + _i
                    _i3 = _str1.find('s') + _i
                    trigwinstr = line[_i2+1:_i3].strip()
                    self.trigger_window = float(trigwinstr)
                _i = line.find("average window")
                if _i > -1:
                    _str1 = line[_i:_i+24]
                    _i2 = _str1.find('=') + _i
                    _i3 = _str1.find('s') + _i
                    self.average_window = float(line[_i2+1:_i3])
        return None

    def __str__(self):
        #str = "aeffile: %s %r " % (self.path, os.path.exists(self.path))
        str = "aeffile: %s " % self.path
        #str = "\n\texists: " + str(os.path.exists(self.path))
        if self.trigger_window:
            str += "\n\ttrigger window: %f" % self.trigger_window
        if self.average_window:
            str += "\n\taverage window: %f" % self.average_window
        for aefrow in self.aefrows:
            str += aefrow.__str__()
        return str
            
    def cat(self):
        print(cat(self.path))
        return None

    def parse_aefline(line):
        aefrow = {'station':None, 'channel':None, 'amplitude':None, 'energy':None, 'ssam':None, 'maxf':None}
        a_idx = line[15:22].find('A')+15 # where amplitude info starts
        #try:
        station = line[6:10].strip()
        channel = line[11:14].strip()
        
        aefrow['id'] = '.%s..%s' % (station, channel)
        if station[0:2]=='MB':            
            aefrow['fixed_id'] = correct_nslc(aefrow['id'], 100.0, shortperiod=False)
        else:
            aefrow['fixed_id'] = correct_nslc(aefrow['id'], 100.0, shortperiod=True)
        aefrow['amplitude'] = float(line[a_idx+1:a_idx+9].strip())
        aefrow['energy'] = float(line[a_idx+11:a_idx+19].strip())
        aefrow['ssam'] = AEFfile.parse_F(line,aefrow['energy'],a_idx+21)
        if a_idx < 20:
            aefrow['maxf'] = float(line[73:79].strip())
        #except:
        #    print(aefrow)
        return aefrow
    
    def parse_F(line, energy, startindex):
        F = {"frequency_bands":[0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 30.0], 
             "percentages": [], "energies": []}
        while startindex < 79 and len(F["percentages"])<12:
            ssamstr = line[startindex:startindex+3].strip()
            if not "." in ssamstr:
                thisval = int(ssamstr)
                F["percentages"].append(thisval)
                F["energies"].append(thisval/100.0 * energy)
            startindex += 3
        return F