"""This module contains helper functions to build different earthquake catalogs"""
import datetime

def georss_header(head):
  georss=''
  # Write xml header
  georss="""<?xml version=\"1.0\" encoding=\"utf-8\"?>
<feed xmlns=\"http://www.w3.org/2005/Atom\" 
  xmlns:georss=\"http://www.georss.org/georss\">
"""
  georss+="   <title>%s</title>\n" % (head['title'])
  georss+="   <subtitle>%s</subtitle>\n" % (head['subtitle'])
  georss+="    <link> href=\"%s\"</link>\n" % (head['link'])
  georss+="    <updated>%sZ</updated>\n" % (datetime.datetime.now().isoformat("T"))
  georss+="    <author>\n      <name>%s</name>\n      <email>%s</email>\n    </author>\n" % (head['name'],head['email'])
  georss+="    <id>%s</id>\n" % (head['id'])
  return georss

def georss_footer():
  return "</feed>"

def xmlcat_header(head):
  xml=''
  # Write xml header
  xml="""<?xml version=\"1.0\" encoding=\"utf-8\"?>
<?xml-stylesheet type="text/xsl" href="okeqcatalog.xsl"?>
<okeqcatalog>
"""
  xml+="   <title>%s</title>\n" % (head['title'])
  xml+="   <subtitle>%s</subtitle>\n" % (head['subtitle'])
  xml+="    <link>%s</link>\n" % (head['link'])
  xml+="    <updated>%sZ</updated>\n" % (datetime.datetime.now().isoformat("T"))
  xml+="    <author>\n      <name>%s</name>\n      <email>%s</email>\n    </author>\n" % (head['name'],head['email'])
  xml+="    <id>%s</id>\n" % (head['id'])
  return xml

def xmlcat_footer():
  return "</okeqcatalog>"

def dict2xml(val,indent):
#This helper function returns an xml string with all keys mapping to xml tags
#except the entry element which defines the container element 
#All keys and values need to be strings!
#indent determines the indent
  istr=" "*indent
  xml="%s<%s>\n" % (istr,val['entry'])
  for key in list(val.keys()):
    if key!='entry':
      xml+="%s%s<%s>%s</%s>\n" % (istr,istr,key,val[key],key)
  xml+="%s</%s>\n" % (istr,val['entry'])
  return xml

def roman_intensity(intensity):
  conversion=['','I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII']
  return conversion[intensity]

def isdst(utctime):
  begin={'2001':datetime.datetime(2001,4,1,8,0),
    '2002':datetime.datetime(2002,4,7,8,0),
    '2003':datetime.datetime(2003,4,6,8,0),
    '2004':datetime.datetime(2004,4,4,8,0),
    '2005':datetime.datetime(2005,4,3,8,0),
    '2006':datetime.datetime(2006,4,2,8,0),
    '2007':datetime.datetime(2007,3,11,8,0),
    '2008':datetime.datetime(2008,3,9,8,0),
    '2009':datetime.datetime(2009,3,8,8,0),
    '2010':datetime.datetime(2010,3,14,8,0),
    '2011':datetime.datetime(2011,3,13,8,0),
    '2012':datetime.datetime(2012,3,11,8,0),
    '2013':datetime.datetime(2013,3,10,8,0),
    '2014':datetime.datetime(2014,3,9,8,0),
    '2015':datetime.datetime(2015,3,8,8,0)}
  end={'2001':datetime.datetime(2001,10,28,8,0),
    '2002':datetime.datetime(2002,10,27,7,0),
    '2003':datetime.datetime(2003,10,26,7,0),
    '2004':datetime.datetime(2004,10,31,7,0),
    '2005':datetime.datetime(2005,10,30,7,0),
    '2006':datetime.datetime(2006,10,29,7,0),
    '2007':datetime.datetime(2007,11,4,7,0),
    '2008':datetime.datetime(2008,11,2,7,0),
    '2009':datetime.datetime(2009,11,1,7,0),
    '2010':datetime.datetime(2010,11,7,7,0),
    '2011':datetime.datetime(2011,11,6,7,0),
    '2012':datetime.datetime(2012,11,4,7,0),
    '2013':datetime.datetime(2013,11,3,7,0),
    '2014':datetime.datetime(2014,11,2,7,0),
    '2015':datetime.datetime(2015,11,1,7,0)}
  flag=False
  year="%i" % (utctime.year)
  if utctime >= begin[year] and utctime < end[year]:
    flag=True
  return flag

def calc_localtime(utc_time):
  if isdst(utc_time):
    ltime= utc_time - datetime.timedelta(hours=5.0)
  else:
    ltime= utc_time - datetime.timedelta(hours=6.0)
  return ltime
