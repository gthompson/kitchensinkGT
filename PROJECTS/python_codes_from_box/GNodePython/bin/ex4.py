str1 = """g fmnc wms bgblr rpylqjyrc gr zw fylb. rfyrq ufyr
amknsrcpq ypc dmp. bmgle gr gl zw fylb gq glcddgagclr ylb rfyr'q
ufw rfgq rcvr gq qm jmle. sqgle qrp.rpylqjyrc() gq pcamkkclbcb!"""

str2 = ""

for c in xrange(len(str1)):
    ascii = ord(str1[c])
    if (ascii >= ord('a') and (ascii <= ord('z'))):
        ascii += 2
    if ascii > ord('z'):
        ascii -= 26
    str2 += chr(ascii)

print str2

# use str.translate() !
