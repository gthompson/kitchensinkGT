#!/usr/bin/env python
# coding: utf-8

# In[54]:


"""
This notebook explores different functions to compute twos complement, 
and determine which correct method is the most efficient.

The algorithm needs to work for an arbitrary number of bits, expressed
by the NBITS variable. 

The desired output is best illustrated by example:
1. Assume the input is an 8-bit unsigned integer.
2. Values from 0-127 should be unmodified.
3. Values from 128-255 should be replaced by their twos complement.
4. Therefore the sequence [255 1 128 85 171 0] should be replaced by...
5. [-1 1 -128 85 -85 0], i.e. integers in range 128-255 have had 2**8 (256) subtracted
6. 2**8 can also be expressed as (1 << 8) which is binary 1 bit shifted 8 positions left
7. Generally, 2**NBITS can be expressed as 1<<NBITS

Credit for the first three functions goes to Gary Bastin. 
I modified them to allow array input, rather than scalars, 
and to modify the input array, rather than return a new array.

Glenn Thompson 2021/10
"""
# imports
import timeit
import numpy as np

# constants
NBITS = 32
ITERATIONS = 100000
TWONBITSMINUS1 = 1<<(NBITS-1)
TWONBITS = 1<<NBITS

def twos_comp_1(vals, nbits):
    for i in range(len(vals)):
        if (vals[i] & (1 << (nbits - 1))) != 0:
            vals[i] = vals[i] - (1 << nbits)
"""
def twos_comp_2(vals, nbits):
    for i in range(len(vals)):
        if vals[i] < 0: # convert back
            vals[i] += (1 << nbits)
        elif (vals[i] & (1 << (nbits - 1))) != 0:
            vals[i] = vals[i] - (1 << nbits) 
"""          
def twos_comp_2(vals, nbits):
    for i in range(len(vals)):
        if (vals[i] & (1 << (nbits - 1))):
            vals[i] = -= (1 << nbits)  
            
def twos_comp_np(vals, nbits):
    vals[vals & (1<<(nbits-1)) != 0] -= (1<<nbits)
    
def twos_comp_gt(vals, nbits):
    vals[vals>=(1<<(nbits-1))] -= (1<<nbits)
    
def twos_comp_simple(vals, nbits):
    vals[vals>=2**(nbits-1)] -= (2**nbits)
    
def twos_comp_gt_v2(vals, nbits):
    vals[vals>=TWONBITSMINUS1] -= TWONBITS  

"""
def twos_comp_gt_v3(vals, nbits):
    for i in range(len(vals)):
        if (vals[i] & TWONBITSMINUS1) != 0:
            vals[i] = vals[i] - TWONBITS   
"""      
def twos_comp_gt_v3(vals, nbits):
    for i in range(len(vals)):
        if (vals[i] > TWONBITSMINUS1):
            vals[i] -= TWONBITS   
            
def twos_comp_gt_v4(vals, nbits):
    for i in range(len(vals)):
        if (vals[i] & TWONBITSMINUS1):
            vals[i] -= TWONBITS               

# create some values
tot = 0
for c in range(0,NBITS,2):
    tot += 2**c
unsigned_list = [(1<<NBITS)-1, 1, 1<<(NBITS-1), tot, (1<<NBITS)-tot, 0]
#unsigned = np.array([(1<<NBITS)-1, 1, 1<<(NBITS-1), tot, (1<<NBITS)-tot, 0], dtype='i8') #dtype=np.uint64)


# Run tests
print('\nRESULTS:')
for func in [twos_comp_1, twos_comp_2, twos_comp_np, twos_comp_gt, twos_comp_simple, twos_comp_gt_v2, twos_comp_gt_v3, twos_comp_gt_v4]:
    for values in [unsigned_list, np.array(unsigned_list)]:
        try:
            # Test for correct results
            v = values.copy()
            print('\nFunction: ', func)
            print('Input: ', v)
            func(v, NBITS)
            print('Output: ', v)

            # Test for speed
            v = values.copy()
            print('Seconds: ', timeit.timeit('func(v, NBITS)', number=ITERATIONS, globals=globals()))
        except Exception as e:
            print(e)


# In[48]:


# Speed test - see if creating a map is faster
def twos_comp_gt_v5(val):
    if (val > TWONBITSMINUS1):
        val -= TWONBITS 
    return val

print('Seconds: ', timeit.timeit('list(map(twos_comp_gt_v5, values))', number=ITERATIONS, globals=globals()))


# In[53]:


# Speed test - see if creating a lambda function is faster
print('Seconds: ', timeit.timeit('''
ans=[None]*len(values)
for i,x in enumerate(values):
    def res(x): return x-TWONBITS
    ans[i] = res(x)
''', number=ITERATIONS, globals=globals()))


# In[63]:


import numpy.distutils.system_info as sysinfo
print(sysinfo.platform_bits)

import platform
print(platform.architecture())

from platform import python_version
print(python_version())
print(np.version.version)


# In[ ]:




