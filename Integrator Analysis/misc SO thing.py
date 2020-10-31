# -*- coding: utf-8 -*-
"""
Created on Thu Jul  2 16:05:32 2020

@author: Jurriez
"""

import numpy as np

x = [10, 11, 12]
a = [[11, 10, 12], [12, 15, 20], [11, 14, 16]]

x = np.array(x)
a = np.array(a)

# numpy arrays allow the use of element-wise comparison
logic = x <= a
print("logic selection matrix:")
print(logic)

# flag entries that don't fully meet the conditions as dictated by the logic matrix
flags = np.sum(logic,axis=1) != 3

# counter of false entries
c = np.sum(flags)
print (f"final value of counter is {c}")

mean = np.mean(a[flags == False][:,-1])
print (f"found mean of entries is {mean}")



x = [10, 11, 12]
a = [[11, 10, 12], [12, 15, 20], [11, 14, 16]]

x = np.array(x)
a = np.array(a)

lastitems = []
for lst in a:
    if np.all(x <= lst):
        lastitems.append(lst[-1])
    else:
        c = np.sum(x <= lst)
        print(f"found {c} entries smaller than x array in list")
        
print(f"list of last_items: {lastitems}")
mean = np.mean(lastitems)
print(f"mean of last items: {mean}")