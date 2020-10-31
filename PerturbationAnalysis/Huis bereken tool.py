# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 16:20:35 2020

@author: USER
"""

bod = 360e3

KK = bod*0.02
notaris = 1500
roelie = 2500
verb = 0

totaalprijs = bod + KK + notaris + roelie+ verb

#%% financia

Eigen = 150e3

lening_min = totaalprijs - Eigen

lening = 220e3

deficiet = lening_min - lening
# lening_min = 230e3

aflossing = 0.5*lening
aflossing_pm = aflossing/(30*12)

per_rente = 1.09/100
rente = per_rente * lening

rente_pm = rente/12 + aflossing_pm


print(f"Met een bod van {bod} is de rente per maand $ {rente_pm}, en het deficiet {deficiet}")