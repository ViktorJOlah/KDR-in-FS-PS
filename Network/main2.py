#%%

####################################### MAIN RUN CYCE ################################
import numpy as np
from neuron import h, gui
h.load_file("L5PC.hoc")

import matplotlib.pyplot as plt

h.vshift_kdr = 0

h.connect_pcs()
h.connect_pvs()
h.connect_fs_gaps()
h.connect_pvs_to_pcs()
h.connect_pcs_to_pvs()
#h.runandplot()

def savelfp(num):
    avg1 = []

    for i in range(200000):
        counter = 0
        for k in range(220):
            counter += h.recvec[k].x[i]
        if i%1000 == 0:
            print(i)
        avg1.append(counter/220)

    np.savetxt(f"LFP{num}.txt", avg1)



for i in range(20):
    if i < 10:
        h.vshift_kdr = 0
    else:
        h.vshift_kdr = 10
    h.runandplot()
    savelfp(i)







# %%

def savelfp(num):
    avg1 = []

    for i in range(200000):
        counter = 0
        for k in range(220):
            counter += h.recvec[k].x[i]
        if i%1000 == 0:
            print(i)
        avg1.append(counter/220)

    np.savetxt(f"LFP_shift{num}.txt", avg1)


for i in range(10):
    h.vshift_kdr = i
    h.runandplot()
    savelfp(i)


