#%%

####################################### MAIN RUN CYCE ################################

from neuron import h, gui
h.load_file("thalamocort.hoc")

import matplotlib.pyplot as plt

h.vshift_kdr = 5

h.connect_pcs()
h.connect_pvs()
h.connect_fs_gaps()
h.connect_pvs_to_pcs()
h.connect_pcs_to_pvs()
h.runandplot()
#%%
h.quit()
# %%

####################################### raster plot generation ################################

import matplotlib.pyplot as plt

for k in range(200,220):
    plt.plot(h.aprec[k], h.ivecs[k], 'bo')
plt.axis([0, 1000, 200, 220])
plt.show()


for k in range(0,200):
    plt.plot(h.aprec[k], h.ivecs[k], 'bo')
plt.axis([0, 1000, 0, 220])
plt.show()

"""
binvec = []
for i in range(220):
    for k in range(h.aprec[i].size()):
        binvec.append(h.aprec[i].x[k])
plt.hist(binvec, bins=500, range=(0,500))
plt.show()
"""
# %%

####################################### extract spike times ################################

import numpy as np

vmi1 = [[] for i in range(220)]

maxsize = 0
for k in range(220):
    if h.aprec[k].size()>maxsize:
        maxsize = h.aprec[k].size()

print(maxsize)

for k in range(220):
    templist = []
    for i in range(h.aprec[k].size()):
        templist.append(h.aprec[k][i])
    for i in range(maxsize-h.aprec[k].size()):
        templist.append(0)
    vmi1[k] = np.asarray(templist)


np.savetxt("vmi1.txt", vmi1)

# %%

import numpy as np

V_list = [[] for i in range(20)]

for k in range(20):
    for i in range(h.recvec[0].size()-1):
        if k < 10:
            V_list[k].append(h.recvec[k].x[i])
        else:
            V_list[k].append(h.recvec[200+k].x[i])


np.savetxt("V_traces_5shift.txt", np.transpose(np.asarray(V_list)))


# %%

####################################### calculate spike correlations ################################

import numpy as np

def calc_corr(cell1, cell2):
    npts = 10000
    x = np.linspace(0, 50, npts)

    cell1_1 = np.zeros(10000)
    cell2_1 = np.zeros(10000)
    for ap in cell1:
        cell1_1[int(ap*10)] = 1
    for ap in cell2:
        cell2_1[int(ap*10)] = 1

    y1 = cell1_1
    y2 = cell2_1

    lags = np.arange(-npts + 1, npts)
    ccov = np.correlate(y1 - y1.mean(), y2 - y2.mean(), mode='full')
    ccor = ccov / (npts * y1.std() * y2.std())
    """
    fig, axs = plt.subplots(nrows=2)
    fig.subplots_adjust(hspace=0.4)
    ax = axs[0]
    ax.plot(x, y1, 'b', label='y1')
    ax.plot(x, y2, 'r', label='y2')
    ax.set_ylim(-10, 10)
    ax.legend(loc='upper right', fontsize='small', ncol=2)

    ax = axs[1]
    ax.plot(lags, ccor)
    ax.set_ylim(-1.1, 1.1)
    ax.set_ylabel('cross-correlation')
    ax.set_xlabel('lag of y1 relative to y2')
    """
    maxlag = lags[np.argmax(ccor)]
    print("max correlation is at lag %d" % maxlag)
    return ccor 

corr_fin = [[] for i in range(40000)]
counter = 0
for k in range(200):
    for j in range(200):
        print(h.aprec[k].size(), h.aprec[j].size())
        if h.aprec[k].size() > 0 and h.aprec[j].size() > 0:
            corr_fin[counter] = calc_corr(h.aprec[k], h.aprec[j])
        else:
            corr_fin[counter] = np.asarray([0 for l in range(19999)])
        counter += 1
        print(counter)

vmi1 = np.average(corr_fin, axis=0)
lags = np.arange(-10000 + 1, 10000)
plt.plot(lags, vmi1)

np.savetxt("vmi2.txt", vmi1)
# %%
################################### get firing freq of all cells############################

for i in range(220):
    print(h.aprec[i].size())
# %%

################################## get pv subthreshold correlations ###########################

from scipy.stats.stats import pearsonr 


vmi3 = []

for i in range(200,220):
    for j in range(200,220):
        if i != j:
            vmi3.append(pearsonr(h.recvec[i], h.recvec[j])[0])

np.savetxt("vmi3.txt", vmi3)

# %%

################################## get pc subthreshold correlations ###########################

from scipy.stats.stats import pearsonr 


vmi3 = []

for i in range(200):
    for j in range(200):
        if i != j:
            vmi3.append(pearsonr(h.recvec[i], h.recvec[j])[0])

np.savetxt("vmi3.txt", vmi3)

#%%
################################## calc LFP from averages ######################

import numpy as np

avg1 = []

for i in range(200000):
    counter = 0
    for k in range(220):
        counter += h.recvec[k].x[i]
    if i%1000 == 0:
        print(i)
    avg1.append(counter/220)

np.savetxt("LFP_5_2.txt", avg1)
#%%
from scipy.fftpack import fft
import numpy as np

avg1 = np.asarray(avg1)

xf = fft(avg1-np.mean(avg1))
Sxx = 2 * 0.005 ** 2 / 1000 * (xf * xf.conj())
Sxx = Sxx[:int(len(avg1) / 2)]

df = 1 / 1000
fNQ = 1 / 0.005 / 2 
faxis = np.arange(0,fNQ,df) 

plt.semilogy(faxis, Sxx.real)                 # Plot spectrum vs frequency
#plt.xlim([0, 100])                        # Select frequency range
plt.xlim([0, 0.1])                        # Select frequency range
plt.xlabel('Frequency [Hz]')              # Label the axes
plt.ylabel('Power [$\mu V^2$/Hz]')
plt.show()

#%%
from scipy import signal


freqs, psd = signal.welch(avg1, fs=200, nperseg=10000)

plt.figure(figsize=(5, 4))
plt.semilogy(freqs, psd)
plt.title('PSD: power spectral density')
plt.xlabel('Frequency')
plt.ylabel('Power')
plt.tight_layout()

np.savetxt("psd.txt", psd)
# %%

import numpy as np

exc_spikes = []
inh_spikes = []

for i in range(220):
    if i < 200:
        for k in range(h.aprec[i].size()):
            exc_spikes.append([i, h.aprec[i].x[k]])
    else:    
        for k in range(h.aprec[i].size()):
            inh_spikes.append([i,h.aprec[i].x[k]])

np.savetxt("exccells.txt", exc_spikes, fmt='%10.5f')
np.savetxt("inhcells.txt", inh_spikes, fmt='%10.5f')


# %%
############################### calculate LFP    #############################################


import numpy as np
import matplotlib.pyplot as plt
from matplotlib import pyplot as plt
plt.switch_backend('agg')


# 1. read data in vectors

N = 220  # nb of cells to consider
Ne = 200  # nb of excitatory cells
Ni = 20  # nb of inhibitory cells

tmin = 100  # min time (to skip)
tmax = 1000  # max time

dtype = {"names": ["cellid", "time"], "formats": ["i4", "f8"]}
inh_cells = np.loadtxt("inhcells.txt", dtype=dtype)
exc_cells = np.loadtxt("exccells.txt", dtype=dtype)

# adjust time and convert to ms
inh_cells["time"] = inh_cells["time"]
exc_cells["time"] = exc_cells["time"]
# for inhibitory cells, ids start from Ne
inh_cells["cellid"] += Ne



# %%

# 2. draw raster

Nstp = 10  # step cell to draw
tick_size = 5

fig, axes = plt.subplots(2, 1, figsize=(8, 6), sharex=True)

axes[0].plot(exc_cells[::Nstp]["time"], exc_cells[::Nstp]["cellid"], ".", ms=tick_size)
axes[0].plot(inh_cells[::Nstp]["time"], inh_cells[::Nstp]["cellid"], ".", ms=tick_size)



# 3. distribute cells in a 2D grid

xmax = 0.05  # size of the array (in mm)
ymax = 0.05

X, Y = np.random.rand(2, N) * np.array([[xmax, ymax]]).T


# 4. calculate LFP
#
# Table of respective amplitudes:
# Layer   amp_i    amp_e
# deep    -2       -1.6
# soma    30       4.8
# sup     -12      2.4
# surf    3        -0.8
#

dt = 0.1  # time resolution
npts = int(tmax / dt)  # nb points in LFP vector

xe = xmax / 2
ye = ymax / 2  # coordinates of electrode

va = 200  # axonal velocity (mm/sec)
lambda_ = 0.2  # space constant (mm)
dur = 100  # total duration of LFP waveform
nlfp = int(dur / dt)  # nb of LFP pts
amp_e = 0.7  # uLFP amplitude for exc cells
amp_i = -3.4  # uLFP amplitude for inh cells
sig_i = 2.1  # std-dev of ihibition (in ms)
sig_e = 1.5 * sig_i  # std-dev for excitation

# amp_e = -0.16	# exc uLFP amplitude (deep layer)
# amp_i = -0.2	# inh uLFP amplitude (deep layer)

amp_e = 0.48  # exc uLFP amplitude (soma layer)
amp_i = 3  # inh uLFP amplitude (soma layer)

# amp_e = 0.24	# exc uLFP amplitude (superficial layer)
# amp_i = -1.2	# inh uLFP amplitude (superficial layer)

# amp_e = -0.08	# exc uLFP amplitude (surface)
# amp_i = 0.3	# inh uLFP amplitude (surface)

dist = np.sqrt((X - xe) ** 2 + (Y - ye) ** 2)  # distance to  electrode in mm
delay = 10.4 + dist / va  # delay to peak (in ms)
amp = np.exp(-dist / lambda_)
amp[:Ne] *= amp_i
amp[Ne:] *= amp_e

s_e = 2 * sig_e * sig_e
s_i = 2 * sig_i * sig_i


lfp_time = np.arange(npts) * dt


def f_temporal_kernel(t, tau):
    #function defining temporal part of the kernel
    return np.exp(-(t ** 2) / tau)


def calc_lfp(cells):
    #Calculate LFP from cells

    # this is a vectorised computation and as such it might be memory hungry
    # for long LFP series/large number of cells it may be more efficient to calculate it through looping

    spt = inh_cells["time"]
    cid = inh_cells["cellid"]
    print(np.max(cid-400))
    print(len(amp))
    asdf = amp[None, cid-400]
    kernel_contribs = amp[None, cid-400] * f_temporal_kernel(
        lfp_time[:, None] - delay[None, cid-400] - spt[None, :], s_i
    )
    lfp = kernel_contribs.sum(1)
    return lfp


lfp_inh = calc_lfp(inh_cells)
lfp_exc = calc_lfp(exc_cells)

total_lfp = lfp_inh + lfp_exc
total_lfp = lfp_exc

axes[1].plot(lfp_time, total_lfp)
axes[1].set_xlabel("time, ms")
axes[1].set_xlim(0, tmax)

# prettify graph
axes[0].spines["top"].set_visible(False)
axes[0].spines["right"].set_visible(False)
axes[1].spines["top"].set_visible(False)
axes[1].spines["right"].set_visible(False)
fig.savefig('demo.png')
plt.show()

# %%

#np.savetxt("lfp.txt", total_lfp)

import matplotlib.pyplot as plt
import numpy as np

sampling_rate = 10.0

time = np.arange(0, 1000, 1/sampling_rate)

data = np.loadtxt("lfp.txt")

fourier_transform = np.fft.rfft(data)

abs_fourier_transform = np.abs(fourier_transform)

power_spectrum = np.square(abs_fourier_transform)

frequency = np.linspace(0, sampling_rate/2, len(power_spectrum))

plt.plot(frequency, power_spectrum)

# %%
