import numpy as np
import pandas as pd
import math
import os
import subprocess
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TkAgg')
# plt.rcParams.update({'font.size': 18,
#                      'lines.linewidth': 2,
#                      'font.family':  'sans-serif'})
# plt.figure(figsize=(8,8))


filename = 'PVA-Backbone_300K.dump'

data = pd.read_csv(filename, sep = ' ', header = 8, names= ['id', 'mol', 'type', 'x', 'y', 'z', 'xu', 'yu', 'zu', 'nan', 'nan1'])
data = data.sort_values(['mol', 'id'], ascending=[True,True])


id = np.array(data.id)
mol = np.array(data.mol)
xu = np.array(data.xu)
yu = np.array(data.yu)
zu = np.array(data.zu)

dist = []
for i in range(0,len(id)-1):
    if mol[i+1] == mol[i]:
        vx = xu[i+1] - xu[i]
        vy = yu[i+1] - yu[i]
        vz = zu[i+1] - zu[i]

        dist.append(np.sqrt(vx**2 + vy**2 + vz**2))

# print(angle)
bins = np.linspace(2.2,3,81)
count = 0

######################################################
# Simulation Histograms
######################################################
hist = []
for i in range(len(bins)-1):
    for j in dist:
        if bins[i] < j < bins[i+1]:
            count = count+1
    hist.append(count)
    count = 0

count_total = 0

for i in range(len(hist)):
    count_total = count_total + hist[i]

out = []
for i in range(len(hist)):
    out.append(hist[i]/count_total)

######################################################
# Theoretical Values
######################################################
Etotal = []
for i in range(len(bins)-1):
    Etotal.append( -np.log(out[i]) ) #

Etotal = Etotal - np.min(Etotal)

Eguess = []
for i in range(len(bins)-1):
    Eguess.append( (bins[i]-2.62)**2*140 )

######################################################
# CG
######################################################

filename = 'PVA-ER.dump'

data = pd.read_csv(filename, sep = ' ', header = 8, names= ['id', 'mol', 'type', 'x', 'y', 'z', 'xu', 'yu', 'zu', 'nan', 'nan1'])
data = data.sort_values(['mol', 'id'], ascending=[True,True])


id = np.array(data.id)
mol = np.array(data.mol)
xu = np.array(data.xu)
yu = np.array(data.yu)
zu = np.array(data.zu)

dist = []
for i in range(0,len(id)-1):
    if mol[i+1] == mol[i]:
        vx = xu[i+1] - xu[i]
        vy = yu[i+1] - yu[i]
        vz = zu[i+1] - zu[i]

        dist.append(np.sqrt(vx**2 + vy**2 + vz**2))

# print(angle)
cgbins = np.array(bins)
count = 0

######################################################
# Simulation Histograms
######################################################
cghist = []
for i in range(len(cgbins)-1):
    for j in dist:
        if cgbins[i] < j < cgbins[i+1]:
            count = count+1
    cghist.append(count)
    count = 0

count_total = 0

for i in range(len(cghist)):
    count_total = count_total + cghist[i]

cgout = []
for i in range(len(cghist)):
    cgout.append(cghist[i]/count_total)

######################################################
# Plot
######################################################
plt.rcParams.update({'font.size': 20,
                     'lines.linewidth': 3,
                     'font.family':  'sans-serif'})
fig, ax1 = plt.subplots(figsize=(9,7))
ax2 = ax1.twinx()

l1,=ax1.plot(bins[:80], out, linestyle='none',marker='.', markersize=14, markerfacecolor='none', color = 'black')
ax1.set_ylabel('Probability')
ax1.set_ylim(0, 0.1)

l3,=ax1.plot(cgbins[:80], cgout, linestyle='none',marker='^', markersize=14, markerfacecolor='none', color = 'blue')
ax1.set_ylabel('Probability')
ax1.set_ylim(0, 0.1)

l2,=ax2.plot(bins[:80], Etotal, linestyle='none',marker='.',markersize=14,  markerfacecolor='none', color = 'red')
ax2.plot(bins[:80], Eguess, color = 'red')
ax2.set_ylabel('Potential [k$_B$T]', color = 'red')
ax2.tick_params(axis='y', color = 'red', labelcolor='red')
ax2.set_ylim(0, 10)

ax1.set_xlabel('CHOH Distance [$\AA$]')

plt.legend([l1, l2, l3],['AA','Fit', 'CG'])
plt.show()
