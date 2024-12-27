import numpy as np
import pandas as pd
import math
import os
import subprocess
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TkAgg')

filename = 'PVA-Backbone_300K.dump'

data = pd.read_csv(filename, sep = ' ', header = 8, names= ['id', 'mol', 'type', 'x', 'y', 'z', 'xu', 'yu', 'zu', 'nan', 'nan1'])
data = data.sort_values(['mol', 'id'], ascending=[True,True])

id = np.array(data.id)
mol = np.array(data.mol)
xu = np.array(data.xu)
yu = np.array(data.yu)
zu = np.array(data.zu)


angle = []
for i in range(1,len(id)-1):
    if mol[i-1] == mol[i+1]:
        vx1 = xu[i+1] - xu[i]
        vy1 = yu[i+1] - yu[i]
        vz1 = zu[i+1] - zu[i]

        vx2 = xu[i-1] - xu[i]
        vy2 = yu[i-1] - yu[i]
        vz2 = zu[i-1] - zu[i]

        dot = vx1*vx2+vy1*vy2+vz1*vz2
        mag1 = np.sqrt(vx1**2+vy1**2+vz1**2)
        mag2 = np.sqrt(vx2**2+vy2**2+vz2**2)

        step = np.arccos( dot / (mag1*mag2) )

        angle.append(180.0*step/np.pi)

# print(angle)
bins = np.linspace(0,181,181)
count = 0

kbT= 550.*1.38*10**(-23) # J
conv = kbT*(6.022*10**23)/4184.

######################################################
# Simulation Histograms
######################################################
hist = []
for i in range(len(bins)-1):
    for j in angle:
        if bins[i] < j < bins[i+1]:
            count = count+1
    hist.append(count)
    count = 0

bins = bins[0:-1]
# hist = hist[0:-2]

count_total = 0

for i in range(len(hist)):
    count_total = count_total + hist[i]

out = []
for i in range(len(hist)):
    out.append(hist[i]/count_total)

Etotal = []
for i in range(len(bins)):
    if out[i]<1E-3:
        Etotal.append(np.nan)
    else:
        E = -conv*np.log( out[i] /np.sin( (bins[i])/180.0*np.pi ))
        Etotal.append( E )

minE = np.nanmin(Etotal)

Etotal = Etotal - minE

#########################################
# Theoretical Values
######################################################
def skew_gauss(x, amp, x0, sigma, skew):
    _gauss = amp * np.exp(-(x - x0) ** 2. / (2. * sigma ** 2.))
    return _gauss + skew * (x0 - x) / sigma * _gauss * np.sqrt(2 * np.pi)


E1 = skew_gauss(bins,0.045, 171, 7, 0.0)
E2 = skew_gauss(bins,0.018, 127, 6, 0.0)
E3 = skew_gauss(bins,0.002, 94, 7, 0.0)
E4 = skew_gauss(bins,0.003, 140, 24, 0.0)
Egauss = E1 + E2 + E3 + E4

######################################################
# Tabulated
######################################################

data = pd.read_csv('../CG-REAL-ER.PVA', sep = ' ', header = 8, names= ['i', 'angle', 'energy', 'force', 'nan'])

print(data)

angle = np.array(data.angle)
energy = np.array(data.energy)

######################################################
# Plot
######################################################
plt.rcParams.update({'font.size': 20,
                     'lines.linewidth': 3,
                     'font.family':  'sans-serif'})
fig, ax1 = plt.subplots(figsize=(9,7))
ax2 = ax1.twinx()

ax1.set_xlim(60, 180)

ax1.plot(bins, Egauss, linestyle='-', color = 'black')
ax1.plot(bins, out, linestyle='none',marker='.', markersize=14,  markerfacecolor='none', color = 'black')
ax1.set_ylabel('Probability')
ax1.set_ylim(0, 0.1)

ax2.plot(bins, Etotal, linestyle='none',marker='.', markersize=14, markerfacecolor='none', color = 'red')
ax2.plot(angle, energy, color = 'red')
ax2.set_ylabel('Potential [kcal mol$^{-1}$]', color = 'red')
ax2.tick_params(axis='y', color = 'red', labelcolor='red')
ax2.set_ylim(0, 10)

ax1.set_xlabel('CHOH Angle [$^o$]')

######################################################
# Plot
######################################################

filename = 'PVA-ER.dump'

data = pd.read_csv(filename, sep = ' ', header = 8, names= ['id', 'mol', 'type', 'x', 'y', 'z', 'xu', 'yu', 'zu', 'nan', 'nan1'])
data = data.sort_values(['mol', 'id'], ascending=[True,True])

id = np.array(data.id)
mol = np.array(data.mol)
xu = np.array(data.xu)
yu = np.array(data.yu)
zu = np.array(data.zu)


angle = []
for i in range(1,len(id)-1):
    if mol[i-1] == mol[i+1]:
        vx1 = xu[i+1] - xu[i]
        vy1 = yu[i+1] - yu[i]
        vz1 = zu[i+1] - zu[i]

        vx2 = xu[i-1] - xu[i]
        vy2 = yu[i-1] - yu[i]
        vz2 = zu[i-1] - zu[i]

        dot = vx1*vx2+vy1*vy2+vz1*vz2
        mag1 = np.sqrt(vx1**2+vy1**2+vz1**2)
        mag2 = np.sqrt(vx2**2+vy2**2+vz2**2)

        step = np.arccos( dot / (mag1*mag2) )

        angle.append(180.0*step/np.pi)

# print(angle)
bins = np.linspace(0,181,181)
count = 0

######################################################
# Simulation Histograms
######################################################
hist = []
for i in range(len(bins)-1):
    for j in angle:
        if bins[i] < j < bins[i+1]:
            count = count+1
    hist.append(count)
    count = 0

bins = bins[0:-1]
# hist = hist[0:-2]

count_total = 0

for i in range(len(hist)):
    count_total = count_total + hist[i]

out = []
for i in range(len(hist)):
    out.append(hist[i]/count_total)




ax1.plot(bins, out, linestyle='none',marker='^',markersize=14,  markerfacecolor='none', color = 'blue')



plt.show()
