"""
Test of python PFF implementation against reference MATLAB solution.
"""

import sys
import numpy as np
import scipy.io

import matplotlib.pyplot as plt


sys.path.append('..')

from pff import pff_analysis, plot_pff_results


#################
# Load Reference Data from MATLAB Solution


data = scipy.io.loadmat('../data/Sample_LMS_Data.mat')

# Manually copied from reading the data structure into python and it getting
# scrambled.
params = {'freq' : 153.9, # Hz
          'band' : 10, # Hz
          'order' : 3,
          'tlast' : 2.0
          }

tstart = data['t'][np.argmax(data['x'])] + 0.05
tend = 2.0 # params['tlast']

reference = scipy.io.loadmat('../data/Sample_LMS_Results.mat')


#################
# PFF Calls


freq_rad_s, damp_frac_crit, report_t, report_amp, inter_data = pff_analysis(
            data['t'][:, 0], data['x'], params['freq'],
            tstart, tend, params['band']/params['freq'],
            remove_end=20)


plot_pff_results(freq_rad_s, damp_frac_crit, report_t, report_amp,
                 inter_data, base='')


#################
# Extra Plots to Compare to MATLAB results


fig, axs = plt.subplots(2)

report_t = [reference['Time'][:, 0]] + report_t
freq_rad_s = [reference['Freq'][:, 0]*2*np.pi] + freq_rad_s
damp_frac_crit = [reference['Damp'][:, 0]] + damp_frac_crit
report_amp = [reference['Amp'][:, 0]] + report_amp

report_amp = [np.abs(amp) for amp in report_amp]

min_t = report_t[0][0]
max_t = report_t[0][-1]

label = ['Reference', 'This Code']

for dir_ind in range(len(freq_rad_s)):

    axs[0].plot(report_t[dir_ind], freq_rad_s[dir_ind]/2/np.pi, 'o',
                markerfacecolor='none', label=label[dir_ind])


    axs[1].plot(report_t[dir_ind], damp_frac_crit[dir_ind], 'o',
                markerfacecolor='none', label=label[dir_ind])

    min_t = np.minimum(min_t, report_t[dir_ind].min())
    max_t = np.maximum(max_t, report_t[dir_ind].max())

extra_t = np.pi/freq_rad_s[dir_ind].mean()

for ax in axs:
    ax.tick_params(bottom=True, top=True, left=True, right=True,
                   direction="in")
    
    ax.set_xlim((min_t-extra_t, max_t+extra_t))
    
# only label bottom time axis
axs[0].xaxis.set_tick_params(labelbottom=False)
axs[1].set_xlabel('Time [s]')

axs[0].set_ylabel('Frequency [Hz]')
axs[1].set_ylabel('Fraction Critical Damping')

fig.tight_layout()
fig.subplots_adjust(hspace=0.04)

axs[0].legend()

plt.show()

###############
# Amplitude v. freq + damping plot

fig, axs = plt.subplots(2)


min_t = report_t[0][0]
max_t = report_t[0][-1]

label = ['Reference', 'This Code']

for dir_ind in range(len(freq_rad_s)):

    axs[0].plot(report_amp[dir_ind], freq_rad_s[dir_ind]/2/np.pi, 'o',
                markerfacecolor='none', label=label[dir_ind])


    axs[1].plot(report_amp[dir_ind], damp_frac_crit[dir_ind], 'o',
                markerfacecolor='none', label=label[dir_ind])

    min_t = np.minimum(min_t, report_amp[dir_ind].min())
    max_t = np.maximum(max_t, report_amp[dir_ind].max())

extra_t = np.pi/freq_rad_s[dir_ind].mean()

for ax in axs:
    ax.tick_params(bottom=True, top=True, left=True, right=True,
                   direction="in")
    
    ax.set_xlim((min_t-extra_t, max_t+extra_t))
    
# only label bottom time axis
axs[0].xaxis.set_tick_params(labelbottom=False)
axs[1].set_xlabel('Amplitude')

axs[0].set_ylabel('Frequency [Hz]')
axs[1].set_ylabel('Fraction Critical Damping')

fig.tight_layout()
fig.subplots_adjust(hspace=0.04)

axs[0].legend()

plt.show()


