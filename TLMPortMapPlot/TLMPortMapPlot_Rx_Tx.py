# -*- coding: utf-8 -*-
"""
Created on Mon Dec 18 16:40:42 2023

@author: jmitchell
"""
# imports
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt;

plt.rcParams['font.size'] = 12
import matplotlib.backends.backend_pdf
from matplotlib import cm, colors
from scipy.stats import norm
import os
import glob
import copy
import csv
import json
import time
from pylab import *
# import seaborn as sns
from matplotlib.markers import MarkerStyle
from TLMPortMapPlotFunctions import *

plt.close('all')

TLM_Type = 'Tx' #'Rx'

# set-up
file_path = r'.\2024-07-01_17-16-37_MCR1_Rig1_cal_QR00169_1_20_3820_Bias_0_Override_7_29.50_45C'
if TLM_Type == 'Rx':
    map_tlm_df = pd.read_csv(r'.\20240227_tlm_map_plotter\20221019_TLMCalInputs\Mrk1_S2000_TLM_RX_ArrayGeometry_V20062022_CalInfo_Orig.csv',header=1)
    freq_list = ['19.20']
elif TLM_Type == 'Tx':
    map_tlm_df = pd.read_csv(r'.\20240227_tlm_map_plotter\20221019_TLMCalInputs\Mrk1_S2000_TLM_TX_ArrayGeometry_V20062022_CalInfo.csv',header=1)
    freq_list = ['29.50']#['27.50', '28.00', '28.50', '29.00', '29.50', '30.00', '30.50', '31.00']
align = True
beam_list = [2]
it = 1

# run
for gain_phase in ['gain', 'phase']:

    for delta_pol in [False, True]:

        if delta_pol == True:
            if gain_phase == 'gain':
                vmax_std = 6.0
                v_OP = 8.0
                v_RFA = v_OP * 1.0
                v_OP_step = 0.01
                v_RFA_step = v_OP_step * 1.0
                v_spread_step = 0.01
                tick_step = 1
            if gain_phase == 'phase':
                vmax_std = 90
                v_OP = 360
                v_RFA = v_OP * 1.0  # 3.0
                v_OP_step = 0.01
                v_RFA_step = v_OP_step * 1.0
                v_spread_step = 0.01
                tick_step = 45.0

        if delta_pol == False:
            if gain_phase == 'gain':
                vmax_std = 6.0
                v_OP = 25.0
                v_RFA = v_OP * 1.0
                v_OP_step = 0.01
                v_RFA_step = v_OP_step * 1.0
                v_spread_step = 0.01
                tick_step = 2
            if gain_phase == 'phase':
                vmax_std = 90
                v_OP = 360
                v_RFA = v_OP * 1.0  # 3.0
                v_OP_step = 0.01
                v_RFA_step = v_OP_step * 1.0
                v_spread_step = 0.01
                tick_step = 45.0

        for f_type in ['OP_2', 'RFA']:

            # initialise out_array
            count = 0
            out_array = np.zeros([len(map_tlm_df) * 2, len(freq_list) * 2 * len(beam_list)])


            for beam in beam_list:

                for freq_set in freq_list:
                    # find measurement files and load one of them to make an empty matrix
                    measFiles, files = find_measFiles(file_path, f_type, beam, freq_set, it)
                    meas_info, meas_params, meas_array, meas_frequencies, meas_array_gain, meas_array_phase, paramName, i = load__measFiles(measFiles[0])
                    meas_array_tot = meas_array_gain * 0.0
                    meas_array_list = []

                    # cycle through meas files to make average and stdev arrays
                    for measFile in measFiles:
                        meas_info, meas_params, meas_array, meas_frequencies, meas_array_gain, meas_array_phase, paramName, i = load__measFiles(measFile)
                        if gain_phase == 'gain':
                            meas_array_list.append(meas_array_gain)
                        if gain_phase == 'phase':
                            meas_array_list.append(meas_array_phase)
                    meas_array_av = np.median(meas_array_list, axis=0)
                    meas_array_std = np.std(meas_array_list, axis=0)

                    # initialise figure
                    fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(40, 9))

                    # make plots for OP or RFA
                    if 'OP' in f_type:
                        Z, Z_trim, map_tlm_df = plot_tlm_map(meas_array_av,
                                     f'({f_type}) {gain_phase} (average), N = {len(measFiles)} \n Freq = {freq_set} GHz, Beam {beam}, Delta Pol={delta_pol}, Lens Align={align}',
                                     -v_OP, v_OP, v_OP_step, float(meas_params['f_c']), 0, tick_step=tick_step,
                                     delta_pol=delta_pol,TLM_Type=TLM_Type, axs = axs,meas_frequencies=meas_frequencies, align=align)
                    elif 'RFA' in f_type:
                        Z, Z_trim, map_tlm_df = plot_tlm_map(meas_array_av,
                                     f'({f_type}) {gain_phase} (average), N = {len(measFiles)} \n Freq = {freq_set} GHz, Beam {beam}, Delta Pol={delta_pol}, Lens Align={align}',
                                     -v_RFA, v_RFA, v_RFA_step, float(meas_params['f_c']), 0, tick_step=tick_step,
                                     delta_pol=delta_pol,TLM_Type=TLM_Type, axs = axs,meas_frequencies=meas_frequencies,  align=align)

                    # add to out_array
                    out_array[:, count] = Z

                    # replace the decimal point in the frequency for filenames
                    freq_str = freq_set.replace('.', 'g')

                    # plot the stdev
                    Z, Z_trim, map_tlm_df = plot_tlm_map(meas_array_std,
                                 f'({f_type}) {gain_phase} (stdev), N = {len(measFiles)} \n Freq = {freq_set} GHz, Beam {beam}, Delta Pol={delta_pol}, Lens Align={align}',
                                 0.0, vmax_std, v_spread_step, float(meas_params['f_c']), 1, tick_step=tick_step,
                                 delta_pol=delta_pol,TLM_Type=TLM_Type, axs = axs,meas_frequencies=meas_frequencies,  align=align)

                    # plot the cartesian
                    # for measFile in measFiles:
                    meas_info, meas_params, meas_array, meas_frequencies, meas_array_gain, meas_array_phase, paramName, i = load__measFiles(measFile)
                    col = np.argmin((meas_frequencies - float(freq_set)) ** 2)
                    axs[2].plot(np.linspace(1, len(meas_array), num=len(meas_array)), meas_array_av[:, col])
                    axs[2].set_xlabel('port')
                    axs[2].set_ylabel(f'{gain_phase}')
                    axs[2].set_xlim([0, 500])
                    axs[2].set_xticks(np.linspace(0, 500, num=int(500 / 50) + 1))
                    axs[2].set_title(
                        f'({f_type}) {gain_phase} (average), N = {len(measFiles)} \n Freq = {freq_set} GHz, Beam {beam}, Lens Align={align}')
                    axs[2].grid('on')

                    # add the stdev to the outarray
                    out_array[:, count + 1] = Z

                    # add to count
                    count = count + 2

                    # chck the directory
                    isExist = os.path.exists(f'{file_path}\\analysis')
                    if not isExist:
                        os.makedirs(f'{file_path}\\analysis')

                    # save the figure
                    plt.savefig(
                        f'{file_path}\\analysis\\{f_type}_{freq_str}_beam{beam}_delta{delta_pol}_{gain_phase}.png',
                        dpi=200)

            # make the output excel
            freq_list_av = [s + '_av' for s in freq_list]
            freq_list_std = [s + '_std' for s in freq_list]
            header_list = []
            for beam in beam_list:
                for i in range(len(freq_list)):
                    header_list.append(freq_list_av[i] + f'_beam{beam}');
                    header_list.append(freq_list_std[i] + f'_beam{beam}')
            df = pd.DataFrame(out_array, columns=header_list)
            lens_list = (list(map_tlm_df['Lens no.']) + list(map_tlm_df['Lens no.']))
            lens_list.sort()
            df['Lens no.'] = lens_list
            col = df.pop("Lens no.")
            df.insert(0, col.name, col)
            pol = []
            for i in range(len(map_tlm_df)):
                pol.append('Odd')
                pol.append('Even')
            df['pol'] = pol.copy()
            col = df.pop("pol")
            df.insert(0, col.name, col)
            df.index += 1
            df.to_excel(f'{file_path}\\analysis\\{f_type} {gain_phase}, N = {len(measFiles)}.xlsx', index=True)