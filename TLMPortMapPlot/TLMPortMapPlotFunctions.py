import matplotlib.backends.backend_pdf
from matplotlib import cm, colors
from scipy.stats import norm
import os
import glob
import copy
import csv
import json
import time
import pandas as pd
from pylab import *
# import seaborn as sns
from matplotlib.markers import MarkerStyle
# definitions
def find_measFiles(path, fileString, beam, freq_set, it):
    files = []
    for root, directories, file in os.walk(path):
        for file in file:
            if (file.endswith(".csv")) == True:
                files.append(os.path.join(root, file))
    measFiles = []
    for i in range(len(files)):
        if fileString in files[i] and 'eam' + str(beam) in files[i] and f'{freq_set}_GHz_4' in files[i] and f'teration_{it}' in files[i]:
            measFiles.append(files[i])
    return measFiles, files

def load__measFiles(file_path):
    if os.path.getsize(file_path) > 2000:
        meas_params = {}
        meas_info = []

        # meas_info, array and measurement frequencies
        with open(file_path, 'r') as file:
            filecontent = csv.reader(file, delimiter=',')
            for row in filecontent:
                meas_info.append(row)
            index_start = [index for index in range(len(meas_info)) if 'barcodes' in meas_info[index]][0] + 2
            meas_info = meas_info[0:index_start]
            meas_array = np.genfromtxt(file_path, delimiter=',', skip_header=index_start)
            meas_array_gain = meas_array[:, ::2]
            meas_array_phase = meas_array[:, 1:][:, ::2]
            meas_frequencies = np.array(meas_info[index_start - 1])[::2].astype(float)


    # meas_params
    for i in range(len(meas_info) - 1):
        if len(meas_info[i]) > 1:
            paramName = meas_info[i][0]
            if paramName[0:2] == '# ':
                paramName = paramName[2:]
            meas_params[paramName] = meas_info[i][1]
            if len(meas_params[paramName]) > 0:
                if meas_params[paramName][0] == ' ':
                    meas_params[paramName] = meas_params[paramName][1:]
    return meas_info, meas_params, meas_array, meas_frequencies, meas_array_gain, meas_array_phase, paramName, i

def plot_tlm_map(array_in, title, cmin, cmax, cstp, f_set, plot_no, tick_step, delta_pol, TLM_Type,axs,meas_frequencies,  align=False, col_map='jet'):


    # map tlm as a df (and roated if needed)
    if TLM_Type == 'Rx':
        map_tlm_df = pd.read_csv(r'.\20240227_tlm_map_plotter\20221019_TLMCalInputs\Mrk1_S2000_TLM_RX_ArrayGeometry_V20062022_CalInfo_Orig.csv',header=1)
    elif TLM_Type == 'Tx':
        map_tlm_df = pd.read_csv(r'.\20240227_tlm_map_plotter\20221019_TLMCalInputs\Mrk1_S2000_TLM_TX_ArrayGeometry_V20062022_CalInfo.csv',header=1)

    if align == True:
        map_tlm_df[' Feed x [mm] shift'] = map_tlm_df[' Feed x [mm]'] - map_tlm_df[' Lens x [mm]']
        map_tlm_df[' Feed y [mm] shift'] = map_tlm_df[' Feed y [mm]'] - map_tlm_df[' Lens y [mm]']
        map_tlm_df[' Dual-Pol Probe rotation [deg]']
        angle = np.array(map_tlm_df[' Feed x [mm] shift']) * 0.0
        if TLM_Type == 'Rx':
            for i in range(len(map_tlm_df[' Feed x [mm] shift'])):
                if map_tlm_df['Lens no.'][i] == 1:
                    angle[i] = 0.0
                    if map_tlm_df['Lens no.'][i] == 2:
                        angle[i] = -77.25
                        if map_tlm_df['Lens no.'][i] == 3:
                            angle[i] = -154.5
        if TLM_Type == 'Tx':
            for i in range(len(map_tlm_df[' Feed x [mm] shift'])):
                if map_tlm_df['Lens no.'][i] == 1:
                    angle[i] = 0.0
                    if map_tlm_df['Lens no.'][i] == 2:
                        angle[i] = -102.5
                        if map_tlm_df['Lens no.'][i] == 3:
                            angle[i] = -205
        map_tlm_df['angle'] = angle
        map_tlm_df['x_rot'] = map_tlm_df[' Feed x [mm] shift'] * cos(map_tlm_df['angle'] * np.pi / 180.0) - map_tlm_df[
            ' Feed y [mm] shift'] * sin(map_tlm_df['angle'] * np.pi / 180.0)
        map_tlm_df['y_rot'] = map_tlm_df[' Feed y [mm] shift'] * cos(map_tlm_df['angle'] * np.pi / 180.0) + map_tlm_df[
            ' Feed x [mm] shift'] * sin(map_tlm_df['angle'] * np.pi / 180.0)
        map_tlm_df['x_new'] = map_tlm_df['x_rot'] + map_tlm_df[' Lens x [mm]']
        map_tlm_df['y_new'] = map_tlm_df['y_rot'] + map_tlm_df[' Lens y [mm]']
        print(map_tlm_df)

    # plot rfics
    if TLM_Type == 'Rx':
        map_rfic = pd.read_csv(r'.\20240227_tlm_map_plotter\20221019_TLMCalInputs\MK1_RX_TLM_RFIC_Patch_Feed_mapping_RF.csv')
    elif TLM_Type == 'Tx':
        map_rfic = pd.read_csv(r'.\20240227_tlm_map_plotter\20221019_TLMCalInputs\MK1_TX_TLM_RFIC_Patch_Feed_Mapping.csv')
    rfics = list(set(list(map_rfic['RFIC Number'])))
    print(rfics)
    for rfic in rfics:
        map_rfic_cut = map_rfic[map_rfic['RFIC Number'] == rfic]
        print(map_rfic_cut)
        patches = map_rfic_cut['Patch Number']
        for lens in [0, 1, 2]:
            x_rfic = [];
            y_rfic = []
            for patch in patches:
                if align == True:
                    x_rfic.append(map_tlm_df['x_new'][patch - 1 + lens * float(len(map_tlm_df)) / 3])
                    y_rfic.append(map_tlm_df['y_new'][patch - 1 + lens * float(len(map_tlm_df)) / 3])
                    axs[plot_no].text(map_tlm_df['x_new'][patch - 1 + lens * float(len(map_tlm_df)) / 3],
                                      map_tlm_df['y_new'][patch - 1 + lens * float(len(map_tlm_df)) / 3],
                                      patch, fontsize=3)
                else:
                    x_rfic.append(map_tlm_df[' Feed x [mm]'][patch - 1 + lens * float(len(map_tlm_df)) / 3])
                    y_rfic.append(map_tlm_df[' Feed y [mm]'][patch - 1 + lens * float(len(map_tlm_df)) / 3])
                    axs[plot_no].text(map_tlm_df[' Feed x [mm]'][patch - 1 + lens * float(len(map_tlm_df)) / 3],
                                      map_tlm_df[' Feed y [mm]'][patch - 1 + lens * float(len(map_tlm_df)) / 3],
                                      patch, fontsize=3)
            axs[plot_no].plot(x_rfic, y_rfic, 'm-', linewidth=1.0, alpha=0.7)

    # select the column of data
    col = np.argmin((meas_frequencies - f_set) ** 2)
    Z = array_in[:, col]
    np.median(Z)
    # odd ports
    Z_trim = Z[::2]
    Z_trim_pol1 = Z_trim.copy()

    # marker for odd ports
    m = MarkerStyle('D', fillstyle='left')
    m._transform.rotate_deg(map_tlm_df[' Dual-Pol Probe rotation [deg]'][i] + 45)
    #print('QQQQQqqqq',rot[i])
    # scatter plot for odd pol
    v = np.linspace(cmin, cmax, int((cmax - cmin) / cstp), endpoint=True)
    cmap_chosen = matplotlib.colormaps.get_cmap(col_map)
    #print(v)
    cntr = axs[plot_no].scatter(map_tlm_df[' Feed x [mm]'],map_tlm_df[' Feed y [mm]'], c=Z_trim, marker=m, s=200, edgecolors='black', linewidths=0.5, cmap=cmap_chosen,vmin=min(v), vmax=max(v), alpha=1.0)

    # even pol
    Z_trim = Z[1::2]
    Z_trim_pol2 = Z_trim.copy()

    # marker for even ports
    m = MarkerStyle('D', fillstyle='right')
    m._transform.rotate_deg(map_tlm_df[' Dual-Pol Probe rotation [deg]'][i] + 45)

    # scatter plot for even pol
    cntr = axs[plot_no].scatter(map_tlm_df[' Feed x [mm]'],map_tlm_df[' Feed y [mm]'] , c=Z_trim, marker=m, s=200, edgecolors='black', linewidths=0.5, cmap=cmap_chosen,vmin=min(v), vmax=max(v), alpha=1.0)

    # scatter plot (on-top) for delta pol
    if delta_pol == True:
        m = MarkerStyle('s')
        m._transform.rotate_deg(map_tlm_df[' Dual-Pol Probe rotation [deg]'][i])
        cntr = axs[plot_no].scatter(map_tlm_df[' Feed x [mm]'],map_tlm_df[' Feed y [mm]'] , c=(Z_trim_pol1 - Z_trim_pol2), marker=m, s=200, edgecolors='black',linewidths=0.5, vmin=min(v), vmax=max(v), cmap=cmap_chosen, alpha=1.0)

    cbar = plt.colorbar(cntr)
    cbar.set_ticks(np.arange(min(v), max(v) + tick_step, tick_step))
    axs[plot_no].set_xlabel('X [mm]');
    axs[plot_no].set_ylabel('Y [mm]')
    axs[plot_no].set_title(title)
    return Z, Z_trim, map_tlm_df