import math
from scipy import optimize
from matplotlib.ticker import AutoMinorLocator
from matplotlib import gridspec
import matplotlib.ticker as ticker
import random
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import statistics
from lmfit import models
from scipy.signal import find_peaks_cwt
import peakutils

pat = #'..../medidas/'


######################################################################
#crate a list of the directories
######################################################################
dir_list = os.listdir(pat)

listdir=[]
for d in dir_list:
    listdir.append(pat + d)

listdir.sort()

tb_lisdir = pd.DataFrame(data=listdir)

#f1 = open('tb_lisdir.txt', 'w')
#tb_lisdir[0].to_csv(f1, sep='\t')
#f1.close()

#
print('path:', pat)

def dirplot(path):
    plotDir = path + 'plots/'
    if not os.path.exists(plotDir):
        os.makedirs(plotDir)
    return plotDir
def dirdata(path):
    pathTab = path + 'data/'
    if not os.path.exists(pathTab):
        os.makedirs(pathTab)
    return pathTab

def files_list(path, extension):
    fills=[]
    for i in os.listdir(path):
        if i.endswith(extension):
            fills.append(path + str(i))
    fills.sort()
    return fills

def read_infos(file, info=''): #info = 'T' for Temperature; 'time'; 'I0'; 'RingCurr'
    if info == 'T':
        infos=np.loadtxt(file, skiprows=1, usecols=(2,), unpack=True)
    elif info == 'time':
        infos = []
#            infos.append(0)
        with open(file) as infofile:
            for line in infofile:
                if line != '\n':
                    if (line.find("#")==-1):
                        if (line.find("i")==-1):
                            infos.append(float(line.split()[5]))
    elif info == 'I0':
        infos=np.loadtxt(file, skiprows=1, usecols=(3,), unpack=True)
    elif info == 'RingCurr':
        infos=np.loadtxt(file, skiprows=1, usecols=(3,), unpack=True)
    else:
        print('info = T for Temperature; time; I0; RingCurr')
    return infos

def read_files(files):
    diffra_map=[]

    for i in range(0,len(files)):
        diffra_map.append(np.loadtxt(files[i],skiprows=0,usecols=(1,),unpack=True))

    return diffra_map

def q2theta(q, wavelenght): #A-1 0.15406
    two_theta=[]
    for n in range (0,len(q)):
        x=math.degrees(10*2*np.arcsin((q[n]*wavelenght)/(4*math.pi))) #*10 because q is in A-1
        two_theta.append(x)
    return two_theta


def normBaseline(tab):
    baseline=[]
    for b in range(0,len(tab)):
        baseline.append(peakutils.baseline(tab[b]))
    return baseline
def subBaseline(baseline, tab):
    sub=[]
    for a in range(0, len(baseline)):
        sub.append(tab[a]-baseline[a])
    return sub

def plot_map_temp(time_min, Temperature, sample_name, diffra_map, two_theta_start, two_theta_final, plotDir, save_eps=False,save_png=False):
    fig=plt.figure(figsize=(8,8))
    #-------------------------------------plot temperature vs time
    gs = gridspec.GridSpec(2, 1, height_ratios=[1, 7])
    ax0 = plt.subplot(gs[0])
    ax0.plot(time_min, Temperature, linewidth=1,linestyle='-',  color='C0', marker='o',markersize=3,markerfacecolor='C0')
    ax0.axvline(x=time_min[40], c='black', linestyle='dotted')
    ax0.annotate('Anti-solvent', xy=(time_min[40], 6000), xytext=(time_min[60], 2500),
            arrowprops=dict(arrowstyle="-|>",
                            connectionstyle="angle3,angleA=0,angleB=-90"))
    plt.xlim(0, time_total_min)
    plt.ylim(-1000, 6800)
    plt.ylabel('rpm', size=16)
    plt.xticks(())
#    plt.title(sample_name)
    plt.yticks(fontsize=16)

    fig.subplots_adjust(left=0.15, right=0.85, hspace=0.03)
    #------------------------------------ plot map#------------------------------------
    ax1 = plt.subplot(gs[1])
    cont=ax1.imshow(diffra_map, origin='lower', cmap='jet',aspect='auto', interpolation='nearest', extent=(0, time_total_min, two_theta_start, two_theta_final))
    ax1.axvline(x=time_min[40], c='black', linestyle='dotted')
    plt.xlabel('Time / min', size=16) # Label X
    plt.ylabel('q / $\AA^{-1}$', size=16) # Label Y
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)

    box = ax1.get_position()
    pad, width = 0.02, 0.02
    cax = fig.add_axes([box.xmax + pad, box.ymin, width, box.height])
    fig.colorbar(cont, cax=cax)
    plt.yticks(fontsize=16)


    if save_png:
        plt.savefig(plotDir + sample_name + '_dif_map.png',dpi=500)
    if save_eps:
        plt.savefig(plotDir + sample_name + '_dif_map.pdf')

def plot_map(time_min, sample_name, diffra_map, two_theta_start, two_theta_final, plotDir, save_eps=False,save_png=False):
    fig=plt.figure(figsize=(8,8))

    cont=plt.imshow(diffra_map, origin='lower', cmap='jet',aspect='auto', interpolation='nearest', extent=(1, time_total_min, two_theta_start, two_theta_final))
    plt.xlabel('Time / min', size=16) # Label X
    plt.ylabel('q / $\AA^{-1}$', size=16) # Label Y
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)


    if save_png:
        plt.savefig(plotDir + sample_name + '_dif_map.png',dpi=500)
    if save_eps:
        plt.savefig(plotDir + sample_name + '_dif_map.pdf')

def plot_two_xrd(first, second, two_theta, diffra_map):

    plt.figure(figsize=(8,4), dpi=80)
    plt.plot(two_theta,diffra_map[first],linewidth=0.5,linestyle='-',  color='C0', marker='o',markersize=3,markerfacecolor='C0', label=time1_leg)
    plt.plot(two_theta,diffra_map[second],linewidth=0.5,linestyle='-',  color='C1', marker='o',markersize=3,markerfacecolor='C1', label=time2_leg)
    #plt.plot(fhand2[0],diffra_map_2[y2],linewidth=0.5,linestyle='-',  color='b', marker='o',markersize=3,markerfacecolor='r', label=time2_leg)
    plt.xlim(.7, 1.2)
    #plt.ylim(0.5, 4)
    plt.xlabel('2Theta / Degrees', size=16) # Label X
    plt.ylabel('Intensity / a.u', size=16) # Label Y
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.title(sample_name)
    plt.legend(fontsize=16)


def generate_model(spec):
    composite_model = None
    params = None
    x = spec['x']
    y = spec['y']
    x_min = np.min(x)
    x_max = np.max(x)
    x_range = x_max - x_min
    y_max = np.max(y)
    for i, basis_func in enumerate(spec['model']):
        prefix = f'm{i}_'
        model = getattr(models, basis_func['type'])(prefix=prefix)
        if basis_func['type'] in ['GaussianModel', 'LorentzianModel', 'VoigtModel']: # for now VoigtModel has gamma constrained to sigma
            model.set_param_hint('sigma', min=1e-6, max=x_range)
            model.set_param_hint('center', min=x_min, max=x_max)
            model.set_param_hint('height', min=0, max=y_max)
            model.set_param_hint('amplitude', min=1e-6)
            # default guess is horrible!! do not use guess()
            default_params = {
                prefix+'center': x_min + x_range * random.random(),
                prefix+'height': y_max,
                prefix+'sigma': x_range * random.random()
            }
        else:
            raise NotImplemented(f'model {basis_func["type"]} not implemented yet')
        if 'help' in basis_func:  # allow override of settings in parameter
            for param, options in basis_func['help'].items():
                model.set_param_hint(param, **options)
        model_params = model.make_params(**default_params, **basis_func.get('params', {}))
        if params is None:
            params = model_params
        else:
            params.update(model_params)
        if composite_model is None:
            composite_model = model
        else:
            composite_model = composite_model + model
    return composite_model, params

def update_spec_from_peaks(spec, model_indicies, peak_widths=(10, 25), **kwargs):
    x = spec['x']
    y = spec['y']
    x_range = np.max(x) - np.min(x)
    peak_indicies = find_peaks_cwt(y, peak_widths)
    np.random.shuffle(peak_indicies)
    for peak_indicie, model_indicie in zip(peak_indicies.tolist(), model_indicies):
        model = spec['model'][model_indicie]
        if model['type'] in ['GaussianModel', 'LorentzianModel', 'VoigtModel']:
            params = {
                'height': y[peak_indicie],
                'sigma': x_range / len(x) * np.min(peak_widths),
                'center': x[peak_indicie]
            }
            if 'params' in model:
                model.update(params)
            else:
                model['params'] = params
        else:
            raise NotImplemented(f'model {basis_func["type"]} not implemented yet')
    return peak_indicies

def plot_find(x, y, peaks_found):
    fig, ax = plt.subplots()
    ax.scatter(x, y, s=4)
    for i in peaks_found:
        ax.axvline(x=x[i], c='black', linestyle='dotted')


def plot_decv(spec,output):
    fig, ax = plt.subplots()
    ax.scatter(spec['x'], spec['y'], s=4)
    components = output.eval_components(x=spec['x'])
    print(len(spec['model']))
    for i, model in enumerate(spec['model']):
       ax.plot(spec['x'], components[f'm{i}_'])

def read_col(files, column):
    x = np.loadtxt(files,skiprows=0,usecols=(0,1),unpack=True)
    xx = x[column]
    return xx

def tab_params(output):

    for i in output.fit_report().split('\n'):
        if i.find('m0_amplitude:')==4:
            results[sample_name]['m0_amplitude:'].append(float(i.split(' ')[6]))
        elif i.find('m0_sigma:')==4:
            results[sample_name]['m0_sigma:'].append(float(i.split(' ')[10]))
        elif i.find('m0_center:')==4:
            results[sample_name]['m0_center:'].append(float(i.split(' ')[9]))
        elif i.find('m0_fwhm:')==4:
            results[sample_name]['m0_fwhm:'].append(float(i.split(' ')[11]))
        elif i.find('m0_height:')==4:
            results[sample_name]['m0_height:'].append(float(i.split(' ')[9]))
    return results

def infoSave(tab_par, pathTab, sample_name):
    infos=[]
    cols3=['#m0_amplitude', 'm0_sigma', 'm0_center', 'm0_fwhm', 'm0_height']
    keys=[]
    for i in tab_par[sample_name].keys():
        keys.append(i)

    infos.append(tab_par[sample_name][keys[0]])
    infos.append(tab_par[sample_name][keys[1]])
    infos.append(tab_par[sample_name][keys[2]])
    infos.append(tab_par[sample_name][keys[3]])
    infos.append(tab_par[sample_name][keys[4]])
    tab_infos=pd.DataFrame(data=np.transpose(infos), columns=cols3) #cira um dataframe das informações sobre as medidas

    f3 = open(pathTab + sample_name + '_params' + '.txt', 'w')
    tab_infos.to_csv(f3, sep='\t')
    f3.close()
##########################################################################
##Read Data
#########################################################################
for c in range(0,len(listdir)):

    path = pat+listdir[c].split('/')[-1] + '/'
    sample_name = path.split('/')[-2][0:-5]

    print(path)

    pathTab=dirdata(path)
    plotDir=dirplot(path)

    files = files_list(path, '.dat')


    files2 = files_list(path, '.info') #infos

    files3 = files_list(path, '.csv')


    spin = np.loadtxt(files3[0],skiprows=0,usecols=(0,),unpack=True)

    time = read_infos(files2[0], info='time')
    time_min=[]
    for i in time:
        time_min.append(i/60)


    time_total_s=np.max(time)
    time_total_min=float(time_total_s)/60

    diffra_map = read_files(files)

    base_diffra_map = normBaseline(diffra_map)

    sub_diffra_map = subBaseline(base_diffra_map, diffra_map)



    diffra_map_2_transp=np.transpose(diffra_map)
    sub_diffra_map_2_transp=np.transpose(sub_diffra_map)

    q = np.loadtxt(files[-1],skiprows=0,usecols=(0,),unpack=True)

    wavelenght=0.15406
    q_start=q[0]
    q_final=q[-1]
    two_theta_start=math.degrees(10*2*np.arcsin((q_start*wavelenght)/(4*math.pi)))
    two_theta_final=math.degrees(10*2*np.arcsin((q_final*wavelenght)/(4*math.pi)))


    map_temp = plot_map_temp(time_min, spin, sample_name, diffra_map_2_transp, two_theta_start, two_theta_final, plotDir, save_eps=True)
    sub_map_temp = plot_map_temp(time_min, spin, sample_name + '_sub', sub_diffra_map_2_transp, two_theta_start, two_theta_final, plotDir, save_eps=True)
    q_map_temp = plot_map_temp(time_min, spin, sample_name + '_q', diffra_map_2_transp, q_start, q_final, plotDir, save_eps=True)
    sub_q_map_temp = plot_map_temp(time_min, spin, sample_name + '_q_sub', sub_diffra_map_2_transp, q_start, q_final, plotDir, save_eps=True)

    q_i=.9
    q_f=1.10

    N_q_i_peak_2=[]
    N_q_f_peak_2=[]

    step = (q[-1]-q[0])/len(q)
    N_q_i_peak_2 = int((q_i-q[0])/step)
    N_q_f_peak_2 = int((q_f-q[0])/step)

    part_map=[]
    for i in range(40,len(files)):

        part_map.append(diffra_map[i][N_q_i_peak_2:N_q_f_peak_2])


    part_map_temp = plot_map(time_min, sample_name + '_q_part', np.transpose(part_map),q_i, q_f, plotDir, save_eps=True)

    time1=0
    time2= int(time[-1]/60)
    time1_leg=str(time1) +' min.'
    time2_leg=str(time2) +' min.'



    y1=int(time1/(time_total_min/len(files)))

    y2=int(time2/(time_total_min/len(files)))-1



    wavelenght=0.15406
    two_theta = q2theta(q, wavelenght)



#    xrds=plot_two_xrd(y1, y2, q, sub_diffra_map)


    fhand = read_col(files[y2], 1)
    q = read_col(files[y2], 0)




    results={}
    results[sample_name]={}
    results[sample_name]['m0_amplitude:']=[]
    results[sample_name]['m0_sigma:']=[]
    results[sample_name]['m0_center:']=[]
    results[sample_name]['m0_fwhm:']=[]
    results[sample_name]['m0_height:']=[]


    q_initial=.97
    q_final=1.04

    N_q_initial_peak_2=[]
    N_q_final_peak_2=[]

    step = (q[-1]-q[0])/len(q)
    N_q_initial_peak_2 = int((q_initial-q[0])/step)
    N_q_final_peak_2 = int((q_final-q[0])/step)

    q1=q[N_q_initial_peak_2:N_q_final_peak_2]

    for i in range(0,len(files)):

        fhand1=sub_diffra_map[i][N_q_initial_peak_2:N_q_final_peak_2]

        spec = {
            'x': q1,
            'y': fhand1,
            'model': [
                {'type': 'GaussianModel'},

            ]
        }

        peaks_found = update_spec_from_peaks(spec, [0], peak_widths=(15,))


    #    plot_find(spec['x'], spec['y'], peaks_found)



        model, params = generate_model(spec)
        output = model.fit(spec['y'], params, x=spec['x'])

    #    fig = output.plot(data_kws={'markersize':  1})

    #    plot_decv(spec, output)

    #    out=output.fit_report().split('\n')

        tab_par=tab_params(output)

    #plt.plot(time[40:],tab_par[sample_name]['m0_amplitude:'][39:], linestyle=' ',marker='o',markersize=3)



    infoSave(tab_par, pathTab, sample_name)
