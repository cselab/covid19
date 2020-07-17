import sys
sys.path.append('../../')
from epidemics.cantons.data.canton_population import CANTON_LIST, CANTON_LIST_SHORT
import os
from subprocess import call
import argparse
import matplotlib.pyplot as plt 
import json
import numpy as np

def create_folder(name):
    if not os.path.exists(name):
        os.makedirs(name)

def get_samples_data(path):

    configFile = path + '/gen00000000.json'

    with open(configFile) as f: js = json.load(f)
    configRunId = js['Run ID']

    resultFiles = [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f)) and f.startswith('gen')]
    resultFiles = sorted(resultFiles)

    genList = { } 

    for file in resultFiles:
        with open(path + '/' + file) as f:
            genJs = json.load(f)
            solverRunId = genJs['Run ID']

        if (configRunId == solverRunId):
            curGen = genJs['Current Generation']
            genList[curGen] = genJs
    del genList[0]

    lastGen = 0
    for i in genList: 
     if genList[i]['Current Generation'] > lastGen:
      lastGen = genList[i]['Current Generation']

    return genList[lastGen]

def plot_histogram(ax,ax_i,ax_j,theta,var_idx):
    dim = theta.shape[1]
    num_bins = 30
    
    ax_loc = ax[ax_i,ax_j]

    hist, bins, _ = ax_loc.hist(theta[:, var_idx], num_bins, density=True, color='lightgreen', ec='black')
    
    hist = hist / np.max(hist) * (ax_loc.get_xlim()[1] - ax_loc.get_xlim()[0])
    bottom = ax_loc.get_xlim()[0]
    widths = np.diff(bins)
    ax_loc.cla()
    ax_loc.bar(bins[:-1], hist, widths, color='lightgreen', ec='black', bottom=bottom)
    ax_loc.set_ylim(ax_loc.get_xlim())
    ax_loc.set_yticklabels([])
    ax_loc.tick_params(axis='both', which='both', length=0)

    return np.max(bins), np.min(bins)

def plot_comparison_single_model(phase,model,variable,save_dir):
    print('Phase {}'.format(phase))
    if phase == '1':
        path = 'data/'+model+'/phase_1_results/'
        file_path = '/_korali_samples'
    elif phase == '3':
        path = 'data/'+model+'/phase_3_results/'
        file_path = '/_korali_samples'
        file_path = '/_korali_samples_phase_3'

    cantons = [folder for folder in os.listdir(path)]
    cantons_folders = [path+canton+'/'+model+ file_path for canton in cantons]
    n_cantons = len(cantons)

    ref_data = get_samples_data(cantons_folders[0])
    variables = [ref_data['Variables'][i]['Name'] for i in range(len(ref_data['Variables']))]
    var_idx = variables.index(variable)

    nrows = 6
    ncols = 5
    fig, ax = plt.subplots(nrows=nrows,ncols=ncols,figsize=(16,12))

    bin_max = 0
    bin_min = 10

    for i in range(n_cantons):
        canton = cantons[i]
        canton_path = cantons_folders[i]
        print('({}/{}) {}'.format(i+1,n_cantons,canton))
        data = get_samples_data(canton_path)
        numdim = len(data['Variables'])
        samples = data['Solver']['Sample Database']
        numentries = len(samples)
        samplesTmp = np.reshape(samples,(numentries,numdim))

        # Get subplot index
        ax_i = i//ncols
        ax_j = i%ncols

        bin_max_i, bin_min_i = plot_histogram(ax,ax_i,ax_j,samplesTmp,var_idx)
        bin_max = np.max([bin_max,bin_max_i])
        bin_min = np.min([bin_min,bin_min_i])
        ax[ax_i,ax_j].set_title(canton)

    for i in range(n_cantons):
        ax_i = i//ncols
        ax_j = i%ncols   
        ax[ax_i,ax_j].set_xlim([bin_min,bin_max])
        ax[ax_i,ax_j].set_yticklabels([])

    ax[5][1].set_visible(False)
    ax[5][2].set_visible(False)
    ax[5][3].set_visible(False)
    ax[5][4].set_visible(False)

    plt.savefig(save_dir+'/'+variable+'_phase_'+str(phase)+'.pdf')

if __name__ == "__main__":  

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--phase', '-ph', default='1', help='Which phases to plor.')
    parser.add_argument('--model', '-m', default='sir_int.nbin', help='Model type')
    parser.add_argument('--variable', '-v', default='R0', help='Model type')
    parser.add_argument('--save_dir', '-sd', default='data/sir_int.nbin/figures/_comparison/', help='Model type')

    args = parser.parse_args()

    if args.phase == 'all':
        plot_comparison_single_model('1',args.model,args.variable,args.save_dir)
        plot_comparison_single_model('3',args.model,args.variable,args.save_dir)

    else:
        plot_comparison_single_model(args.phase,args.model,args.variable,args.save_dir)



