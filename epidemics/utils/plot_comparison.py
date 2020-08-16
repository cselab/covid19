import sys
sys.path.append('../../')
import matplotlib
matplotlib.use('Agg')
# from epidemics.cantons.data.canton_population import CANTON_LIST, CANTON_LIST_SHORT
import os
from subprocess import call
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt 
import json
import numpy as np
import seaborn as sns
import matplotlib.patches as mpatches

tags = {'australia':   'AU',
        'canada':      'CA',
        'china':       'CN',
        'france':      'FR',
        'germany':     'DE',
        'italy':       'IT',
        'japan':       'JP',
        'russia':      'RU',
        'south korea': 'KR',
        'spain':       'ES',
        'switzerland': 'CH',
        'uk':          'UK',
        'us':          'US'
        }

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


def set_axis_style(ax, labels):
    ax.get_xaxis().set_tick_params(direction='out')
    ax.xaxis.set_ticks_position('bottom')
    ax.set_xticks(np.arange(1, len(labels) + 1))
    ax.set_xticklabels(labels)
    ax.set_xlim(0.25, len(labels) + 0.75)
    # ax.set_xlabel('Sample name')

def get_prior(ref_data,variable):

    n_var = len(ref_data['Variables'])

    idx = [ref_data['Variables'][i]['Distribution Index'] for i in range(n_var) 
            if ref_data['Variables'][i]['Name'] == variable][0]

    distrib_type = ref_data['Distributions'][idx]['Type']

    if distrib_type == 'Univariate/Uniform':
        samples = np.random.uniform(low=ref_data['Distributions'][idx]['Minimum'],
                                    high=ref_data['Distributions'][idx]['Maximum'],
                                    size=10000)

    elif distrib_type == 'Univariate/Gaussian':
        samples = np.random.normal(loc=25,scale=10,size=10000)

    return samples

def plot_parameters_comparison(folder,models,countries,variable,save_dir):

    line_color = 'gray'
    face_colors = ['#d53e4f','#99d594','#3288bd']
    face_colors = ['#d53e4f','#1a9850','#3288bd']
    alpha = 0.7
    med_width_factor = 2

    plot_medians = True
    plot_centers = False
    plot_prior = True

    # Get data
    data_all = {}
    for model in models:

        print(model)
        countries_folders = [folder+'/'+country+'/'+model+ '/_korali_samples' for country in countries]
        n_countries = len(countries)

        ref_data = get_samples_data(countries_folders[0])
        variables = [ref_data['Variables'][i]['Name'] for i in range(len(ref_data['Variables']))]
        var_idx = variables.index(variable)
        prior_samples = get_prior(ref_data,variable)
        
        bin_max = 0
        bin_min = 10

        if plot_prior:
            data_to_plot = [prior_samples]
        else:
            data_to_plot = []
        for i in range(n_countries):
            country = countries[i]
            country_path = countries_folders[i]
            print('   ({}/{}) {}'.format(i+1,n_countries,country))
            data = get_samples_data(country_path)
            numdim = len(data['Variables'])
            samples = data['Solver']['Posterior Sample Database']
            numentries = len(samples)
            samplesTmp = np.reshape(samples,(numentries,numdim))

            data_to_plot.append(samplesTmp[:,var_idx])

        data_all[model] = data_to_plot

    print('Plotting {}'.format(variable))

    # Plotting 

    # Labels
    common = os.path.commonprefix(models)
    unique = [model.replace('country.reparam.','') for model in models]

    fig, ax = plt.subplots(nrows = 1, ncols = 1,figsize =(18, 9))
    ax.grid(which='minor', axis='y', linestyle='--')
    # Matplotlib
    # def add_label(violin, label):
    #     color = violin["bodies"][0].get_facecolor().flatten()
    #     labels.append((mpatches.Patch(color=color), label))
    labels = []
    for i, model in enumerate(models):
        violins = plt.violinplot(data_all[model],vert=True,showmedians=True)

        # Make all the violin statistics marks a specific color:
        for partname in ('cbars','cmins','cmaxes'):
            vp = violins[partname]
            vp.set_edgecolor(line_color)
            vp.set_linewidth(1)

        # Plot medians
        if plot_medians or plot_centers:
            vp = violins['cmedians']
            segments = vp.get_segments()

            if plot_medians:
                med_width = segments[0][1][0]-segments[0][0][0]
                med_width *= med_width_factor
                for j in range(len(segments)):
                    segments[j][0][0] -=0.5*med_width
                    segments[j][1][0] +=0.5*med_width
                vp.set_segments(segments)
                vp.set_edgecolor(face_colors[i])
                vp.set_linewidth(1)
            else:
                vp.set_alpha(0)

            if plot_centers:
                centers = [[(segment[0][0]+segment[1][0])/2,segment[0][1]] for segment in segments]
                plt.plot([c[0] for c in centers],[c[1] for c in centers],color=face_colors[i])
                

        for vp in violins['bodies']:
            vp.set_facecolor(face_colors[i])
            vp.set_edgecolor(line_color)
            vp.set_linewidth(1)
            vp.set_alpha(alpha)

        labels.append((mpatches.Patch(color=face_colors[i],alpha=alpha), unique[i]))
    
    plt.legend(*zip(*labels), loc=2)
    plt.title(common[:-1])
    
    x_labels = [tags[country] for country in countries]
    if plot_prior:
        x_labels = ['Prior'] + x_labels

    set_axis_style(ax, x_labels)
    plt.ylabel(variable)


    create_folder(save_dir+'/_figures/')

    output = save_dir+'/_figures/'+variable+'_'+'-'.join(unique)+'.pdf'
    print("Creating output {}".format(output))
    plt.savefig(output)

    # plt.savefig(save_dir+'/'+variable+'_phase_'+str(phase)+'.pdf')
    # #Seaborn101
    # colors = ['green','blue','orange']
    # for i, model in enumerate(models):
    #     sns.violinplot(data=[d for d in data_all[model]],color=colors[i],inner=None,saturation=0.7)
    #     for violin, alpha in zip(ax.collections[::2], [0.8,0.6,0.4,0.2]):
    #         violin.set_alpha(alpha)

if __name__ == "__main__":  

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--folder', '-df', default='./data', help='Main results folder')
    parser.add_argument('--models', '-m', default='country.reparam.sir_int.nbin', type=str, nargs='+', help='Model type')
    parser.add_argument('--variables', '-v', default='R0', type=str, nargs='+', help='Model type')
    parser.add_argument('--countries', '-c', default=['canada'], type=str, nargs='+', help='Model type')
    parser.add_argument('--save_dir', '-sd', default='./', help='Model type')

    args = parser.parse_args()

    args.models = ['country.reparam.sir_int.nbin']
    args.variables = ['R0']
    args.folder = '/scratch/wadaniel/covid19/intervention/data/g9'

    countries = ['canada','china','france','germany','italy',
                 'japan','russia', 'uk','us']
    args.countries = countries #['canada', 'france', 'germany', 'italy']

    for variable in args.variables:
        plot_parameters_comparison(args.folder,args.models,args.countries,variable,args.save_dir)