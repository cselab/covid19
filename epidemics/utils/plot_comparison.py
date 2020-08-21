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
from scipy.special import gamma

################################################################################

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

 vdict = {
        'R0':       'Basic Reproduction Number',
        'D':        'Sumptomatic Infectious Period',
        'Y':        'Presymptomatic Infectious Period',
        'Z':        'Latency Period',
        'Zl':       'Latency Period (before Presymptomatic)',
        'alpha':    'Reporting Rate',
        'eps':      'Mortality Rate',
        'mu':       'Reduction Factor',
        'kbeta':    'Intervention Reduction Factor',
        'tact':     'Intervention Time',
        'dtact':    'Intervention Duration',
        }

# Plotting Options
line_color = 'gray'
face_colors = ['#d53e4f','#99d594','#3288bd']
face_colors = ['#d53e4f','#1a9850','#3288bd']
alpha = 0.7
med_width_factor = 2

plot_medians = True
plot_centers = False
plot_prior = True

################################################################################

# ---------------------------------------------------------------------------- #
#                                Helpers                                       #
# ---------------------------------------------------------------------------- #


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

# ---------------------------------------------------------------------------- #
#                               Get data                                       #
# ---------------------------------------------------------------------------- #

def get_prior_info(ref_data,variable):

    n_var = len(ref_data['Variables'])

    idx = [ref_data['Variables'][i]['Distribution Index'] for i in range(n_var) 
            if ref_data['Variables'][i]['Name'] == variable][0]

    distrib = ref_data['Distributions'][idx]
    distrib.pop('Range', None)
    distrib.pop('Random Seed', None)
    
    return distrib

def get_data(folder,models,countries,variable):

    # Get data
    samples = {}
    prior_info = {}
    for model in models:

        print(model)
        countries_folders = [folder+'/'+country+'/'+model+ '/_korali_samples' for country in countries]
        n_countries = len(countries)

        ref_data = get_samples_data(countries_folders[0])
        variables = [ref_data['Variables'][i]['Name'] for i in range(len(ref_data['Variables']))]
        var_idx = variables.index(variable)
        
        prior_info[model] = get_prior_info(ref_data,variable)

        data_to_plot = []

        for i in range(n_countries):
            country = countries[i]
            country_path = countries_folders[i]
            print('   ({}/{}) {}'.format(i+1,n_countries,country))
            
            country_data = get_samples_data(country_path)
            numdim = len(country_data['Variables'])
            country_samples = country_data['Solver']['Posterior Sample Database']
            numentries = len(country_samples)
            samplesTmp = np.reshape(country_samples,(numentries,numdim))

            data_to_plot.append(samplesTmp[:,var_idx])

        samples[model] = data_to_plot

    # Check that the priors for each model are the same
    val = list(prior_info.values())
    assert(all(x==val[0] for x in val))
    prior_info = val[0]

    return {'samples':samples,'prior_info':prior_info,'countries':countries,'variable':variable,'models':models}
        
def get_prior(data,scale=True):
    prior_info = data['prior_info']
    y_max = np.max([np.max([np.max(data['samples'][model][i]) for i in range(len(data['countries']))]) for model in data['models']])
    y_min = np.min([np.min([np.min(data['samples'][model][i]) for i in range(len(data['countries']))]) for model in data['models']])

    n = 1000
    if prior_info['Type'] ==  'Univariate/Uniform':
        y = np.linspace(prior_info['Minimum'],prior_info['Maximum'],n)  
        x = np.ones(n)*1/(prior_info['Maximum']-prior_info['Minimum'])

    elif prior_info['Type'] == 'Univariate/Gamma':
        k = prior_info['Shape']
        theta = prior_info['Scale']
        y = np.linspace(0,np.max([y_max,5*k]),n)
        x = 1/(gamma(k)*theta**k)*y**(k-1)*np.exp(-y/theta)

    if scale:
        x = x/x.max()*0.3

    return y, x

# ---------------------------------------------------------------------------- #
#                                  Plots                                       #
# ---------------------------------------------------------------------------- #

def plot_violin_style(data,save_dir):

    models = data['models']
    variable = data['variable']
    countries = data['countries']
    print('Plotting {}'.format(variable))

    # Labels
    common = os.path.commonprefix(models)
    unique = [model.replace('country.reparam.','') for model in models]

    fig, ax = plt.subplots(nrows = 1, ncols = 1,figsize =(18, 9))
    ax.grid(which='minor', axis='y', linestyle='--')

    labels = []
    for i, model in enumerate(models):
        
        labels.append((mpatches.Patch(color=face_colors[i],alpha=alpha), unique[i]))

        # Get posterior samples and plot violins
        samples = data['samples'][model]
        positions = np.arange(1, len(samples) + 1)
        if plot_prior:
            positions += 1
        violins = plt.violinplot(dataset=samples,positions=positions,vert=True,showmedians=True)

        # Make all the violin statistics marks a specific color:
        for partname in ('cbars','cmins','cmaxes'):
            vp = violins[partname]
            vp.set_edgecolor(line_color)
            vp.set_linewidth(1)

        for vp in violins['bodies']:
            vp.set_facecolor(face_colors[i])
            vp.set_edgecolor(line_color)
            vp.set_linewidth(1)
            vp.set_alpha(alpha)

        # Median or center options
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
    
    if plot_prior:        
        y_prior, x_prior = get_prior(data)
        plt.plot(x_prior+1,y_prior,color=line_color)
        plt.fill_betweenx(y_prior,x_prior+1,np.ones(len(x_prior)),color=line_color,alpha=alpha)

    plt.title(common[:-1])
    plt.legend(*zip(*labels), loc=2)
    plt.legend(*zip(*labels), loc=2)
    plt.title("Comparison {}".format(vdict[variable]))
    
    x_labels = [tags[country] for country in countries]
    if plot_prior:
        x_labels = ['Prior'] + x_labels
    set_axis_style(ax, x_labels)
    plt.ylabel(variable)

    create_folder(save_dir+'/_figures/')
    print("Creating output {}".format(save_dir+'/_figures/'+variable+'_'+'-'.join(unique)+'.pdf'))
    plt.savefig(save_dir+'/_figures/'+variable+'_'+'-'.join(unique)+'_violin.pdf')



def plot_ridge_style(data,save_dir):

    models = data['models']
    variable = data['variable']
    countries = data['countries']
    N = len(countries)+1

    # Legends
    common = os.path.commonprefix(models)
    unique = [model.replace('country.reparam.','') for model in models]

    print('Plotting {}'.format(variable))
    fig, ax = plt.subplots(nrows = N, ncols = 1,figsize =(9, 18),constrained_layout=False)

    labels = []
    l_idx = 0
    for i, model in enumerate(models):
        labels.append((mpatches.Patch(color=face_colors[i],alpha=alpha), unique[i]))

        # Get posterior samples and plot violins
        samples = data['samples'][model]
        for j in range(1,N):
            p = sns.kdeplot(data=samples[j-1], ax=ax[j], shade=True, color=face_colors[i],gridsize=1000,  bw=1, legend=False)
            x,y = p.get_lines()[l_idx].get_data()
            # ax[j].plot(x,y,c='k')
            get_y = lambda val : (y[np.where([x<val])[-1][-1]] + y[np.where([x>val])[-1][0]])/2

            median = np.median(samples[j-1])
            q10 = np.quantile(samples[j-1],0.1)
            q90 = np.quantile(samples[j-1],0.9)

            y_median = get_y(median)
            y_q10 = get_y(q10)
            y_q90 = get_y(q90)

            ax[j].plot([median,median],[0,y_median],color=face_colors[i],alpha=0.9)
            # ax[j].plot([q10,q10],[0,y_q10],color=face_colors[i],linestyle='--')
            # ax[j].plot([q5,q5],[0,y_q5],color=face_colors[i],linestyle='--')

            if j < (N - 1): 
                ax[j].set_xticks([])

            if i == 0 and j == 1:
                n_lines = len(p.get_lines())                    
        l_idx += n_lines

    x_prior, y_prior = get_prior(data,scale=False)
    ax[0].plot(x_prior,y_prior,color=line_color)
    ax[0].fill_between(x_prior,y_prior,color=line_color,alpha=alpha)
    ax[0].set_xticks([])

    for j in range(N):
        ax[j].set_yticks([])
        
        ax[j].spines['left'].set_visible(False)
        ax[j].spines['right'].set_visible(False)
        ax[j].spines['top'].set_visible(False)
        ax[j].spines['bottom'].set_edgecolor('#444444')
        ax[j].spines['bottom'].set_linewidth(1)
        ax[j].set_ylim([0,ax[j].get_ylim()[1]])

        if j == 0:
            ax[j].text(0.02, 0.05, 'Prior', fontsize=17, transform = ax[j].transAxes) 
        else:
            ax[j].text(0.02, 0.05, countries[j-1].capitalize(), fontsize=17,transform = ax[j].transAxes) 

    # plt.legend(*zip(*labels), loc=2)
    ax[0].set_title('Model comparison for {}'.format(vdict[variable]),fontsize=17,fontweight='bold',pad=10)
    plt.legend(*zip(*labels),loc='upper center', bbox_to_anchor=(0.5, -0.15),
          fancybox=False, shadow=False, ncol=len(models),frameon=False,fontsize='x-large')

    create_folder(save_dir+'/_figures/')
    print("Creating output {}".format(save_dir+'/_figures/'+variable+'_'+'-'.join(unique)+'.pdf'))
    plt.subplots_adjust(left=0.1, right=0.9, top=0.92, bottom=0.08)
    plt.savefig(save_dir+'/_figures/'+variable+'_'+'-'.join(unique)+'_ridge.pdf')

    # x_max = np.max([np.max([np.max(data['samples'][model][i]) for i in range(len(data['countries']))]) for model in data['models']])
    # x_min = np.min([np.min([np.min(data['samples'][model][i]) for i in range(len(data['countries']))]) for model in data['models']])
    # x_max_prior = np.max(x_prior)
    # x_min_prior = np.max(x_prior)

    # x_lim = [np.min([x_min,x_min_prior]),np.max([x_max_prior,x_max])]


    # for j in range(N):
    #     ax[j].set_xlim(x_lim)
    # plt.title(common[:-1])
    # plt.legend(*zip(*labels), loc=2)
    
    # x_labels = [tags[country] for country in countries]
    # if plot_prior:
    #     x_labels = ['Prior'] + x_labels
    # set_axis_style(ax, x_labels)
    # plt.ylabel(variable)

    # def ax_settings(ax):
    #     # ax.set_xlim(x_min,x_max)
    #     ax.set_yticks([])
        
    #     ax.spines['left'].set_visible(False)
    #     ax.spines['right'].set_visible(False)
    #     ax.spines['top'].set_visible(False)
        
    #     ax.spines['bottom'].set_edgecolor('#444444')
    #     ax.spines['bottom'].set_linewidth(1)
    #     ax.set_ylim([0,ax.get_ylim()[1]])

    #     # ax.text(0.02, 0.05, var_name, fontsize=17, fontweight="bold", transform = ax.transAxes) 
    #     return None

if __name__ == "__main__":  

    # python3 plot_comparison.py -df /scratch/wadaniel/covid19/intervention/data/g9/ -m country.reparam.sird_int.nbin country.reparam.seirud_int.nbin country.reparam.seird_int.geo -v R0 -c canada china france germany italy japan russia switzerland uk us 

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--folder', '-df', default='./data', help='Main results folder')
    parser.add_argument('--models', '-m', default='country.reparam.sir_int.nbin', type=str, nargs='+', help='Model type')
    parser.add_argument('--variables', '-v', default='R0', type=str, nargs='+', help='Model type')
    parser.add_argument('--countries', '-c', default=['canada'], type=str, nargs='+', help='Model type')
    parser.add_argument('--save_dir', '-sd', default='./', help='Model type')

    args = parser.parse_args()

    for variable in args.variables:

        data = get_data(args.folder,args.models,args.countries,variable)

        # plot_violin_style(data,args.save_dir)
        plot_ridge_style(data,args.save_dir)
