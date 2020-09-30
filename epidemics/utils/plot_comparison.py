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

vdict = {   'R0':       'Basic Reproduction Number ($R_0$)',
            'D':        'Symptomatic Infectious Period (D)',
            'Y':        'Presymptomatic Infectious Period (Y)',
            'Z':        'Latency Period (Z)',
            'Zl':       'Latency Period \n (before Presymptomatic) ($Z_P$)',
            'alpha':    'Reporting Rate ($\alpha$)',
            'eps':      'Mortality Rate (f)',
            'mu':       'Reduction Factor ($\mu$)',
            'kbeta':    'Intervention Reduction Factor',
            'tact':     'Intervention Time',
            'dtact':    'Intervention Duration',
            'delay':    'Isolation Period',
            'R0_after': 'Basic Reproduction Number \n After Intervention',
            'Re':       'Effective Reproduction Number ($R_e$)',
            'r':        'Dispersion'
        }

population = {
        'canada':      37057765,
        'china':       1392730000,
        'france':      66977107,
        'germany':     82905782,
        'italy':       60421760,
        'japan':       126529100,
        'russia':      144478050,
        'switzerland': 8513227,
        'uk':          66460344,
        'us':          326687501
        }

model_names = {
        'country.reparam.sirdelay_int.nbin':'SIRD',
        'country.reparam.seirdelay_int.nbin':'SEIRD',
        'country.reparam.seiirdelay_int.nbin':'SEIIRD',
        'country.reparam.saphiredelay_int.nbin':'SAPHIRED',
        'country.reparam.seirudelay_int.nbin':'SEIRUD',
        }

# Plotting Options
line_color = 'gray'
# face_colors = ['#d53e4f','#99d594','#3288bd']
# face_colors = ['#d53e4f','#1a9850','#3288bd']
#face_colors = ['#fbb4ae','#b3cde3','#ccebc5','#decbe4','#fed9a6']
face_colors = ['#66c2a5','#fc8d62','#8da0cb','#e78ac3','#a6d854']

face_colors = { 'SIRD': '#66c2a5',
                'SEIRD':'#fc8d62' ,
                'SEIIRD': '#8da0cb',
                'SAPHIRED':'#e78ac3',
                'SEIRUD': '#a6d854'
              }

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

def get_data(folder,models,countries,variable,get_prior_data=True):

    # Get data
    samples = {}
    prior_info = {}
    for model in models:

        print(model)
        countries_folders = [folder+'/'+country+'/'+model+ '/_korali_samples' for country in countries]
        n_countries = len(countries)

        ref_data = get_samples_data(countries_folders[0])
        variables = [ref_data['Variables'][i]['Name'] for i in range(len(ref_data['Variables']))]
        print('Variables: {}'.format(variables))

        if variable == 'R0_after':
            get_prior_data = False
        elif variable == 'Re':
            get_prior_data = False
        else:
            var_idx = variables.index(variable)
            if get_prior_data:
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

            if variable == 'R0_after':
                R0_idx = variables.index('R0')
                kbeta_idx = variables.index('kbeta')
                data_to_plot.append(np.multiply(samplesTmp[:,R0_idx],samplesTmp[:,kbeta_idx]))
            
            elif variable == 'Re':

                if model_names[model] == 'SIRD' or model_names[model] == 'SEIRD':
                    R0 = samplesTmp[:,variables.index('R0')]
                    data_to_plot.append(R0)
                elif model_names[model] == 'SEIIRD':
                    R0 = samplesTmp[:,variables.index('R0')]
                    mu = samplesTmp[:,variables.index('mu')]
                    alpha = samplesTmp[:,variables.index('alpha')]
                    data_to_plot.append(R0*alpha+R0*(1-alpha)*mu)

                elif model_names[model] == 'SAPHIRED':
                    R0 = samplesTmp[:,variables.index('R0')]
                    mu = samplesTmp[:,variables.index('mu')]
                    alpha = samplesTmp[:,variables.index('alpha')]
                    Y = samplesTmp[:,variables.index('Y')]
                    D = samplesTmp[:,variables.index('D')]

                    data_to_plot.append(mu*R0*Y/D+R0*(1-alpha)*mu+alpha*R0)
                elif model_names[model] == 'SEIRUD':
                    R0 = samplesTmp[:,variables.index('R0')]
                    Y = samplesTmp[:,variables.index('Y')]
                    D = samplesTmp[:,variables.index('D')]
                    alpha = samplesTmp[:,variables.index('alpha')]
                    data_to_plot.append(R0*Y/D+R0*(1-alpha))
            else:
                data_to_plot.append(samplesTmp[:,var_idx])
            

        samples[model] = data_to_plot

    # Check that the priors for each model are the same
    if get_prior_data:
        val = list(prior_info.values())
        assert(all(x==val[0] for x in val))
        prior_info = val[0]
    else:
        prior_info = {'Type':None}

    return {'samples':samples,'prior_info':prior_info,'countries':countries,'variable':variable,'models':models}
        
def get_prior(data,scale=True):
    prior_info = data['prior_info']
    y_max = np.max([np.max([np.max(data['samples'][model][i]) for i in range(len(data['countries']))]) for model in data['models']])
    y_min = np.min([np.min([np.min(data['samples'][model][i]) for i in range(len(data['countries']))]) for model in data['models']])

    n = 1000
    if prior_info['Type'] ==  'Univariate/Uniform':
        x = np.linspace(prior_info['Minimum'],prior_info['Maximum'],n)  
        y = np.ones(n)*1/(prior_info['Maximum']-prior_info['Minimum'])

        x_median = (prior_info['Maximum']+prior_info['Minimum'])/2
        y_median = get_median_y(x,y,x_median)
        median = [x_median,y_median]
    elif prior_info['Type'] == 'Univariate/Gamma':
        k = prior_info['Shape']
        theta = prior_info['Scale']
        x = np.linspace(0,np.max([y_max,5*k]),n)
        y = 1/(gamma(k)*theta**k)*x**(k-1)*np.exp(-x/theta)

        samples = np.random.gamma(k, scale=theta, size=10000)
        x_median = np.median(samples)
        y_median = get_median_y(x,y,x_median)
        median = [x_median,y_median]

    if scale:
        y = y/y.max()*0.3

    return x, y, median

def get_median_y(x,y,val):
    return (y[np.where([x<val])[-1][-1]] + y[np.where([x>val])[-1][0]])/2


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
        y_prior, x_prior, _ = get_prior(data)
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



def plot_ridge_style(data,save_dir,plot_prior=True):

    models = data['models']
    variable = data['variable']
    countries = data['countries']
    N = len(countries)+1

    # Legends
    common = os.path.commonprefix(models)
    unique = [model_names[model] for model in models]

    print('Plotting {}'.format(variable))
    fig, ax = plt.subplots(nrows = N, ncols = 1,figsize =(6, 18),constrained_layout=False)

    labels = []
    l_idx = 0
    for i, model in enumerate(models):
        labels.append((mpatches.Patch(color=face_colors[model_names[model]],alpha=alpha), unique[i]))

        # Get posterior samples and plot violins
        samples = data['samples'][model]
        for j in range(1,N):

            if data['prior_info']['Type'] == 'Univariate/Uniform':
                clip = [data['prior_info']['Minimum'],data['prior_info']['Maximum']]
            else:
                clip = [0,100]
            bw=0.1

            if model == 'country.reparam.saphired_int.nbin':
                if j == 2 and variable == 'eps':
                    bw = 0.1

            p = sns.kdeplot(data=samples[j-1], ax=ax[j], 
                clip = clip,
                shade=True, color=face_colors[model_names[model]],gridsize=1000,  bw=bw, legend=False)
            x,y = p.get_lines()[l_idx].get_data()

            median = np.median(samples[j-1])
            q10 = np.quantile(samples[j-1],0.1)
            q90 = np.quantile(samples[j-1],0.9)

            y_median = get_median_y(x,y,median)
            y_q10 = get_median_y(x,y,q10)
            y_q90 = get_median_y(x,y,q90)

            ax[j].plot([median,median],[0,y_median],color=face_colors[model_names[model]],alpha=0.9)
            # ax[j].plot([q10,q10],[0,y_q10],color=face_colors[i],linestyle='--')
            # ax[j].plot([q5,q5],[0,y_q5],color=face_colors[i],linestyle='--')

            if j < (N - 1): 
                ax[j].set_xticks([])

            if i == 0 and j == 1:
                n_lines = len(p.get_lines())                    
        l_idx += n_lines

    if plot_prior:
        x_prior, y_prior, prior_median = get_prior(data,scale=False)
        ax[0].plot(x_prior,y_prior,color=line_color)
        ax[0].fill_between(x_prior,y_prior,color=line_color,alpha=alpha)
        ax[0].set_xticks([])
        ax[0].plot([prior_median[0],prior_median[0]],[0,prior_median[1]],color=line_color,alpha=0.9)
        start_idx = 0
    else:
        start_idx = 1
        ax[0].set_xticks([])

    x_min = np.min([ax[i].get_xlim()[0] for i in range(start_idx,N)])
    x_max = np.max([ax[i].get_xlim()[1] for i in range(start_idx,N)])
    print(x_min,x_max)

    if data['prior_info']['Type'] == 'Univariate/Uniform':
        x_range = ax[0].get_xlim()
        x_min = x_range[0]
        x_max = x_range[1]

    for j in range(0,N):
        ax[j].set_yticks([])
        
        ax[j].spines['left'].set_visible(False)
        ax[j].spines['right'].set_visible(False)
        ax[j].spines['top'].set_visible(False)
        ax[j].spines['bottom'].set_edgecolor('#444444')
        ax[j].spines['bottom'].set_linewidth(1)
        ax[j].set_ylim([0,ax[j].get_ylim()[1]])

        ax[j].set_xlim([x_min,x_max])
        
        if variable == 'Y':
            ax[j].set_xlim([0,8])
        elif variable == 'R0_after':
            ax[j].set_xlim([0,3])
        # elif variable == 'Re':
        #     ax[j].set_xlim([0,20])

        if j == 0:
            if plot_prior:
                ax[j].text(0.02, 0.05, 'Prior', fontsize=17, transform = ax[j].transAxes) 
        else:
            ax[j].text(0.02, 0.05, tags[countries[j-1]], fontsize=17,transform = ax[j].transAxes) 

    # plt.legend(*zip(*labels), loc=2)
    ax[0].set_title('{}'.format(vdict[variable]),fontsize=17,fontweight='bold',pad=10)
    plt.legend(*zip(*labels),loc='upper center', bbox_to_anchor=(0.5, -0.15),
          fancybox=False, shadow=False, ncol=3,frameon=False,fontsize='x-large')
    plt.subplots_adjust(left=0.1, right=0.9, top=0.92, bottom=0.08)

    create_folder(save_dir+'/_figures/')
    print("Creating output {}".format(save_dir+'/_figures/'+variable+'_'+'-'.join(unique)+'.pdf'))
    plt.savefig(save_dir+'/posteriors/'+variable+'_'+'-'.join(unique)+'_ridge.pdf')


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

        if variable == 'R0_after':
            plot_prior = False
        elif variable == 'Re':
            plot_prior = False
        else:
            plot_prior = True
        data = get_data(args.folder,args.models,args.countries,variable,plot_prior)

        # plot_violin_style(data,args.save_dir)
        plot_ridge_style(data,args.save_dir,plot_prior)
