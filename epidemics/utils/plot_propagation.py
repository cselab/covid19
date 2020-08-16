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

def get_stat(data, nt, pct):
    mean   = np.zeros((nt,1))
    median = np.zeros((nt,1))
    quant  = np.zeros((2,nt))
    for k in range(nt):
        median[k]  = np.quantile( data[:,k],0.5)
        quant[0,k] = np.quantile( data[:,k],0.5+pct/2)
        quant[1,k] = np.quantile( data[:,k],0.5-pct/2)
        mean[k]    = np.mean( data[:,k] )

    return mean, median, quant



def plot_samples_data(paths, output, pct=0.90, ndraws=100):
    
    fig = plt.figure(figsize=(12, 8))
    #fig.suptitle(self.modelDescription + '  (' + self.country + ')')

    ax  = fig.subplots(2,2)

    ax_daily_incidence = ax[0][0]
    ax_cumul_incidence = ax[1][0]

    ax_daily_deaths = ax[0][1]
    ax_cumul_deaths = ax[1][1]


    for path in paths :
        print("Reading {} ..".format(path))
        with open(path) as f: 
            genJs = json.load(f)
            results = genJs["Samples"]
            
            ns = len(results)
            nt = results[0]["Saved Results"]["Length of Variables"]

            all_incidence = np.zeros((ns*ndraws,nt))
            all_deaths    = np.zeros((ns*ndraws,nt))

            for k in range(ns):
                nv = results[k]["Saved Results"]["Number of Variables"]
                dispI = np.array(results[k]["Saved Results"]["Dispersion Daily Deaths"])
                dispD = np.array(results[k]["Saved Results"]["Dispersion Daily Incidence"])
                
                death     = None
                incidence = None
                for variable in results[k]["Saved Results"]["Variables"]:
                    if variable["Name"] == "Daily Incidence":
                        incidence = np.array(variable["Values"])
                    if variable["Name"] == "Daily Deaths":
                        death     = np.array(variable["Values"])

                pi = incidence/(incidence+dispI)
                pd = death/(death+dispD)
                xi = np.asarray([ np.random.negative_binomial(dispI,1-pi) for _ in range(ndraws) ])
                xd = np.asarray([ np.random.negative_binomial(dispD,1-pd) for _ in range(ndraws) ])
                
                all_incidence[k*ndraws:(k+1)*ndraws,:] = xi
                all_deaths[k*ndraws:(k+1)*ndraws,:]    = xd
  
            all_incidence_cum = np.cumsum(all_incidence, axis=1)
            all_deaths_cum    = np.cumsum(all_deaths, axis=1)

            meani, mediani, quanti = get_stat(all_incidence, nt, pct)
            meanic, medianic, quantic = get_stat(all_incidence_cum, nt, pct)
            meand, mediand, quantd = get_stat(all_deaths, nt, pct)
            meandc, mediandc, quantdc = get_stat(all_deaths_cum, nt, pct)


            ax_daily_incidence.plot( range(nt), meani, '--', lw=2, label='Mean', color='black')
            ax_daily_incidence.plot( range(nt), mediani, '-', lw=2, label='Median', color='black')
            ax_daily_incidence.fill_between( range(nt), quanti[0,:], quanti[1,:],  alpha=0.3)

            ax_cumul_incidence.plot( range(nt), meanic, '--', lw=2, label='Mean', color='black')
            ax_cumul_incidence.plot( range(nt), medianic, '-', lw=2, label='Median', color='black')
            ax_cumul_incidence.fill_between( range(nt), quantic[0,:], quantic[1,:],  alpha=0.3)
 
            #ax_daily_incidence.plot( 0:nt, quanti[0,:], '-', lw=2, label='Median', color='black')
            #ax_daily_incidence.plot( 0:nt, quanti[0,:], '-', lw=2, label='Median', color='black')
            #ax.plot( self.data['Propagation']['x-data'], median, '--', lw=2, label='Median', color='black')

            ax_daily_deaths.plot( range(nt), meand, '--', lw=2, label='Mean', color='black')
            ax_daily_deaths.plot( range(nt), mediand, '-', lw=2, label='Median', color='black')
            ax_daily_deaths.fill_between( range(nt), quantd[0,:], quantd[1,:],  alpha=0.3)
  
            ax_cumul_deaths.plot( range(nt), meandc, '--', lw=2, label='Mean', color='black')
            ax_cumul_deaths.plot( range(nt), mediandc, '-', lw=2, label='Median', color='black')
            ax_cumul_deaths.fill_between( range(nt), quantdc[0,:], quantdc[1,:],  alpha=0.3)
 
            #ax.legend(loc='upper left')
            #ax.set_ylabel( ylabel )
            #x = range( np.ceil( max( self.data['Propagation']['x-data'] )+1 ).astype(int) )
            #ax.set_xticks( x[0:-1:14] )
    
    create_folder(output+'/_figures/')
    output = output+'/_figures/propagation.pdf'
    print("Creating output {}".format(output))
    plt.savefig(output)

def set_axis_style(ax, labels):
    ax.get_xaxis().set_tick_params(direction='out')
    ax.xaxis.set_ticks_position('bottom')
    ax.set_xticks(np.arange(1, len(labels) + 1))
    ax.set_xticklabels(labels)
    ax.set_xlim(0.25, len(labels) + 0.75)
    # ax.set_xlabel('Sample name')

def plot_propagation(folder,models,countries,save_dir):

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

    for country in countries:
        n_models = len(models)
        models_files = [folder+'/'+country+'/'+model+ '/_korali_propagation/latest' for model in models]

        ref_data = plot_samples_data(models_files, save_dir)

if __name__ == "__main__":  

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--folder', '-df', default='./data', help='Main results folder')
    parser.add_argument('--models', '-m', default='country.reparam.sir_int.nbin', type=str, nargs='+', help='Model type')
    parser.add_argument('--countries', '-c', default=['canada'], type=str, nargs='+', help='Model type')
    parser.add_argument('--save_dir', '-sd', default='./', help='Model type')

    args = parser.parse_args()
    plot_propagation(args.folder,args.models,args.countries, args.save_dir)
