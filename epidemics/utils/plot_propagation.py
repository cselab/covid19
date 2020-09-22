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



def plot_samples_data(paths, models, samplespath, country, output, pct=0.90, ndraws=5):
    models = [model.replace('country.reparam.','') for model in models]
    models = [model.replace('_int','') for model in models]
  
    line_color = 'gray'
    face_colors = ['#d53e4f','#1a9850','#3288bd']
    if (len(paths) > len(face_colors)):
        print("Too many models, not enough colors defined in script. Exit..")
        sys.exit()

    alpha = 0.4
    med_width_factor = 2

    plot_mean    = False
    plot_medians = True

   
    fig = plt.figure(figsize=(12, 8))

    ax  = fig.subplots(2,3)

    ax_daily_incidence = ax[0][0]
    ax_cumul_incidence = ax[1][0]
 
    ax_daily_unreported = ax[0][1]
    ax_cumul_unreported = ax[1][1]

    ax_daily_deaths = ax[0][2]
    ax_cumul_deaths = ax[1][2]
 
    incidences = []
    incidences_cum = []
    deaths = []
    deaths_cum = []
 
    ax_daily_incidence.set_title('Daily reported infections')
    ax_cumul_incidence.set_title('Total reported infections')
  
    ax_daily_unreported.set_title('Daily unreported infections')
    ax_cumul_unreported.set_title('Total unreported infections')
 
    ax_daily_deaths.set_title('Daily deaths')
    ax_cumul_deaths.set_title('Total deaths')

    with open(samplespath) as f: 
        genJs   = json.load(f)
        refdata = np.array(genJs["Problem"]["Reference Data"])
        ndata = len(refdata)
        mid   = int(ndata/2)
        incidences = refdata[:mid] 
        incidences_cum = np.cumsum(incidences)
        deaths = refdata[mid:] 
        deaths_cum = np.cumsum(deaths)

        ax_daily_incidence.plot( range(len(incidences)), incidences, 'o', markersize=2, color='black')
        ax_cumul_incidence.plot( range(len(incidences_cum)), incidences_cum, 'o', markersize=2, color='black')
        ax_daily_deaths.plot( range(len(deaths)), deaths, 'o', markersize=2, color='black')
        ax_cumul_deaths.plot( range(len(deaths_cum)), deaths_cum, 'o', markersize=2, color='black')
            


    for idx, path in enumerate(paths) :
        print("Reading {} ..".format(path))
        with open(path) as f: 
            genJs = json.load(f)
            results = genJs["Samples"]
            
            ns = len(results)
            nt = results[0]["Saved Results"]["Length of Variables"]

            all_incidence  = np.zeros((ns*ndraws,nt))
            all_unreported = np.zeros((ns*ndraws,nt))
            all_deaths     = np.zeros((ns*ndraws,nt))

            plotUnreported = False
            for k in range(ns):
                nv = results[k]["Saved Results"]["Number of Variables"]
                dispI = np.array(results[k]["Saved Results"]["Dispersion Daily Incidence"])
                dispD = np.array(results[k]["Saved Results"]["Dispersion Daily Deaths"])
                
                death      = None
                incidence  = None
                unreported = None
                for variable in results[k]["Saved Results"]["Variables"]:
                    if variable["Name"] == "Daily Incidence":
                        incidence = np.array(variable["Values"])
                    if variable["Name"] == "Daily Deaths":
                        death     = np.array(variable["Values"])
                    if variable["Name"] == "Daily Unreported":
                        plotUnreported = True
                        unreported = np.array(variable["Values"])

                pi = incidence/(incidence+dispI)
                pd = death/(death+dispD)
                xi = np.asarray([ np.random.negative_binomial(dispI,1-pi) for _ in range(ndraws) ])
                xd = np.asarray([ np.random.negative_binomial(dispD,1-pd) for _ in range(ndraws) ])
                
                all_incidence[k*ndraws:(k+1)*ndraws,:]  = xi
                all_deaths[k*ndraws:(k+1)*ndraws,:]     = xd
  
                if plotUnreported:
                    pu = unreported/(unreported+dispI)
                    xu = np.asarray([ np.random.negative_binomial(dispI,1-pu) for _ in range(ndraws) ])
                    all_unreported[k*ndraws:(k+1)*ndraws,:] = xu

            all_incidence_cum  = np.cumsum(all_incidence, axis=1)
            all_deaths_cum     = np.cumsum(all_deaths, axis=1)

            meani, mediani, quanti    = get_stat(all_incidence, nt, pct)
            meanic, medianic, quantic = get_stat(all_incidence_cum, nt, pct)
            meand, mediand, quantd    = get_stat(all_deaths, nt, pct)
            meandc, mediandc, quantdc = get_stat(all_deaths_cum, nt, pct)

            ax_daily_incidence.fill_between( range(nt), quanti[0,:], quanti[1,:],  alpha=alpha, color=face_colors[idx], label=models[idx])
            ax_cumul_incidence.fill_between( range(nt), quantic[0,:], quantic[1,:],  alpha=alpha, color=face_colors[idx])
            
            ax_daily_deaths.fill_between( range(nt), quantd[0,:], quantd[1,:],  alpha=alpha, color=face_colors[idx])
            ax_cumul_deaths.fill_between( range(nt), quantdc[0,:], quantdc[1,:],  alpha=alpha, color=face_colors[idx])
  
            if plotUnreported:
                all_unreported_cum        = np.cumsum(all_unreported, axis=1)
                meanu, medianu, quantu    = get_stat(all_unreported, nt, pct)
                meanuc, medianuc, quantuc = get_stat(all_unreported_cum, nt, pct)
               
                ax_daily_unreported.fill_between( range(nt), quantu[0,:], quantu[1,:],  alpha=alpha, color=face_colors[idx], label=models[idx])
                ax_cumul_unreported.fill_between( range(nt), quantuc[0,:], quantuc[1,:],  alpha=alpha, color=face_colors[idx])

            if plot_medians:
                ax_daily_incidence.plot( range(nt), mediani, '-', lw=1, color=line_color)
                ax_cumul_incidence.plot( range(nt), medianic, '-', lw=1, color=line_color)
                
                ax_daily_deaths.plot( range(nt), mediand, '-', lw=1, color=line_color)
                ax_cumul_deaths.plot( range(nt), mediandc, '-', lw=1, color=line_color)

                if plotUnreported:
                    ax_daily_unreported.plot( range(nt), medianu, '-', lw=1, color=line_color)
                    ax_cumul_unreported.plot( range(nt), medianuc, '-', lw=1, color=line_color)
  
            if plot_mean:
                ax_daily_incidence.plot( range(nt), meani, '--', lw=1, color=line_color)
                ax_cumul_incidence.plot( range(nt), meanic, '--', lw=1, color=line_color)
 
                ax_daily_deaths.plot( range(nt), meand, '--', lw=1, color=line_color)
                ax_cumul_deaths.plot( range(nt), meandc, '--', lw=1, color=line_color)

                if plotUnreported:
                   ax_daily_unreported.plot( range(nt), meanu, '--', lw=1, color=line_color)
                   ax_cumul_unreported.plot( range(nt), meanuc, '--', lw=1, color=line_color)

            _, imax = ax_daily_incidence.get_ylim()
            _, umax = ax_daily_unreported.get_ylim()
            iumax = max(imax, umax)
 
            _, icmax = ax_cumul_incidence.get_ylim()
            _, ucmax = ax_cumul_unreported.get_ylim()
            iucmax = max(icmax, ucmax)

            ax_daily_incidence.set_ylim(ymin=0.0, ymax=iumax)
            ax_daily_unreported.set_ylim(ymin=0.0, ymax=iumax)
            
            ax_cumul_incidence.set_ylim(ymin=0.0, ymax=iucmax)
            ax_cumul_unreported.set_ylim(ymin=0.0, ymax=iucmax)
 
            ax_daily_incidence.legend(loc="upper right")
            
    output = output+'/_figures/'
    create_folder(output)
    model_str = '-'.join(models)
    plot = output+"/propagation_{}_{}.pdf".format(country, model_str)
    print("Creating output {}".format(plot))
    plt.savefig(plot)

def plot_propagation(folder,models,countries,save_dir):

    # Get data
    data_all = {}

    for country in countries:
        n_models = len(models)
        models_files = [folder+'/'+country+'/'+model+ '/_korali_propagation/latest' for model in models]
        samples_file = folder+'/'+country+'/'+models[0]+'/_korali_samples/latest'
        ref_data = plot_samples_data(models_files, models, samples_file, country, save_dir)

if __name__ == "__main__":  

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--folder', '-df', default='./data', help='Main results folder')
    parser.add_argument('--models', '-m', default='country.reparam.sird_int.nbin', type=str, nargs='+', help='Model type')
    parser.add_argument('--countries', '-c', default=['canada'], type=str, nargs='+', help='Model type')
    parser.add_argument('--save_dir', '-sd', default='./', help='Model type')

    args = parser.parse_args()
    plot_propagation(args.folder, args.models, args.countries, args.save_dir)
