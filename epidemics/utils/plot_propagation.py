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

vdict = {   'R0':       'Basic Reproduction Number (R0)',
            'D':        'Symptomatic Infectious Period (D)',
            'Y':        'Presymptomatic Infectious Period (Y)',
            'Z':        'Latency Period (Z)',
            'Zl':       'Latency Period \n (before Presymptomatic) (Z)',
            'alpha':    r'Reporting Rate ($\alpha$)',
            'eps':      'Mortality Rate (f)',
            'mu':       'Reduction Factor ($\mu$)',
            'kbeta':    'Intervention Reduction Factor',
            'tact':     'Intervention Time',
            'dtact':    'Intervention Duration',
            'delay':    'Isolation Period',
            'R0_after': 'Basic Reproduction Number \n After Intervention',
            'Re':       'Effective Reproduction Number (Re)'

        }

model_names = {
        'country.reparam.sirdelay_int.nbin':'SIRD',
        'country.reparam.seirdelay_int.nbin':'SEIRD',
        'country.reparam.seiirdelay_int.nbin':'SEIIRD',
        'country.reparam.saphiredelay_int.nbin':'SAPHIRE',
        'country.reparam.seirudelay_int.nbin':'SEIRUD',
        }

# Plotting Options
line_color = 'gray'
lw_daily = 3
lw_cumul = 3
face_colors = { 'SIRD': '#66c2a5',
                'SEIRD':'#fc8d62' ,
                'SEIIRD': '#8da0cb',
                'SAPHIRE':'#e78ac3',
                'SEIRUD': '#a6d854'
              }
alpha = 0.4
fontsize=14

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

def setup_axis(ax):
    ax.spines['left'].set_visible(True)
    ax.spines['left'].set_edgecolor('#444444')

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_edgecolor('#444444')
    ax.spines['bottom'].set_linewidth(1)


def plot_samples_data(paths, models, samplespath, country, output, pct=0.90, ndraws=5):
    # models = [model.replace('country.reparam.','') for model in models]
    # models = [model.replace('_int','') for model in models]
  

    med_width_factor = 2

    plot_mean    = False
    plot_medians = True


    fig = plt.figure(figsize=(8, 12))

    ax  = fig.subplots(3,2)

    ax_daily_incidence = ax[0][0]
    ax_cumul_incidence = ax[0][1]
 
    ax_daily_unreported = ax[1][0]
    ax_cumul_unreported = ax[1][1]

    ax_daily_deaths = ax[2][0]
    ax_cumul_deaths = ax[2][1]
 
    incidences = []
    incidences_cum = []
    deaths = []
    deaths_cum = []
 
    ax_daily_incidence.set_title('Daily reported infections',fontsize=fontsize)
    setup_axis(ax_daily_incidence)

    ax_cumul_incidence.set_title('Cumulative reported infections',fontsize=fontsize)
    setup_axis(ax_cumul_incidence)

    ax_daily_unreported.set_title('Daily unreported infections',fontsize=fontsize)
    setup_axis(ax_daily_unreported)

    ax_cumul_unreported.set_title('Cumulative unreported infections',fontsize=fontsize)
    setup_axis(ax_cumul_unreported)

    ax_daily_deaths.set_title('Daily deaths',fontsize=fontsize)
    setup_axis(ax_daily_deaths)

    ax_cumul_deaths.set_title('Cumulative deaths',fontsize=fontsize)
    setup_axis(ax_cumul_deaths)

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
            

    labels = []
    tags = [model_names[model] for model in models]
    for idx, path in enumerate(paths) :
        print("Reading {} ..".format(path))

        model = models[idx]
        model_tag = model_names[model]
        # labels.append(model_tag)
        labels.append((mpatches.Patch(color=face_colors[model_tag],alpha=alpha), model_tag))

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

            fill_medians = True
            if fill_medians:
                ax_daily_incidence.fill_between( range(nt), quanti[0,:], quanti[1,:],  alpha=alpha, facecolor=face_colors[model_tag],edgecolor=None, label=model_tag)
                ax_cumul_incidence.fill_between( range(nt), quantic[0,:], quantic[1,:],  alpha=alpha, facecolor=face_colors[model_tag],edgecolor=None)
                
                ax_daily_deaths.fill_between( range(nt), quantd[0,:], quantd[1,:],  alpha=alpha, facecolor=face_colors[model_tag],edgecolor=None)
                ax_cumul_deaths.fill_between( range(nt), quantdc[0,:], quantdc[1,:],  alpha=alpha, facecolor=face_colors[model_tag],edgecolor=None)
            else:
                ax_daily_incidence.plot( range(nt), quanti[1,:],  color=face_colors[model_tag], label=model_tag)
                ax_daily_incidence.plot( range(nt), quanti[1,:],  color=face_colors[model_tag], label=model_tag)

                ax_cumul_incidence.fill_between( range(nt), quantic[0,:], quantic[1,:],  alpha=alpha, color=face_colors[model_tag])
                
                ax_daily_deaths.fill_between( range(nt), quantd[0,:], quantd[1,:],  alpha=alpha, color=face_colors[model_tag])
                ax_cumul_deaths.fill_between( range(nt), quantdc[0,:], quantdc[1,:],  alpha=alpha, color=face_colors[model_tag])
            

            if plotUnreported:
                all_unreported_cum        = np.cumsum(all_unreported, axis=1)
                meanu, medianu, quantu    = get_stat(all_unreported, nt, pct)
                meanuc, medianuc, quantuc = get_stat(all_unreported_cum, nt, pct)
               
                ax_daily_unreported.fill_between( range(nt), quantu[0,:], quantu[1,:],  alpha=alpha, facecolor=face_colors[model_tag],edgecolor=None, label=model_tag)
                ax_cumul_unreported.fill_between( range(nt), quantuc[0,:], quantuc[1,:],  alpha=alpha, facecolor=face_colors[model_tag],edgecolor=None)

            if plot_medians:
                ax_daily_incidence.plot( range(nt), mediani, '-', lw=lw_daily, color=face_colors[model_tag])
                ax_cumul_incidence.plot( range(nt), medianic, '-', lw=lw_cumul, color=face_colors[model_tag])
                
                ax_daily_deaths.plot( range(nt), mediand, '-', lw=lw_daily, color=face_colors[model_tag])
                ax_cumul_deaths.plot( range(nt), mediandc, '-', lw=lw_cumul, color=face_colors[model_tag])

                if plotUnreported:
                    ax_daily_unreported.plot( range(nt), medianu, '-', lw=lw_daily, color=face_colors[model_tag])
                    ax_cumul_unreported.plot( range(nt), medianuc, '-', lw=lw_cumul, color=face_colors[model_tag])
  
            if plot_mean:
                ax_daily_incidence.plot( range(nt), meani, '--', lw=lw_daily, color=face_colors[model_tag])
                ax_cumul_incidence.plot( range(nt), meanic, '--', lw=lw_cumul, color=face_colors[model_tag])
 
                ax_daily_deaths.plot( range(nt), meand, '--', lw=lw_daily, color=face_colors[model_tag])
                ax_cumul_deaths.plot( range(nt), meandc, '--', lw=lw_cumul, color=face_colors[model_tag])

                if plotUnreported:
                   ax_daily_unreported.plot( range(nt), meanu, '--', lw=lw_daily, color=face_colors[model_tag])
                   ax_cumul_unreported.plot( range(nt), meanuc, '--', lw=lw_cumul, color=face_colors[model_tag])

            _, imax = ax_daily_incidence.get_ylim()
            _, umax = ax_daily_unreported.get_ylim()
            iumax = max(imax, umax)
 
            _, icmax = ax_cumul_incidence.get_ylim()
            _, ucmax = ax_cumul_unreported.get_ylim()
            iucmax = max(icmax, ucmax)

            # ax_daily_incidence.set_ylim(ymin=0.0, ymax=iumax)
            # #ax_daily_unreported.set_ylim(ymin=0.0, ymax=iumax)
            
            # ax_cumul_incidence.set_ylim(ymin=0.0, ymax=iucmax)
            # #ax_cumul_unreported.set_ylim(ymin=0.0, ymax=iucmax)
 
            # ax_daily_incidence.legend(loc="upper right")
    
    plt.legend(*zip(*labels),loc='upper center', bbox_to_anchor=(-0.1, -0.15),
      fancybox=False, shadow=False, ncol=3,frameon=False,fontsize='x-large')
    plt.subplots_adjust(left=0.1, right=0.9, top=0.92, bottom=0.1)
    plt.suptitle('{}'.format(country.capitalize()),fontsize=fontsize,fontweight='bold')

    output = output+'/_figures/'
    create_folder(output)
    model_str = '-'.join(tags)
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
