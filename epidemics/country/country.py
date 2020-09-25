import os
import sys
import numpy as np
import datetime
from scipy.stats import truncnorm
from scipy.special import loggamma

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.ioff()
 
from epidemics.country.data.cases import CountryData
from epidemics.data.synthetic import SyntheticData
from epidemics.epidemics import EpidemicsBase
from epidemics.utils.misc import save_file, prepare_folder

class EpidemicsCountry( EpidemicsBase ):

  def __init__( self, **kwargs ):
    
    self.country           = kwargs.pop('country', 'switzerland')
    self.futureDays        = kwargs.pop('futureDays', 3)
    self.nPropagation      = kwargs.pop('nPropagation', 100)
    self.logPlot           = kwargs.pop('logPlot', True)
    self.nValidation       = kwargs.pop('nValidation', 0)
    self.percentages       = kwargs.pop('percentages', [0.5])
    self.plotMeanMedian    = kwargs.pop('plotMeanMedian', False)
    self.up_to_int         = kwargs.pop('up_to_int', False)
    self.useIntervention   = kwargs.pop('useIntervention', False)
    self.useInformedPriors = kwargs.pop('useInformedPriors', False)
    self.preprocess        = kwargs.pop('preprocess')

    self.defaults = { 
            'R0'    : (1.0, 25.0),
            'beta'  : (0.01, 30.0),
            'D'     : (1.0, 25.0),  # recovery period
            'F'     : (1.0, 50.0),  # removal period
            'gamma' : (0.01, 1.0),  # recovery rate
            'Z'     : (0.0, 25.0),  # latency period (latency == incubation period)
            'Zl'    : (0.0, 25.0),  # latency period (incubation period)
            'Y'     : (0.0, 10.0),  # preasymptomatic period
            'mu'    : (0.0, 1.0),   # reduction factor unreported
            'alpha' : (0.05, 1.0),   # reporting rate
            'eps'   : (0.01, 0.25), # case fatality rate
            'tact'  : (0.0, 100.0), # intervention time
            'dtact' : (0.0, 60.0),  # intervention duration
            'kbeta' : (0.0, 1.0),   # reduction factor
            'kexp'  : (0.1, 3.3),   # 90% decay in ~ (1,30) days
            'Sigma' : (0.0, 100.0), # sdev Normal
            'dof'   : (2.0, 100.0), # DoF StudentT
            'cdof'  : (0.0, 100.0), # multiplicator variance StudentT
            'r'     : (0.0, 100.0), # dispersion NB 
            'eps3'  : (0.0, 1.0),   # CZ model, death outside of ICU
            'eps4'  : (0.0, 1.0),   # CZ model, death rate in ICU
            'e0'    : (0.0, 10.0),  # multiplicator e0 init
            'p0'    : (0.0, 10.0),  # multiplicator p0 init
            'iu0'   : (0.0, 10.0),  # multiplicator iu0 init
            'delay' : (0.0, 14.0)   # delay for deaths
        }

    self.constants = {
            'gamma'  : 1.0/5.2,
            'D'      : 5.2, # from nature paper
            'D_sdev' : 2.8, # from nature paper
            'Z'      : 2.7,
            'dtact'  : 14.0,
            'eps'    : 0.01 # dummy
    }

    self.informed_priors = {
        'D_mean'    : 1.0/(np.log(2)/14.0), # median at 14 days, 87.5 pct at 6w (exponential dist.)
#        'D_shape'   : 4.337,  # median at 14 days, 99pct < 6w (gamma dist.)
#        'D_scale'   : 3.970,  # median at 14 days, 99pct < 6w
        'D_shape'   : 3.448979591836735,   # test (taken from Z, TODO: calibrate)
        'D_scale'   : 1.5076923076923074,  # test (taken from Z, TODO: calibrate)
        'Y_shape'   : 32.62105263157895,   # preasymptomatic period started from 2.3 days (0.8, 3.0 95%-CI)
        'Y_scale'   : 0.07050661503710874, # preasymptomatic period
        'Zl_shape'  : 1.1671052631578946,  # simulated and fitted latency period
        'Zl_scale'  : 2.4931393061857263,  # simulated and fitted latency period
        'Z_shape'   : 3.448979591836735,   # mean 5.2, sdev 2.8 (latency == incubation period)
        'Z_scale'   : 1.5076923076923074,  # mean 5.2, sdev 2.8
    }
  
    self.bz_constants = {
            'sigma' : 1.0/2.6,
            'gamma' : 1.0/2.6,
            'omega1': 1.0/5.0,
            'omega2': 1.0/6.0,
            'omega3': 1.0/10.0,
            'omega4': 1.0/11.2,
            'omega5': 1.0/10.5,
            'eps1'  : 0.035,
            'eps2'  : 0.3,
            'eps3'  : 0.35,
            'eps4'  : 0.23,
    }
    
    super().__init__( **kwargs )
  
    if(self.synthetic):
        self.regionalData = SyntheticData( self.datafile, self.useInfections, self.useDeaths )
    else:
        self.regionalData = CountryData(self.country, lastDay=self.lastDay,
                                        preprocess=self.preprocess, up_to_int=self.up_to_int)

    self.process_data()
 
  def save_data_path( self ):
      return ( self.dataFolder, self.country, self.modelName )
 
  def process_data( self ):
  
    N        = self.regionalData.populationSize
    infected = self.regionalData.infected
    deaths   = self.regionalData.deaths
    
    nt = len(self.regionalData.time)
    t  = self.regionalData.time[1:nt-self.nValidation]

    I0 = infected[0]
    S0 = N - I0
    y0 = S0, I0
 
    self.data['Model']['Initial Condition'] = y0
    self.data['Model']['Population Size']   = N
    self.data['Model']['Intervention Day']  = self.regionalData.tact
       
    if self.useInfections:
        incidences = np.diff( infected )
        incidences, t_incidences = self.filter_daily_data('daily incidences', incidences, t)
    else:
        incidences, t_incidences = [], []

    if self.useDeaths:
        deaths = np.diff( deaths )
        deaths, t_deaths = self.filter_daily_data('daily deaths', deaths, t)
    else:
        deaths, t_deaths = [], []
    
    if self.useIntervention:
        self.intday  = self.data["Model"]["Intervention Day"]
    else:
        self.intday  = 0

    tx = list(set(np.concatenate([t_incidences, t_deaths])))
    tx.sort()

    self.data['Model']['x-data'] = tx
    self.data['Model']['y-data'] = np.concatenate([incidences, deaths])
    self.data['Model']['x-infected'] = t_incidences
    self.data['Model']['y-infected'] = incidences
    self.data['Model']['x-deaths']   = t_deaths
    self.data['Model']['y-deaths']   = deaths

    print('[Epidemics] Lengths incidences {} deaths {} total {}'.format(len(incidences), \
            len(deaths),len(np.concatenate([incidences,deaths]))), flush=True)

    T = np.ceil( t[-1-self.nValidation] + self.futureDays )
    self.data['Propagation']['x-data'] = np.linspace(0,T,int(T+1))
    save_file( self.data, self.saveInfo['inference data'], 'Data for Inference', 'pickle' )


  def computational_model( self, s ):

    p  = s['Parameters']
    t  = self.data['Model']['x-data']
    y0 = self.data['Model']['Initial Condition']
    N  = self.data['Model']['Population Size']
    
    T   = np.ceil(t[-1])
    tt  = np.linspace(0, t[-1], int(T+1))
    sol = self.solve_ode(y0=y0,T=t[-1], t_eval = tt, N=N, p=p)

    eps = 1e-12
    # get infected
    infected = np.diff(sol.y) 
    infected[np.isnan(infected-infected)] = eps
    infected[infected < eps] = eps
    
    # get deaths
    if hasattr(sol, 'd'):
        deaths = np.diff(sol.d)
        deaths[np.isnan(infected-infected)] = eps
        deaths[deaths < eps] = eps
    else:
        deaths = []

    # Transform gradients
    if(self.sampler == 'mTMCMC' or self.sampler=='HMC'):
        s["Gradient Mean"] = sol.gradMu
        if self.likelihoodModel == 'Normal' or self.likelihoodModel == 'Positive Normal':
            s["Gradient Standard Deviation"] = sol.gradSig
        elif self.likelihoodModel == 'Negative Binomial':
            s["Gradient Dispersion"] = sol.gradSig
        else:
            print('Gradients not implemented for other likelihood models', Flush=True)
            sys.exit(0)

    y = np.array([])
    if self.useInfections:
        y = np.concatenate( (y, self.getCases( infected, self.data['Model']['x-infected'])) )

    if self.useDeaths:
        y = np.concatenate( (y, self.getCases( deaths, self.data['Model']['x-deaths'])) )
 
    s['Reference Evaluations'] = list(y)
    
    if self.likelihoodModel == 'Normal':
        s['Standard Deviation'] = ( p[-1] * y ).tolist()
    elif self.likelihoodModel == 'Positive Normal':
        s['Standard Deviation'] = ( p[-1] * y ).tolist()
    elif self.likelihoodModel == 'Positive StudentT' and self.modelName.endswith('_alt'):
        var = 1.0+p[-1]*p[-1]*y*y+1e-9
        dof = 2*var/(var-1.0)
        s['Degrees Of Freedom'] = dof.tolist()
    elif self.likelihoodModel == 'Positive StudentT':
        s['Degrees Of Freedom'] = [p[-1]] * len(y)
    elif self.likelihoodModel == 'Negative Binomial':
        s['Dispersion'] = [p[-1]] * len(y)


  def computational_model_propagate( self, s ):

    p  = s['Parameters']
    t  = self.data['Propagation']['x-data']
    y0 = self.data['Model']['Initial Condition']
    N  = self.data['Model']['Population Size']

    tt = [t[0]-1] + t.tolist()
    sol = self.solve_ode(y0=y0,T=t[-1],t_eval=t.tolist(), N=N,p=p)
    
    _, ir0    = y0
    incidences = np.diff(sol.y)
    incidences = np.append(ir0, incidences)
     
    eps = 1e-32
    incidences[incidences < eps] = eps
 
    recovered  = np.array([])
    exposed    = np.array([])
    unreported = np.array([])
    deaths     = np.array([])

    if hasattr(sol, 'r'):
        recovered = np.diff(sol.r)
        recovered = np.append(0, recovered)
        eps = 1e-32
        recovered[recovered < eps] = eps
 
    if hasattr(sol, 'e'):
        exposed = np.diff(sol.e)
        exposed = np.append(0, exposed)
        exposed[exposed < eps] = eps

    if hasattr(sol, 'iu'):
        unreported = np.diff(sol.iu)
        unreported = np.append(0, unreported)
        unreported[unreported < eps] = eps

    if hasattr(sol, 'd'):
        deaths = np.diff(sol.d)
        deaths = np.append(0, deaths)
        deaths[deaths < eps] = eps

    k = 0
    js = {}
    js['Variables'] = []

    js['Variables'].append({})
    js['Variables'][k]['Name'] = 'Daily Incidence'
    js['Variables'][k]['Values'] = list(incidences)
    k += 1

    if recovered.size is not 0:
        js['Variables'].append({})
        js['Variables'][k]['Name'] = 'Daily Recovered'
        js['Variables'][k]['Values'] = list(recovered)
        k += 1

    if exposed.size is not 0:
       js['Variables'].append({})
       js['Variables'][k]['Name'] = 'Daily Exposed'
       js['Variables'][k]['Values'] = list(exposed)
       k += 1

    if unreported.size is not 0:
       js['Variables'].append({})
       js['Variables'][k]['Name'] = 'Daily Unreported'
       js['Variables'][k]['Values'] = list(unreported)
       k += 1

    if deaths.size is not 0:
        js['Variables'].append({})
        js['Variables'][k]['Name'] = 'Daily Deaths'
        js['Variables'][k]['Values'] = list(deaths)
        k += 1

    js['Number of Variables'] = len(js['Variables'])
    js['Length of Variables'] = len(t)

    if self.likelihoodModel == 'Normal':
        js['Standard Deviation Daily Incidence'] = ( p[-1] * incidences ).tolist()
        js['Standard Deviation Daily Deaths']    = ( p[-1] * deaths ).tolist()
    elif self.likelihoodModel == 'Positive Normal':
        js['Standard Deviation Daily Incidence'] = ( p[-1] * incidences ).tolist()
        js['Standard Deviation Daily Deaths']    = ( p[-1] * deaths ).tolist()
    elif self.likelihoodModel == 'Positive StudentT' and self.modelName.endswith('_alt'):
        varI = 1.0+p[-1]*p[-1]*incidences*incidences+1e-9
        dofI = 2*varI/(varI-1.0)
        varD = 1.0+p[-1]*p[-1]*deaths*deaths+1e-9
        dofD = 2*varD/(varD-1.0)
        js['Degrees Of Freedom Daily Incidence'] = dofI.tolist()
        js['Degrees Of Freedom Daily Deaths']    = dofD.tolist()
    elif self.likelihoodModel == 'StudentT':
        js['Degrees Of Freedom Daily Incidence'] = (len(incidences)) * [p[-1]]
        js['Degrees Of Freedom Daily Deaths']    = (len(deaths)) * [p[-1]]
    elif self.likelihoodModel == 'Positive StudentT':
        js['Degrees Of Freedom Daily Incidence'] = (len(incidences)) * [p[-1]]
        js['Degrees Of Freedom Daily Deaths']    = (len(deaths)) * [p[-1]]
    elif self.likelihoodModel == 'Negative Binomial':
        js['Dispersion Daily Incidence'] = (len(incidences)) * [p[-1]]
        js['Dispersion Daily Deaths']    = (len(deaths)) * [p[-1]]

    s['Saved Results'] = js


  
  def plot_intervals( self, ns=10):

    fig = self.new_figure()

    if self.useDeaths:
        ax  = fig.subplots(2,2)
        ax_daily = ax[0][0]
        ax_cumul = ax[1][0]
        
        ax_daily_deaths = ax[0][1]
        ax_cumul_deaths = ax[1][1]


        ax_daily_deaths.plot( self.data['Model']['x-deaths'], self.data['Model']['y-deaths'], 'x', \
                lw=2, label='Daily Deaths (data)', color='black')
        
        cumul_deaths = np.cumsum(self.data['Model']['y-deaths'])
        ax_cumul_deaths.set_xlabel('time in days')
        ax_cumul_deaths.plot( self.data['Model']['x-deaths'], cumul_deaths, 'x', \
                lw=2, label='Cumulative Deaths (data)', color='black')
 
        self.compute_plot_intervals( 'Daily Deaths', ns, ax_daily_deaths, 'Daily Deaths' )
        self.compute_plot_intervals( 'Daily Deaths', ns, ax_cumul_deaths, 'Cumulative number of deaths', cumulate=1)
 
        if 'Daily Deaths' in self.propagatedVariables:
            self.compute_mean_median( 'Daily Deaths', 'black', ns, ax_daily_deaths, 'Daily Deaths')
            self.compute_mean_median( 'Daily Deaths', 'black', ns, ax_cumul_deaths, 'Cumulative number of deaths', cumulate=1)

    else:
        ax  = fig.subplots( 2 )
        ax_daily = ax[0]
        ax_cumul = ax[1]
        ax_cumul.set_xlabel('time in days')


    ax_daily.plot( self.data['Model']['x-infected'], self.data['Model']['y-infected'], 'o', \
            lw=2, label='Daily Infected (data)', color='black')
    
    cumul_infected = np.cumsum(self.data['Model']['y-infected'])
    ax_cumul.set_xlabel('time in days')
    ax_cumul.plot( self.data['Model']['x-infected'], cumul_infected, 'o', \
            lw=2, label='Cumulative Infected (data)', color='black')


    if self.nValidation > 0:
      sys.exit()
      ax[0].plot( self.data['Validation']['x-data'], self.data['Validation']['y-data'], 'x', \
              lw=2, label='Daily Infected (validation data)', color='black')


    self.compute_plot_intervals( 'Daily Incidence', ns, ax_daily, 'Daily Incidence' )
    self.compute_plot_intervals( 'Daily Incidence', ns, ax_cumul, 'Cumulative number of infected', cumulate=1)
 
    if self.plotMeanMedian:

        self.compute_mean_median( 'Daily Incidence', 'blue', ns, ax_daily, 'Daily Incidence' )
        self.compute_mean_median( 'Daily Incidence', 'blue', ns, ax_cumul, 'Daily Incidence', cumulate=1 )

        if 'Daily Recovered' in self.propagatedVariables:
          self.compute_mean_median( 'Daily Recovered', 'green', ns, ax_daily, 'Daily Recovered')
          self.compute_mean_median( 'Daily Recovered', 'green', ns, ax_cumul, 'Cumulative number of recovered', cumulate=1)

        if 'Daily Exposed' in self.propagatedVariables:
          self.compute_mean_median( 'Daily Exposed', 'yellow', ns, ax_daily, 'Daily Exposed')
          self.compute_mean_median( 'Daily Exposed', 'yellow', ns, ax[1], 'Cumulative number of exposed', cumulate=1)

        if 'Daily Unreported' in self.propagatedVariables:
          self.compute_mean_median( 'Daily Unreported', 'orange', ns, ax_daily, 'Daily Unreported')
          self.compute_mean_median( 'Daily Unreported', 'orange', ns, ax_cumul, 'Cumulative number of unreported', cumulate=1)


    file = os.path.join(self.saveInfo['figures'],'prediction.png');
    prepare_folder( os.path.dirname(file) )
    fig.savefig(file)

    if (self.display): 
        plt.show()
    else: 
        print("[Epidemics] Cant show figure, '$DISPLAY' not set..")

    plt.close(fig)


    fig = self.new_figure()
    ax  = fig.subplots( 1 )

    found = False
    if self.parameters[0]['Name'] == 'R0':
      found = True
      ax.hist( self.parameters[0]['Values'], bins=40, density=1)
    elif (self.parameters[0]['Name'] == 'beta' and self.parameters[1]['Name'] == 'gamma') :
      found = True
      ax.hist( self.parameters[0]['Values']/self.parameters[1]['Values'], bins=40, density=1)
    else:
        print("[Epidemics] Cant plot R0, R0 nor (beta, gamma) not found in variables.")

    if (found == True):
        file = os.path.join(self.saveInfo['figures'],'R0.png');
        fig.savefig(file)
        if (self.display): 
            plt.show()
        else:
            print("[Epidemics] Cant show figure, '$DISPLAY' not set..")

    plt.close(fig)
 
 
  def getCases( self, data, tidx):
    """ helper to extract cases """
    cases = []
    for t in tidx:
        cases = cases + [data[t-1]]

    return cases


  def filter_daily_data(self,field_name,field,t):
    """ helper to remove and filter invalid data """

    if ((field < 0).any()):
        print("[Epidemics] Warning, removing negative values from {}!!!".format(field_name))
        valid = field >= 0
        field = field[valid]
        t = t[valid]
    
    if ((field > 1e32).any()):
        print("[Epidemics] Warning, removing extremely large (>1e32) values from {}!!!".format(field_name))
        valid = field <= 1e32
        field = field[valid]
        t = t[valid]

    return field, t
