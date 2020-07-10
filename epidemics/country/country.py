import os
import numpy as np
from scipy.stats import truncnorm
from scipy.special import loggamma

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.ioff()
 
from epidemics.epidemics import EpidemicsBase
from epidemics.data.combined import RegionalData
from epidemics.data.synthetic import SyntheticData
from epidemics.tools.tools import save_file, prepare_folder

class EpidemicsCountry( EpidemicsBase ):

  def __init__( self, **kwargs ):
    
    self.country        = kwargs.pop('country', 'switzerland')
    self.futureDays     = kwargs.pop('futureDays', 2)
    self.nPropagation   = kwargs.pop('nPropagation', 100)
    self.logPlot        = kwargs.pop('logPlot', True)
    self.nValidation    = kwargs.pop('nValidation', 0)
    self.percentages    = kwargs.pop('percentages', [0.5, 0.95, 0.99])
    self.preprocess     = kwargs.pop('preprocess', False)
    self.plotMeanMedian = kwargs.pop('plotMeanMedian', False)
    self.up_to_int      = kwargs.pop('up_to_int', False)
    self.data_fields    = kwargs.pop('data_fields', 'infected')
    print(self.data_fields)


    self.defaults = { 
            'R0'    : (1.0, 15.0),
            'D'     : (1.0, 20.0),
            'Z'     : (1.0, 20.0),
            'mu'    : (0.0, 5.0),
            'alpha' : (0.0, 1.0),
            'eps'    : (0.0, 0.2),
            'tact'  : (0.0, 100.0),
            'dtact' : (0.0, 30.0),
            'kbeta' : (0.0, 1.0),
            'Sigma' : (0.0, 100.0),
            'r'     : (0.0, 100.0)
        }

    self.constants = {
            'gamma' : 1.0/5.2,
            'Z'     : 2.7,
            'dtact' : 14.0
    }
    
    super().__init__( **kwargs )
  
    if(self.synthetic):
        self.regionalData = SyntheticData( self.datafile )
    else:
        self.regionalData = RegionalData( self.country,self.preprocess,self.up_to_int)

    self.process_data()
 
  def save_data_path( self ):
      return ( self.dataFolder, self.country, self.modelName )
 
  def process_data( self ):

    t = []
    incidences = []
    t_incidences = []
    I0 = 0
    if 'infected' in self.data_fields:
        y = self.regionalData.infected
        t_incidences = self.regionalData.time[1:]
        incidences = np.diff( y[0:] )
        incidences, t_incidences = self.filter_daily_data('daily incidences',incidences,t_incidences)
        incidences = np.clip(incidences, a_min=0, a_max=1e32)
        I0 = y[0]
        t = t_incidences


    deaths = []
    t_deaths = []
    if 'deaths' in self.data_fields:
        y = self.regionalData.deaths
        t_deaths = self.regionalData.time[1:]
        deaths = np.diff( y[0:] )
        deaths, t_deaths = self.filter_daily_data('daily deaths',deaths,t_deaths)
        deaths = np.clip(deaths, a_min=0, a_max=1e32)
        if len(t) == 0:
            t = t_deaths


    N = self.regionalData.populationSize
    S0 = N - I0
    y0 = S0, I0

    if self.nValidation == 0:
        self.data['Model']['x-data'] = t
        self.data['Model']['y-data'] = np.concatenate([incidences,deaths])

        # For the SEIRD model
        self.data['Model']['x-eval'] = self.regionalData.time[1:]
        self.data['Model']['x-infected'] = t_incidences
        self.data['Model']['x-deaths'] = t_deaths
        self.data['Model']['y-infected'] = incidences
        self.data['Model']['y-deaths'] = deaths

        print('Lengths incidences {} deaths {} total {}'.format(len(incidences),
                len(deaths),len(np.concatenate([incidences,deaths]))))
        print(len(t_incidences),len(t_deaths))

    else:
        print('TODO, need to take into account diff between t_incidences and t_deaths')
        self.data['Model']['x-data'] = t[:-self.nValidation]
        self.data['Model']['y-data'] = incidences[:-self.nValidation]
        #self.data['Validation']['x-data'] = t[-self.nValidation:]
        #self.data['Validation']['y-data'] = incidences[-self.nValidation-1:]

    self.data['Model']['Initial Condition'] = y0
    self.data['Model']['Population Size']   = N

    if self.nValidation == 0:
        T = np.ceil( t[-1] + self.futureDays )
    else:
        T = np.ceil( t[-self.nValidation] + self.futureDays )
    
    self.data['Propagation']['x-data'] = np.linspace(0,T,int(T+1))
    save_file( self.data, self.saveInfo['inference data'], 'Data for Inference', 'pickle' )

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
 
    recovered  = None
    exposed    = None
    unreported = None
    deaths     = None

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

    eps = 1e-32
    
    k = 0
    js = {}
    js['Variables'] = []

    js['Variables'].append({})
    js['Variables'][k]['Name'] = 'Daily Incidence'
    js['Variables'][k]['Values'] = list(incidences)
    k += 1

    if recovered is not None:
        js['Variables'].append({})
        js['Variables'][k]['Name'] = 'Daily Recovered'
        js['Variables'][k]['Values'] = list(recovered)
        k += 1

    if exposed is not None:
       js['Variables'].append({})
       js['Variables'][k]['Name'] = 'Daily Exposed'
       js['Variables'][k]['Values'] = list(exposed)
       k += 1

    if unreported is not None:
       js['Variables'].append({})
       js['Variables'][k]['Name'] = 'Daily Unreported'
       js['Variables'][k]['Values'] = list(unreported)
       k += 1

    if deaths is not None:
        js['Variables'].append({})
        js['Variables'][k]['Name'] = 'Daily Deaths'
        js['Variables'][k]['Values'] = list(deaths)
        k += 1

    js['Number of Variables'] = len(js['Variables'])
    js['Length of Variables'] = len(t)

    if self.likelihoodModel == 'Normal':
        js['Standard Deviation'] = ( p[-1] * incidences ).tolist()
    elif self.likelihoodModel == 'Positive Normal':
        js['Standard Deviation'] = ( p[-1] * incidences ).tolist()
    elif self.likelihoodModel == 'Negative Binomial':
        js['Dispersion'] = (len(incidences)) * [p[-1]]

    s['Saved Results'] = js


  
  def plot_intervals( self, ns=10):

    fig = self.new_figure()

    if 'deaths' in self.data_fields:
        ax  = fig.subplots(2,2)
        ax_daily = ax[0][0]
        ax_daily_deaths = ax[0][1]
        ax_cumul = ax[1][0]
        ax_cumul_deaths = ax[1][1]
        ax_cumul.set_xlabel('time in days')
        ax_cumul_deaths.set_xlabel('time in days')


        ax_daily_deaths.plot( self.data['Model']['x-deaths'], self.data['Model']['y-deaths'], 'o', lw=2, label='Daily Deaths(data)', color='black')
        z = np.cumsum(self.data['Model']['y-deaths'])
        ax_cumul.plot( self.data['Model']['x-deaths'], z, 'o', lw=2, label='Cumulative Deaths(data)', color='black')

    else:
        ax  = fig.subplots( 2 )
        ax_daily = ax[0]
        ax_cumul = ax[1]
        ax_daily_deaths = ax_daily
        ax_cumul_deatsh = ax_cumul
        ax_cumul.set_xlabel('time in days')


    ax_daily.plot( self.data['Model']['x-infected'], self.data['Model']['y-infected'], 'o', lw=2, label='Daily Infected(data)', color='black')
    z = np.cumsum(self.data['Model']['y-infected'])
    ax_cumul.plot( self.data['Model']['x-infected'], z, 'o', lw=2, label='Cumulative Infected(data)', color='black')



    #if self.nValidation > 0:
    #  ax[0].plot( self.data['Validation']['x-data'], self.data['Validation']['y-data'], 'x', lw=2, label='Daily Infected (validation data)', color='black')


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

        if 'Daily Deaths' in self.propagatedVariables:
          self.compute_mean_median( 'Daily Deaths', 'black', ns, ax_daily_deaths, 'Daily Deaths')
          self.compute_mean_median( 'Daily Deaths', 'black', ns, ax_cumul_deatsh, 'Cumulative number of deaths', cumulate=1)

    file = os.path.join(self.saveInfo['figures'],'prediction.png');
    prepare_folder( os.path.dirname(file) )
    fig.savefig(file)

    if (self.display): 
        plt.show()
    else: 
        print("[Epidemics] Cant show figure, '$DISPLAY' not set..")

    plt.close(fig)

    #----------------------------------------------------------------------------------------------------------------------------------
    #----------------------------------------------------------------------------------------------------------------------------------

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


  def llk_model_nbin ( self, p, t, refy, y0, N ):

    tt = [t[0]-1] + t.tolist()
    sol = self.solve_ode(y0=y0,T=t[-1], t_eval = tt,N=N,p=p)

    # get incidences
    y = np.diff(sol.y)
 
    eps = 1e-32
    y[y < eps] = eps
    refy[refy < 0] = 0
  
    llk = 0.0
    for idx, incident in enumerate(y):
        llk -= loggamma ( refy[idx] + 1. ) #fix
        m    = incident
        r    = incident*p[-1]
        prob = m / (m+r)

        llk += loggamma ( refy[idx] + r )
        llk -= loggamma ( r )
        llk += r * np.log(1-prob)
        llk += refy[idx] * np.log(prob)

    return llk

  def llk_model_tnrm ( self, p, t, refy, y0, N ):

    tt = [t[0]-1] + t.tolist()
    sol = self.solve_ode(y0=y0,T=t[-1], t_eval = tt,N=N,p=p)

    # get incidences
    y = np.diff(sol.y)
 
    eps = 1e-32
    y[y < eps] = eps
    refy[refy < 0] = 0
  
    llk = 0.0
    for idx, incident in enumerate(y):
        std  = incident*p[-1]
        a    = -1.0*incident/std
        b    = np.inf
        llk += truncnorm.logpdf(refy[idx], a, b, incident, std)

    return llk


  def filter_daily_data(self,field_name,field,t):

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
