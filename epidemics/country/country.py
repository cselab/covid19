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
    
    self.country      = kwargs.pop('country', 'switzerland')
    self.futureDays   = kwargs.pop('futureDays', 2)
    self.nPropagation = kwargs.pop('nPropagation', 100)
    self.logPlot      = kwargs.pop('logPlot', True)
    self.nValidation  = kwargs.pop('nValidation', 0)
    self.percentages  = kwargs.pop('percentages', [0.5, 0.95, 0.99])
    self.preprocess   = kwargs.pop('preprocess', False)
    
    self.defaults = { 
            'R0'    : (1.0, 10.0),
            'D'     : (1.0, 10.0),
            'Z'     : (1.0, 11.0),
            'mu'    : (0.0, 1.0),
            'alpha' : (0.0, 1.0),
            'tact'  : (0.0, 100.0),
            'dtact' : (0.0, 14.0),
            'kbeta' : (0.0, 1.0),
            'Sigma' : (0.0, 100.0),
            'r'     : (0.0, 100.0)
        }

    self.constants = {
            'gamma' : 1.0/5.2,
            'Z'     : 2.7
    }
    
    super().__init__( **kwargs )
  
    if(self.synthetic):
        self.regionalData = SyntheticData( self.datafile )
    else:
        self.regionalData = RegionalData( self.country,self.preprocess )

    self.process_data()
 
  def save_data_path( self ):
      return ( self.dataFolder, self.country, self.modelName )
 
  def process_data( self ):
    y = self.regionalData.infected
    t = self.regionalData.time
    N = self.regionalData.populationSize

    I0 = y[0]
    S0 = N - I0
    y0 = S0, I0

    incidents = np.diff( y[0:] )
    if ((incidents < 0).any()):
        print("[Epidemics] Warning, removing negative values from daily infections!!!")

    incidents = np.clip(incidents, a_min=0, a_max=1e32)

    if self.nValidation == 0:
      self.data['Model']['x-data'] = t[1:]
      self.data['Model']['y-data'] = incidents
    else:
      self.data['Model']['x-data'] = t[1:-self.nValidation]
      self.data['Model']['y-data'] = incidents[0:-self.nValidation]
      self.data['Validation']['x-data'] = t[-self.nValidation:]
      self.data['Validation']['y-data'] = incidents[-self.nValidation-1:]

    self.data['Model']['Initial Condition'] = y0
    self.data['Model']['Population Size'] = self.regionalData.populationSize

    T = np.ceil( t[-1] + self.futureDays )
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
    incidents = np.diff(sol.y)
    incidents = np.append(ir0, incidents)
     
    eps = 1e-32
    incidents[incidents < eps] = eps
 
    recovered  = None
    exposed    = None
    unreported = None

    if hasattr(sol, 'r'):
        recovered = np.diff(sol.r)
        recovered = np.append(0, recovered)
 
    if hasattr(sol, 'e'):
        exposed = np.diff(sol.e)
        exposed = np.append(0, exposed)

    if hasattr(sol, 'iu'):
        unreported = np.diff(sol.iu)
        unreported = np.append(0, unreported)

    eps = 1e-32
    
    k = 0
    js = {}
    js['Variables'] = []

    js['Variables'].append({})
    js['Variables'][k]['Name'] = 'Daily Incidence'
    js['Variables'][k]['Values'] = list(incidents)
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

    js['Number of Variables'] = len(js['Variables'])
    js['Length of Variables'] = len(t)

    if self.likelihoodModel == 'Normal':
        js['Standard Deviation'] = ( p[-1] * incidents ).tolist()
    elif self.likelihoodModel == 'Positive Normal':
        js['Standard Deviation'] = ( p[-1] * incidents ).tolist()
    elif self.likelihoodModel == 'Negative Binomial':
        js['Dispersion'] = (len(incidents)) * [p[-1]]

    s['Saved Results'] = js
  
  def plot_intervals( self, ns=10):

    fig = self.new_figure()

    ax  = fig.subplots( 2 )

    ax[0].plot( self.data['Model']['x-data'], self.data['Model']['y-data'], 'o', lw=2, label='Daily Infected(data)', color='black')

    if self.nValidation > 0:
      ax[0].plot( self.data['Validation']['x-data'], self.data['Validation']['y-data'], 'x', lw=2, label='Daily Infected (validation data)', color='black')

    z = np.cumsum(self.data['Model']['y-data'])
    ax[1].plot( self.data['Model']['x-data'], z, 'o', lw=2, label='Cummulative Infected(data)', color='black')

    self.compute_plot_intervals( 'Daily Incidence', ns, ax[0], 'Daily Incidence' )
    self.compute_plot_intervals( 'Daily Incidence', ns, ax[1], 'Cummulative number of infected', cummulate=1)
    
    if 'Daily Recovered' in self.propagatedVariables:
      self.compute_mean_median( 'Daily Recovered', 'green', ns, ax[0], 'Daily Recovered')
      self.compute_mean_median( 'Daily Recovered', 'green', ns, ax[1], 'Cummulative number of recovered', cummulate=1)

    if 'Daily Exposed' in self.propagatedVariables:
      self.compute_mean_median( 'Daily Exposed', 'yellow', ns, ax[0], 'Daily Exposed')
      self.compute_mean_median( 'Daily Exposed', 'yellow', ns, ax[1], 'Cummulative number of exposed', cummulate=1)

    if 'Daily Unreported' in self.propagatedVariables:
      self.compute_mean_median( 'Daily Unreported', 'orange', ns, ax[0], 'Daily Unreported')
      self.compute_mean_median( 'Daily Unreported', 'orange', ns, ax[1], 'Cummulative number of unreported', cummulate=1)


    ax[-1].set_xlabel('time in days')

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

    # get incidents
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

    # get incidents
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
