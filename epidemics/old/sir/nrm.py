#!/usr/bin/env python3
# Author: George Arampatzis
# Date:   27/3/2020
# Email:  garampat@ethz.ch

import os

from scipy.integrate import solve_ivp
import numpy as np

from epidemics.utils.misc import prepare_folder, save_file
from .model_base import ModelBase
import epidemics.ode_solver as solver


class Model( ModelBase ):


  def __init__( self, **kwargs ):

    self.modelName        = 'sir.nrm'
    self.modelDescription = 'Fit SIR on Daily Infected Data with Normal Likelihood'
    self.likelihoodModel  = 'Normal'

    super().__init__( **kwargs )

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

    if self.nValidation == 0:
      self.data['Model']['x-data'] = t[1:]
      self.data['Model']['y-data'] = np.diff( y[0:])
    else:
      self.data['Model']['x-data'] = t[1:-self.nValidation]
      self.data['Model']['y-data'] = np.diff( y[0:-self.nValidation] )
      self.data['Validation']['x-data'] = t[-self.nValidation:]
      self.data['Validation']['y-data'] = np.diff( y[-self.nValidation-1:] )

    self.data['Model']['Initial Condition'] = y0
    self.data['Model']['Population Size'] = self.regionalData.populationSize

    T = np.ceil( t[-1] + self.futureDays )
    self.data['Propagation']['x-data'] = np.linspace(0,T,int(T+1))

    save_file( self.data, self.saveInfo['inference data'], 'Data for Inference', 'pickle' )




  def get_variables_and_distributions( self ):

    p = ['beta','gamma','[Sigma]']
    js = {}
    js['Variables']=[]
    js['Distributions']=[]

    for k,x in enumerate(p):
      js['Variables'].append({})
      js['Variables'][k]['Name'] = x
      js['Variables'][k]['Prior Distribution'] = 'Prior for ' + x

    self.nParameters = len(p)

    k=0
    js['Distributions'].append({})
    js['Distributions'][k]['Name'] = 'Prior for beta'
    js['Distributions'][k]['Type'] = 'Univariate/Uniform'
    js['Distributions'][k]['Minimum'] = 1
    js['Distributions'][k]['Maximum'] = 100.

    k+=1
    js['Distributions'].append({})
    js['Distributions'][k]['Name'] = 'Prior for gamma'
    js['Distributions'][k]['Type'] = 'Univariate/Uniform'
    js['Distributions'][k]['Minimum'] = 1.
    js['Distributions'][k]['Maximum'] = 100.

    k+=1
    js['Distributions'].append({})
    js['Distributions'][k]['Name'] = 'Prior for [Sigma]'
    js['Distributions'][k]['Type'] = 'Univariate/Uniform'
    js['Distributions'][k]['Minimum'] = 0.01
    js['Distributions'][k]['Maximum'] = 10.

    return js




  def computational_model( self, s ):
    p = s['Parameters']
    p[0] = p[0]/p[1]
    t  = self.data['Model']['x-data']
    y0 = self.data['Model']['Initial Condition']
    N  = self.data['Model']['Population Size']

    tt = [t[0]-1] + t.tolist()
    sol = solver.solve_ode(self.sir_rhs,T=t[-1],y0=y0,args=(N,p),t_eval = tt,backend='numpy')    
    y = -(sol.y[0][1:]-sol.y[0][:-1])
    # Get gradients here
    y = solver.to_list(y)

    s['Reference Evaluations'] = y
    s['Standard Deviation'] = ( p[-1] * np.maximum(np.abs(y),1e-4) ).tolist()




  def computational_model_propagate( self, s ):
    p = s['Parameters']
    p[0] = p[0]/p[1]
    t  = self.data['Propagation']['x-data']
    y0 = self.data['Model']['Initial Condition']
    N  = self.data['Model']['Population Size']

    sol = solver.solve_ode(self.sir_rhs,T=t[-1],y0=y0,args=(N,p),t_eval = t.tolist(),backend='numpy')
    y = -(sol.y[0][1:]-sol.y[0][:-1])
    y = solver.append_zero(y)
    y = solver.to_list(y)

    js = {}
    js['Variables'] = []

    js['Variables'].append({})
    js['Variables'][0]['Name']   = 'Daily Incidence'
    js['Variables'][0]['Values'] = y

    js['Number of Variables'] = len(js['Variables'])
    js['Length of Variables'] = len(t)

    js['Standard Deviation'] = ( p[-1] * np.maximum(np.abs(y),1e-4) ).tolist()

    s['Saved Results'] = js


#   def plot_intervals( self, ns=10):

#     fig = self.new_figure()

#     ax  = fig.subplots( 2 )

#     ax[0].plot( self.data['Model']['x-data'], self.data['Model']['y-data'], 'o', lw=2, label='Daily Infected(data)', color='black')

#     if self.nValidation > 0:
#       ax[0].plot( self.data['Validation']['x-data'], self.data['Validation']['y-data'], 'x', lw=2, label='Daily Infected (validation data)', color='black')

#     self.compute_plot_intervals( 'Daily Incidence', ns, ax[0], 'Daily Incidence' )

#     #----------------------------------------------------------------------------------------------------------------------------------
#     z = np.cumsum(self.data['Model']['y-data'])
#     ax[1].plot( self.data['Model']['x-data'], z, 'o', lw=2, label='Cummulative Infected(data)', color='black')

#     self.compute_plot_intervals( 'Daily Incidence', ns, ax[1], 'Cummulative number of infected', cummulate=1)

#     #----------------------------------------------------------------------------------------------------------------------------------
#     ax[-1].set_xlabel('time in days')

#     file = os.path.join(self.saveInfo['figures'],'prediction.png');
#     prepare_folder( os.path.dirname(file) )
#     fig.savefig(file)

#     plt.show()

#     plt.close(fig)
# =======
# >>>>>>> ea23ab34b3e8ab122a582dab48c67ff1e31fa5cd:epidemics/sir/nrm.py
