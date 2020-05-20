import os
import numpy as np

import matplotlib.pyplot as plt
plt.ioff()

from epidemics.epidemics import EpidemicsBase
from epidemics.data.combined import RegionalData
from epidemics.tools.tools import prepare_folder, save_file

import libepidemics #cpp backend

class Object(object):
        pass

class ModelBase( EpidemicsBase ):

  def __init__( self, **kwargs ):

    self.country        = kwargs.pop('country', 'switzerland')

    self.stdModel       = kwargs.pop('stdModel', 0)
    self.futureDays     = kwargs.pop('futureDays', 2)
    self.nPropagation   = kwargs.pop('nPropagation', 100)
    self.logPlot        = kwargs.pop('logPlot', False)
    self.nValidation    = kwargs.pop('nValidation', 0)
    self.percentages    = kwargs.pop('percentages', [0.5, 0.95, 0.99])

    super().__init__( **kwargs )

    self.regionalData = RegionalData( self.country )
    self.propagationData = {}

  def process_data( self ):

    y = self.regionalData.infected
    t = self.regionalData.time
    N = self.regionalData.populationSize

    Ir0 = y[0]

    S0  = N - Ir0
    E0  = 0
    Iu0 = 0
    y0  = S0, E0, Ir0, Iu0

    if self.nValidation == 0:
      self.data['Model']['x-data'] = t[1:]
      self.data['Model']['y-data'] = np.diff(y[0:])
    else:
      self.data['Model']['x-data'] = t[1:-self.nValidation]
      self.data['Model']['y-data'] = np.diff( y[0:-self.nValidation] )
      self.data['Validation']['x-data'] = t[-self.nValidation:]
      self.data['Validation']['y-data'] = np.diff( y[-self.nValidation-1:] )

    self.data['Model']['Initial Condition'] = y0
    self.data['Model']['Population Size'] = self.regionalData.populationSize
    self.data['Model']['Standard Deviation Model'] = self.stdModel

    T = np.ceil( t[-1] + self.futureDays )
    self.data['Propagation']['x-data'] = np.linspace(0,T,int(T+1))

    save_file( self.data, self.saveInfo['inference data'], 'Data for Inference', 'pickle' )


  def solve_ode( self, y0, T, t_eval, N, p ):

    seiir     = libepidemics.country.seiir
    data      = libepidemics.country.ModelData(N=N)
    cppsolver = seiir.Solver(data)

    params = seiir.Parameters(beta=p[0], mu=p[1], alpha=p[2], Z=p[3], D=p[4])
 
    s0, e0, ir0, iu0 = y0
    y0cpp   = (s0, e0, ir0, iu0, 0.0)
    
    initial = seiir.State(y0cpp)
 
    cpp_res = cppsolver.solve_ad(params, initial, t_eval=t_eval, dt = 0.01)
  
    yS      = np.zeros(len(cpp_res))
    gradmu  = []
    gradsig = []

    for idx,entry in enumerate(cpp_res):
        yS[idx] = entry.S().val()
        gradmu.append(np.array([ entry.S().d(0), entry.S().d(1), entry.S().d(2), entry.S().d(3), entry.S().d(4), 0.0 ])) 
        gradsig.append(np.array([ 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 ]))
 
    # Fix bad values
    yS[np.isnan(yS)] = 0
 
    # Create Solution Object
    sol = Object()
    sol.y       = [yS]
    sol.gradMu  = gradmu
    sol.gradSig = gradsig
 
    return sol


  def plot_intervals( self, ns=10):

    fig = self.new_figure()

    ax  = fig.subplots( 2 )

    ax[0].plot( self.data['Model']['x-data'], self.data['Model']['y-data'], 'o', lw=2, label='Daily Infected (data)', color='black')

    if self.nValidation > 0:
      ax[0].plot( self.data['Validation']['x-data'], self.data['Validation']['y-data'], 'x', lw=2, label='Daily Infected (validation data)', color='black')

    self.compute_plot_intervals( 'Daily Reported Incidence', ns, ax[0], 'Daily Reported Incidence' )

    #----------------------------------------------------------------------------------------------------------------------------------
    z = np.cumsum(self.data['Model']['y-data'])
    ax[1].plot( self.data['Model']['x-data'], z, 'o', lw=2, label='Cummulative Infected(data)', color='black')

    self.compute_plot_intervals( 'Daily Reported Incidence', ns, ax[1], 'Cummulative number of reported infected', cummulate=1)

    #----------------------------------------------------------------------------------------------------------------------------------

    ax[-1].set_xlabel('time in days')

    file = os.path.join(self.saveInfo['figures'],'prediction.png');
    prepare_folder( os.path.dirname(file) )
    fig.savefig(file)

    plt.show()

    plt.close(fig)

    #----------------------------------------------------------------------------------------------------------------------------------
    #----------------------------------------------------------------------------------------------------------------------------------

    fig = self.new_figure()

    ax  = fig.subplots( 1 )

    if self.parameters[0]['Name'] == 'R0':
      ax.hist( self.parameters[0]['Values'], bins=40, density=1)
    else:
      ax.hist( self.parameters[0]['Values']/self.parameters[1]['Values'], bins=40, density=1)

    file = os.path.join(self.saveInfo['figures'],'R0.png');
    fig.savefig(file)

    plt.show()

    plt.close(fig)