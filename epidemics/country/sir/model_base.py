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

    self.country      = kwargs.pop('country', 'switzerland')
    self.futureDays   = kwargs.pop('futureDays', 2)
    self.nPropagation = kwargs.pop('nPropagation', 100)
    self.logPlot      = kwargs.pop('logPlot', False)
    self.nValidation  = kwargs.pop('nValidation', 0)
    self.percentages  = kwargs.pop('percentages', [0.5, 0.95, 0.99])
    self.preprocess   = kwargs.pop('preprocess', False)

    super().__init__( **kwargs )

    self.regionalData = RegionalData( self.country,self.preprocess)
    self.propagationData = {}


  def save_data_path( self ):
    if hasattr(self, 'property'):
      return ( self.dataFolder, self.country, self.modelName )
    else:
      return ( self.dataFolder, self.country )


  def solve_ode( self, y0, T, t_eval, N, p ):
    
    sir       = libepidemics.country.sir
    data      = libepidemics.country.ModelData(N=N)
    cppsolver = sir.Solver(data)

    params = sir.Parameters(beta=p[0], gamma=p[1])
    
    s0, i0 = y0
    y0cpp   = (s0, i0, 0.0)
    initial = sir.State(y0cpp)
    
    cpp_res = cppsolver.solve_ad(params, initial, t_eval=t_eval, dt = 0.01)
    
    yS      = np.zeros(len(cpp_res))
    gradmu  = []
    gradsig = []

    for idx,entry in enumerate(cpp_res):
        yS[idx] = entry.S().val()
        gradmu.append(np.array([ entry.S().d(0), entry.S().d(1), 0.0 ])) 
        gradsig.append(np.array([ 0.0, 0.0, 1.0 ]))

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

    ax[0].plot( self.data['Model']['x-data'], self.data['Model']['y-data'], 'o', lw=2, label='Daily Infected(data)', color='black')

    if self.nValidation > 0:
      ax[0].plot( self.data['Validation']['x-data'], self.data['Validation']['y-data'], 'x', lw=2, label='Daily Infected (validation data)', color='black')

    self.compute_plot_intervals( 'Daily Incidence', ns, ax[0], 'Daily Incidence' )

    #----------------------------------------------------------------------------------------------------------------------------------
    z = np.cumsum(self.data['Model']['y-data'])
    ax[1].plot( self.data['Model']['x-data'], z, 'o', lw=2, label='Cummulative Infected(data)', color='black')

    self.compute_plot_intervals( 'Daily Incidence', ns, ax[1], 'Cummulative number of infected', cummulate=1)

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
