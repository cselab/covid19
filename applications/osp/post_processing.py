
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from osp import *
import time
from mpl_toolkits.mplot3d import axes3d
import os
import sys
from plot import Renderer
from matplotlib.pyplot import figure


name = ['AG','AI','AR','BE','BL','BS','FR','GE','GL','GR',\
        'JU','LU','NE','NW','OW','SG','SH','SO','SZ','TG',\
        'TI','UR','VD','VS','ZG','ZH']

CANTON_TO_INDEX = {'AG': 0 , 'AI': 1 , 'AR': 2 , 'BE': 3 , 'BL': 4 , 'BS': 5 ,\
                   'FR': 6 , 'GE': 7 , 'GL': 8 , 'GR': 9 , 'JU': 10, 'LU': 11,\
                   'NE': 12, 'NW': 13, 'OW': 14, 'SG': 15, 'SH': 16, 'SO': 17,\
                   'SZ': 18, 'TG': 19, 'TI': 20, 'UR': 21, 'VD': 22, 'VS': 23,\
                   'ZG': 24, 'ZH': 25}
NUM_CANTONS = len(CANTON_TO_INDEX)


########################
def plot_all_2d(v_list):
########################  
  # v_list = utility[ number of sensors ] [ canton ] [ day ]
  sensors = len(v_list)
  cantons = 26
  days    = len(v_list[0][0])

  v = np.zeros( (sensors,cantons*days) )
  for i in range(sensors):
    for j in range(cantons):
      for k in range(days  ):
        v[i][ j*days + k ] = v_list[i][j][k]

  fig = plt.figure()
  ax1 = fig.add_subplot(111)

  locations = np.arange(0,cantons*days)
  for s in range(sensors):
    ax1.plot(locations,v[s,:], label= str(s+1) + " sensors ")

  locations1 = np.linspace(days/2, cantons*days + days/2 ,cantons+1)
  locations2 = np.linspace(0, cantons*days, cantons+1)

  ax1.set_xticks(locations1)
  ax1.set_xticklabels(name)
  
  legend_x = 1
  legend_y = 0.5
  plt.legend(loc='center left', bbox_to_anchor=(legend_x, legend_y))
  
  plt.grid()
  plt.savefig("cantons2d.eps", format='eps')


########################
def plot_all_3d(v_list):
########################
  sensors = len(v_list)
  cantons = 26
  days    = len(v_list[0][0])

  ca = np.arange(0,cantons)
  da = np.arange(0,days)
  x,y = np.meshgrid( ca ,da )

  fig = plt.figure()
  ax = fig.add_subplot(111, projection='3d')

  locations1 = np.linspace(0, cantons,cantons+1)
  ax.set_xticks(locations1)
  ax.set_xticklabels(name,fontsize='8')
  ax.set_xlabel("Canton")
  ax.set_ylabel("Day")
  ax.set_zlabel("Information Gain")

  for s in range(sensors):
    z = np.zeros((days,cantons))
    for i in ca:
      for j in da:
        z[j,i] = v_list[s][i][j]

    ax.plot_surface(x,y,z)
  plt.savefig("cantons3d.eps", format='eps')


#############################
def plot_all_cantons(v_list):
#############################
  # v_list = utility[ number of sensors ] [ canton ] [ day ]
  sensors = len(v_list)
  cantons = 26
  days    = len(v_list[0][0])

  v = np.zeros( (sensors,cantons*days) )
  for i in range(sensors):
    for j in range(cantons):
      for k in range(days  ):
        v[i][ j*days + k ] = v_list[i][j][k]


  print (days)
  max_v = np.max(v_list)
  locations = np.arange(1,days+1)
  locations_x = [0,20,40,60,80,100]
  
  locations_y = np.arange(0, int(max_v+1),2)

  fig, axs = plt.subplots(6,5)
  axs.titlesize      : xx-small
  for i0 in range (6):
    for i1 in range (5):
      index = i0 * 5 + i1
      print(i0,i1)
      if index > 25:
      	#fig.delaxes(axs[i0][i1])
        break
      else:
        for s in range(sensors):
          lab = str(s+1) + " sensors "
          if s == 0:
          	lab = str(s+1) + " sensor "
          axs[i0,i1].plot(locations,v_list[s][index][:],label=lab)
        axs[i0,i1].grid()
        axs[i0,i1].set_ylim([0.0,max_v])
        axs[i0,i1].text(.5,1.05,name[index],
            horizontalalignment='center',
            transform=axs[i0,i1].transAxes)
        axs[i0,i1].set_xticks(locations_x)
        axs[i0,i1].set_yticks(locations_y)

        axs[i0,i1].set(xlabel='Day', ylabel='Utility')
        if i0 != 4:
        	axs[i0,i1].label_outer()
  

  handles, labels = axs[4,1].get_legend_handles_labels()
  fig.legend(handles, labels, loc='lower center',ncol=sensors,bbox_to_anchor=(0.6, 0.1))
  
  #for ax in axs.flat:

  #for ax in axs.flat:
  #    ax.label_outer()

  i0 = 4
  axs[i0,0].set_xticks(locations_x)
  axs[i0,0].set_xlabel("")
  axs[i0,0].set_xticklabels("")

  for i1 in range(1,5):
        axs[i0,i1].set_xticks(locations_x)
        axs[i0,i1].set_ylabel("")
        axs[i0,i1].set_yticklabels("")




  for i0 in range (6):
      for i1 in range (5):
        index = i0 * 5 + i1
        print(i0,i1)
        if index > 25:
        	fig.delaxes(axs[i0][i1])
          #break
    

  #plt.show()
  

  #fig = plt.gcf()
  fig.set_size_inches(14.5, 10.5)
  fig.savefig('slice.pdf', dpi=100 ,format='pdf')
  #plt.savefig("slice1.eps",format='eps')





#################################
def make_movie(result,utility,n):
#################################
    """
    v = utility[number of sensors][canton][day]
    """
    days = utility.shape[2]

    def frame_callback(rend):
        t = rend.get_frame() * (days - 1) // rend.get_max_frame()
        util = utility[n,:,t]
        res  = result[t,:]
        max_res = np.max(result)
        
        v_u = {}
        v_r1 = {}
        v_r2 = {}
        texts = {}
        for i, c in enumerate(rend.get_codes()):
            i_state = CANTON_TO_INDEX[c]
            v_u[c] = util[i_state]
            v_r1[c] = res [i_state]/max_res
            v_r2[c] = res [i_state]
            tt = (np.asarray(v_r2[c])).astype(int)
            texts[c] = tt.astype(str)

        rend.set_values(v_r1)
        rend.set_texts(texts)
        rend.set_bars(v_u)
        plt.suptitle("Day :" + str(t), fontsize=12)
    rend = Renderer(frame_callback)
    rend.save_movie(frames=days,filename=str(n) + ".mp4")


##########################
if __name__ == '__main__':
##########################  
  parser = argparse.ArgumentParser()

  parser.add_argument('--result',help='utility result.npy file',type=str)
  parser.add_argument('--output',help='model evaluations output.npy file',type=str)
  parser.add_argument('--movie' ,help='make movie or not (1 or 0)',type=int)

  args = vars(parser.parse_args())

  print(args["result"])
  print(args["output"])

  utility = np.load(args["result"])

  plot_all_2d(utility)
  #plot_all_3d(utility)
  plot_all_cantons(utility)
  
  #this fails on Euler!
  if args["movie"] == 1:
     results = np.load(args["output"])
     #res = results[0,:,:]
     res = results
     for n in range(0,utility.shape[0]):
         make_movie(res,utility,n)

  '''
  max_posterior = osp.computePosterior()
  nSensors = args['nSensors']
  R0 = np.arange(osp.Ntheta)
  X, Y = np.meshgrid(R0, R0)
  plt.contourf(X,Y, max_posterior/max_posterior.max(), cmap="hot")
  plt.colorbar()
  plt.xlabel("Possible values")
  plt.ylabel("True value")
  plt.title(str(nSensors) + "  sensors" )
  plt.savefig("objective.png")
  plt.close()
  '''
