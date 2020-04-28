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


  max_v = np.max(v_list)
  locations = np.arange(0,days)
  
  fig, axs = plt.subplots(2,2)
  s = 0
  axs[0,0].plot(locations,v_list[s][0 ][:],label=name[0 ])
  axs[0,0].plot(locations,v_list[s][1 ][:],label=name[1 ])
  axs[0,0].plot(locations,v_list[s][2 ][:],label=name[2 ])
  axs[0,0].plot(locations,v_list[s][3 ][:],label=name[3 ])
  axs[0,0].plot(locations,v_list[s][4 ][:],label=name[4 ])
  axs[0,0].plot(locations,v_list[s][5 ][:],label=name[5 ])

  axs[0,1].plot(locations,v_list[s][6 ][:],label=name[6 ])
  axs[0,1].plot(locations,v_list[s][7 ][:],label=name[7 ])
  axs[0,1].plot(locations,v_list[s][8 ][:],label=name[8 ])
  axs[0,1].plot(locations,v_list[s][9 ][:],label=name[9 ])
  axs[0,1].plot(locations,v_list[s][10][:],label=name[10])
  axs[0,1].plot(locations,v_list[s][11][:],label=name[11])

  axs[1,0].plot(locations,v_list[s][12][:],label=name[12])
  axs[1,0].plot(locations,v_list[s][13][:],label=name[13])
  axs[1,0].plot(locations,v_list[s][14][:],label=name[14])
  axs[1,0].plot(locations,v_list[s][15][:],label=name[15])
  axs[1,0].plot(locations,v_list[s][16][:],label=name[16])
  axs[1,0].plot(locations,v_list[s][17][:],label=name[17])
  axs[1,0].plot(locations,v_list[s][18][:],label=name[18])

  axs[1,1].plot(locations,v_list[s][19][:],label=name[19])
  axs[1,1].plot(locations,v_list[s][20][:],label=name[20])
  axs[1,1].plot(locations,v_list[s][21][:],label=name[21])
  axs[1,1].plot(locations,v_list[s][22][:],label=name[22])
  axs[1,1].plot(locations,v_list[s][23][:],label=name[23])
  axs[1,1].plot(locations,v_list[s][24][:],label=name[24])
  axs[1,1].plot(locations,v_list[s][25][:],label=name[25])

  for i0 in range(2):
    for i1 in range(2):
      axs[i0,i1].grid()
      axs[i0,i1].legend(fontsize='xx-small')
      axs[i0,i1].set_ylim([0.0,max_v])
      axs[i0,i1].set_xlim([0.0,70.0 ])

  for ax in axs.flat:
      ax.set(xlabel='Day', ylabel='Information')
  for ax in axs.flat:
      ax.label_outer()

  plt.savefig("slice.eps", format='eps')


  fig, axs = plt.subplots(5,6)
  axs.titlesize      : xx-small
  for i0 in range (5):
    for i1 in range (6):
      index = i0 * 6 + i1
      if index > 25:
        break
      for s in range(sensors):
        axs[i0,i1].plot(locations,v_list[s][index][:],label=str(s+1) + " sensors ")
      axs[i0,i1].set_title(name[index],fontsize='4')
      axs[i0,i1].grid()
      axs[i0,i1].set_ylim([0.0,max_v])

  for ax in axs.flat:
      ax.set(xlabel='Day', ylabel='Utility')
  for ax in axs.flat:
      ax.label_outer()

  plt.savefig("slice1.eps",format='eps')





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
            texts[c] = str("{:.0f}".format(v_r2[c]))

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

  utility = np.load(args["result"])
  results = np.load(args["output"])

  plot_all_2d(utility)
  plot_all_3d(utility)
  plot_all_cantons(utility)
  
  #this fails on Euler!
  if args["movie"] == 1:
     res = results[0,:,:]
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
