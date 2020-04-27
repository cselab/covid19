import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from osp import *
import time
from mpl_toolkits.mplot3d import axes3d

c = 1


#####################
def plot_all_2d(Ny=1000, Ntheta=1000):
#####################  
  # v_list = utility[ number of sensors ] [ canton ] [ day ]
  
  #v_list = np.load("result_Ny{:05d}_Nt{:05d}.npy".format(Ny,Ntheta))
  v_list = np.load("result.npy")
  
  sensors = len(v_list)
  cantons = 26
  days    = len(v_list[0][0])

  v = np.zeros( (sensors,cantons*days) )
  for i in range(sensors):
    for j in range(cantons):
      for k in range(days  ):
        v[i][ j*days + k ] = v_list[i][j][k]

  name = ['AG','AI','AR','BE','BL','BS','FR','GE','GL','GR',\
          'JU','LU','NE','NW','OW','SG','SH','SO','SZ','TG',\
          'TI','UR','VD','VS','ZG','ZH']
  
  fig = plt.figure()
  ax1 = fig.add_subplot(111)

  locations = c * np.arange(0,cantons*days)
  for s in range(sensors):
    ax1.plot(locations,v[s,:], label= str(s+1) + " sensors ")


  locations1 = np.linspace(c*days/2, c*cantons*days + c*days/2 ,cantons+1)
  locations2 = np.linspace(0, c*cantons*days, cantons+1)

  ax1.set_xticks(locations1)
  ax1.set_xticklabels(name)
  
  legend_x = 1
  legend_y = 0.5
  plt.legend(loc='center left', bbox_to_anchor=(legend_x, legend_y))
  
  plt.grid()
  plt.savefig("cantons2d.eps", format='eps')


########################
def plot_all_3d(Ny=1000, Ntheta=1000):
########################
  
  #v_list = np.load("result_Ny{:05d}_Nt{:05d}.npy".format(Ny,Ntheta))
  v_list = np.load("result.npy")
 
  sensors = len(v_list)
  cantons = 26
  days    = len(v_list[0][0])

  ca = np.arange(0,cantons)
  da = np.arange(0,days)
  x,y = np.meshgrid( ca ,c*da )

  name = ['AG','AI','AR','BE','BL','BS','FR','GE','GL','GR',\
          'JU','LU','NE','NW','OW','SG','SH','SO','SZ','TG',\
          'TI','UR','VD','VS','ZG','ZH']


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


#####################
def plot_all_cantons(Ny=1000, Ntheta=1000):
#####################  
  # v_list = utility[ number of sensors ] [ canton ] [ day ]

  #v_list = np.load("result_Ny{:05d}_Nt{:05d}.npy".format(Ny,Ntheta)) 
  v_list = np.load("result.npy")
  
  sensors = len(v_list)
  cantons = 26
  days    = len(v_list[0][0])

  v = np.zeros( (sensors,cantons*days) )
  for i in range(sensors):
    for j in range(cantons):
      for k in range(days  ):
        v[i][ j*days + k ] = v_list[i][j][k]

  name = ['AG','AI','AR','BE','BL','BS','FR','GE','GL','GR',\
          'JU','LU','NE','NW','OW','SG','SH','SO','SZ','TG',\
          'TI','UR','VD','VS','ZG','ZH']

  max_v = np.max(v_list)
  locations = c*np.arange(0,days)
  
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
  return




  fig, axs = plt.subplots(3,5)

  
  for i0 in range (3):
    for i1 in range (5):
      index = i0 * 5 + i1
      for s in range(sensors):
        axs[i0,i1].plot(locations,v_list[s][index][:],label=str(s+1) + " sensors ")
      axs[i0,i1].set_title(name[index],fontsize='8')
      axs[i0,i1].grid()
      axs[i0,i1].set_ylim([0.0,max_v])

  for ax in axs.flat:
      ax.set(xlabel='Day', ylabel='Utility')
  for ax in axs.flat:
      ax.label_outer()

  plt.savefig("slice1.eps", format='eps')

  fig, axs = plt.subplots(3,4)  
  locations = c*np.arange(0,days)
  for i0 in range (3):
    for i1 in range (4):
      index = 15 + i0 * 4 + i1
      if index > 25:
        break
      for s in range(sensors):
        axs[i0,i1].plot(locations,v_list[s][index][:],label=str(s+1) + " sensors ")
      axs[i0,i1].set_title(name[index],fontsize='8')
      axs[i0,i1].grid()
      axs[i0,i1].set_ylim([0.0,max_v])

  for ax in axs.flat:
      ax.set(xlabel='Day', ylabel='Utility')  
  for ax in axs.flat:
      ax.label_outer()

  plt.savefig("slice2.eps", format='eps')



######################
def post_processing(osp):
######################  
  #plot_all_2d()
  plot_all_cantons()
  #plot_all_3d()
  
  #max_posterior = osp.computePosterior()
  #nSensors = args['nSensors']
  #R0 = np.arange(osp.Ntheta)
  #X, Y = np.meshgrid(R0, R0)
  #plt.contourf(X,Y, max_posterior/max_posterior.max(), cmap="hot")
  #plt.colorbar()
  #plt.xlabel("Possible values")
  #plt.ylabel("True value")
  #plt.title(str(nSensors) + "  sensors" )
  #plt.savefig("objective.png")
  #plt.close()



##########################
if __name__ == '__main__':
##########################  
  parser = argparse.ArgumentParser()

  parser.add_argument('--nSensors'    , help='number of sensors to place'                                      , required=True, type=int)
  parser.add_argument('--nMeasure'    , help='how many numbers describe a measurement taken by a single sensor', required=True, type=int)
  parser.add_argument('--Ntheta'      , help='number of samples', required=True, type=int)
  parser.add_argument('--path'        , help='path to files to perform OSP'                                    , required=True, type=str)
  args = vars(parser.parse_args())

  osp = OSP(path        =args['path']        ,\
            nSensors    =args['nSensors']    ,\
            Ntheta      =args['Ntheta']      ,\
            nMeasure    =args['nMeasure']     )
  
  timer = 0
  timer -= time.time()
  t,s = osp.Sequential_Placement()
  timer += time.time()
  print ("time=",timer)

  post_processing(osp)
