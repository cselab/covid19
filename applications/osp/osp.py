import argparse
import numpy as np
import scipy.stats as sp
import matplotlib.pyplot as plt
import itertools

class OSP:

  ##################################################################################
  def __init__(self, path, nSensors, nMeasure, Ntheta=0, Ny = 100, korali = False, informed_prior = True):
  ##################################################################################
    self.path         = path         # path to output.npy
    self.nSensors     = nSensors     # how many sensors to place
    self.nMeasure     = nMeasure     # how many quantities each sensor measures
    self.Ny = Ny                     # how many samples to draw when evaluating utility

    #output.npy should contain a 3D numpy array [Simulations][Time][Space]
    self.data = np.empty(0)
    self.paramters = np.empty(0)
    if Ntheta == 0:
        self.Ntheta      = self.data.shape[0] # number of simulations
        self.data        = np.load(path+"/output_Ntheta={:05d}.npy".format(Ntheta))
        self.parameters  = np.load(path+"/params_Ntheta={:05d}.npy".format(Ntheta))
    else:
        self.Ntheta      = Ntheta 
        self.data        = np.load(path+"/output_Ntheta={:05d}.npy".format(Ntheta))
        self.parameters  = np.load(path+"/params_Ntheta={:05d}.npy".format(Ntheta))

    self.korali = korali

    self.l = 3
    assert nMeasure == 1
    self.sigma = 0.2 


    self.sigma_mean = self.sigma * np.mean ( np.mean(self.data,axis=0) , axis = 1)


    self.prior = np.zeros(self.Ntheta)
    if informed_prior == True:
#       real_data  = np.load("canton_daily_cases.npy") # Daily cases  [canton][day]
#       model_data = np.load(path + "/reported.npy")    # Model output for daily cases  [Simulations][Days][canton]
#
#       
#       cantons = real_data.shape[0]
#       days    = real_data.shape[1]
#       assert cantons == 26
#
#       E = []
#       m = []
#
#       diag = []
#       for c in range(cantons):
#           for d in range(days):
#              if np.isnan(real_data[c,d]) == False:
#                 E.append(real_data[c,d])
#                 m_ = np.zeros(self.Ntheta)
#                 for theta in range(self.Ntheta):
#                     m_[theta] = model_data[theta,d,c]
#                 diag.append( np.mean(m_) )
#                 m.append(m_)
#
#       E = np.asarray(E)
#       m = np.asarray(m)
#       diag = ( 0.3*np.asarray(diag)**2 ) * np.eye(len(E))
#
#       rv  = sp.multivariate_normal(E , diag ,allow_singular=True)
#       for theta in range(self.Ntheta):
#           self.prior[theta] = rv.pdf( m[:,theta]  )
#
#       print (np.sum(self.prior))
       self.prior = 1.0
    else:
       self.prior = 1.0

#
#
#    reported = np.load(path+"/reported.npy")
#    simulations = reported.shape[0]
#    days        = reported.shape[1]
#    cantons     = reported.shape[2]
#    self.var = []
#    for theta in range (self.Ntheta):
#      print (theta,self.Ntheta)
#      v_theta = []
#      for d in range(0,days):
#        time  = np.repeat(np.arange(d+1),cantons)
#        space = np.tile(np.arange(cantons),d+1)
#        T1,T2 = np.meshgrid(time ,time )
#        X1,X2 = np.meshgrid(space,space)
#        
#        mean = reported[theta,0:d+1,0:cantons]
#        mean = mean.flatten()
#
#        #sig = np.eye(cantons*d)
#        #for c in range(cantons):
#        #  for d1 in range(d):
#        #    i  = c + d1 * cantons
#        #    #sig[i,i] = 1.0 #0.1*np.mean(reported[:,d1,c])
#  
#        block = (self.sigma*mean)**2  * np.exp( -self.distance(T1,X1,T2,X2) ) 
#        cov   = np.kron(np.eye(self.nMeasure), block)
#        v_theta.append(cov)
#
#      self.var.append(v_theta)
#    self.UseReported = False #True

  
  #####################################
  def EvaluateUtility(self,space_time):
  #####################################
    #time is an array containing the time instance of each one of the nSensors measurements
    #space contains the place where each measurement happens
    space = []
    time  = []
    if self.korali:
      st = space_time["Parameters"]
      assert( len(st)%2 == 0 )
      n = int ( len(st)/ 2 )
      for i in range(n*2):
        if i%2 == 0:
          space.append(int(st[-(i+1)]))
        else:
          time.append(int(st[-(i+1)]))
    else:
      n = int ( len(space_time)/2 )
      for i in range(n):
        space.append(space_time[i])
        time.append(space_time[i+n])

    M = self.nMeasure
    N = len(time)
    F = np.zeros((self.Ntheta,  M*N ))
    for s in range(0,N):
      F[:,s*M:(s+1)*M] = self.data[:,time[s], self.nMeasure*space[s] : self.nMeasure*(space[s]+1) ] 


    #Estimate covariance matrix as a function of the sensor locations (time and space)
    T1,T2 = np.meshgrid(time ,time )
    X1,X2 = np.meshgrid(space,space)
    block = np.exp( -self.distance(T1,X1,T2,X2) ) 
    cov   = np.kron(np.eye(self.nMeasure), block)
    #block = self.sigma * self.sigma * np.exp( -self.distance(T1,X1,T2,X2) )

    #if covariance is singular, add small number to diagonal to make it non-singular
    #if np.abs( np.linalg.det(cov) ) < 1e-7:
    #  cov[np.diag_indices_from(cov)] += 1e-3
    ##if this didn't work, just return 0...
    #if np.abs( np.linalg.det(cov) ) < 1e-3:
    #  return 0.0

    sigma_mean = np.zeros(N)
    for i in range(N):
      sigma_mean[i] = self.sigma_mean[time[i]]
    #sig = np.mean(sigma_mean)
    
    #compute utility
    retval = 0.0
    for theta in range(0,self.Ntheta):

      p_theta = 1.0 #self.prior[theta]
      #p_theta = 1.0
      #if self.UseReported:
      #  t_tilde = np.max(time)
      #  if t_tilde > 0:
      #    cantons = 26
      #    rv1 = sp.multivariate_normal(np.zeros((t_tilde+1) * cantons), self.var[theta][t_tilde], allow_singular=True)
      #    p_theta = rv1.pdf(np.zeros((t_tilde+1) * cantons)) + 1e-8

      mean = F[theta][:]




      sig = sigma_mean * np.eye(N)
      rv  = sp.multivariate_normal(np.zeros(n), sig*cov*sig,allow_singular=True)
      y   = np.random.multivariate_normal(mean, sig*cov*sig, self.Ny)  
      s1  = np.mean(rv.logpdf(y-mean))    
      
      #this is a faster way to avoid a second for loop over Ntheta
      m1,m2 = np.meshgrid(y[:,0],F[:,0] )
      m1 = m1.reshape((m1.shape[0]*m1.shape[1],1))
      m2 = m2.reshape((m2.shape[0]*m2.shape[1],1))
      for ns in range(1,n):
        m1_tmp,m2_tmp = np.meshgrid(y[:,ns],F[:,ns] )
        m1_tmp = m1_tmp.reshape(m1_tmp.shape[0]*m1_tmp.shape[1],1)
        m2_tmp = m2_tmp.reshape(m2_tmp.shape[0]*m2_tmp.shape[1],1)
        m1=np.concatenate( (m1,m1_tmp), axis= 1 )
        m2=np.concatenate( (m2,m2_tmp), axis= 1 )
      Pdf_inner = np.empty((self.Ntheta*self.Ny))
      Pdf_inner = rv.pdf(m1-m2)
      Pdf_inner = Pdf_inner.reshape((self.Ntheta,self.Ny))

      s2 = np.mean ( np.log( p_theta*np.mean( Pdf_inner,axis=0) ) )
      retval += (s1-s2)*p_theta

    retval /= self.Ntheta
    print ("Sample ok:",retval)
    if self.korali:
      space_time["F(x)"] = retval
    else:
      return retval

  #############################################
  def distance(self,time1,space1,time2,space2):
  #############################################  
    s = space1 - space2
    s[ s != 0  ] = np.iinfo(np.int32).max #assume different locations in space are uncorrelated
    retval = np.abs(time1-time2)/self.l + s 
    return retval

  ####################
  def index(self,x,d):
  ####################  
    index = 0
    for i in range(0,self.nSensors):
      index += x[i] * d**i 
    return index

  ###############################
  def Sequential_Placement(self):
  ###############################

    np.random.seed(12345)
    t_locations = self.data.shape[1]
    x_locations = int(self.data.shape[2] / self.nMeasure)
    t = np.arange(0,t_locations)
    x = np.arange(0,x_locations)

    t_sensors = []
    x_sensors = []
    v_sensors = []
    self.v_all = []

    counter = 0


    print ("Total evaluations required= ",t_locations*x_locations*self.nSensors )

    for n in range(self.nSensors):

      print ("Placing sensor",n+1,"of",self.nSensors)

      vmax = -1e50
      xmax = -1
      tmax = -1

      t_sensors.append(tmax)
      x_sensors.append(xmax)
      v_sensors.append(vmax)

      v_loc1 = []

      for s_x in x:
        
        v_loc2 = []
        for s_t in t:

          counter += 1

          t_sensors[n] = s_t
          x_sensors[n] = s_x

          v = self.EvaluateUtility(x_sensors + t_sensors)

          print (counter,"x=",s_x,"t=",s_t,"utility=",v)
          
          if v > vmax:
            t_ok = True
            x_ok = True
            for i in range (n-1):
              t_ok = (t_sensors[i] != s_t)
              x_ok = (x_sensors[i] != s_x)
            if x_ok and t_ok:
              vmax = v
              xmax = s_x
              tmax = s_t
      
          v_loc2.append(v)
        v_loc1.append(v_loc2)
      self.v_all.append(v_loc1)

      t_sensors[n] = tmax
      x_sensors[n] = xmax
      v_sensors[n] = vmax

      print ("Placed sensor",n,"at x=",xmax,"t=",tmax,'v=',self.EvaluateUtility(x_sensors + t_sensors))

    max_v = self.EvaluateUtility(x_sensors + t_sensors)

    print ("Maximum utility:", max_v)
    print ("Optimal sensor locations")
    for n in range (self.nSensors):
      print ("Sensor",n,"in location ",x_sensors[n],"at time",t_sensors[n])


    np.save("result_time.npy", t_sensors)
    np.save("result_space.npy", x_sensors)
    np.save("result.npy", self.v_all)
    
    return t_sensors,x_sensors

  ##########################################
  def computePosterior(self):
  ##########################################  
    time = np.load("result_time.npy")
    space = np.load("result_space.npy")
    posterior = np.zeros((self.Ntheta,self.Ntheta))
   
    M = self.nMeasure
    N = self.nSensors
    F = np.zeros((self.Ntheta,  M*N ))
 
    T1,T2 = np.meshgrid(time ,time )
    X1,X2 = np.meshgrid(space,space)
    block = self.sigma * self.sigma * np.exp( -self.distance(T1,X1,T2,X2) )      
    cov = np.kron(np.eye(self.nMeasure,dtype=int), block)
 
    for theta in range(0,self.Ntheta):
      for s in range(0,self.nSensors):
        F[theta][s*M:(s+1)*M] = self.data[theta][time[s]][ self.nMeasure*space[s] : self.nMeasure*(space[s]+1) ] 

    for i in range(0,self.Ntheta):
      mean = F[i][:]
      for j in range(0,self.Ntheta):
        posterior[i,j] = sp.multivariate_normal.pdf(F[j][:] ,mean, cov)

    return posterior      
