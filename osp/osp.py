import argparse
import numpy as np
import scipy.stats as sp
import matplotlib.pyplot as plt
import itertools

class OSP:

  ############################################################################
  def __init__(self, path, nSensors, nMeasure, Ny = 100, korali = False):
  ############################################################################
    self.path         = path         # path to output.npy
    self.nSensors     = nSensors     # how many sensors to place
    self.nMeasure     = nMeasure     # how many quantities each sensor measures
    self.Ny = Ny                     # how many samples to draw when evaluating utility

    #output.npy should contain a 3D numpy array [Simulations][Time][Space]
    self.data         = np.load(path+"/output.npy")
    self.parameters   = np.load(path+"/params.npy")
    self.Ntheta       = self.data.shape[0] # number of simulations 

    self.korali = korali

    self.l = 3
    assert nMeasure == 1
    self.sigma = 0.2 #np.std(self.data.flatten())


    self.sigma_mean = self.sigma * np.mean ( np.mean(self.data,axis=0) , axis = 1)



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

      p_theta = 1.0
      #if self.UseReported:
      #  t_tilde = np.max(time)
      #  if t_tilde > 0:
      #    cantons = 26
      #    rv1 = sp.multivariate_normal(np.zeros((t_tilde+1) * cantons), self.var[theta][t_tilde], allow_singular=True)
      #    p_theta = rv1.pdf(np.zeros((t_tilde+1) * cantons)) + 1e-8

      mean = F[theta][:]



      #sig = (self.sigma*mean)**2 * np.eye(N)
      #sig = (self.sigma)**2 * np.eye(N)

      #sig = sigma_mean**2 * np.eye(N)
      sig = sigma_mean * np.eye(N)


      #rv  = sp.multivariate_normal(np.zeros(n), sig*cov,allow_singular=True)
      #y   = np.random.multivariate_normal(mean, sig*cov, self.Ny)  
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


'''
Slow version of EvaluateUtility and area bombing optimization  

  #####################################
  def EvaluateUtility(self,space_time):
  #####################################
    n = int ( len(space_time)/2 )
    space = []
    time  = []
    for i in range(n):
      space.append(space_time[i])
      time. append(space_time[i+n])


    #time is an array containing the time instance of each one of the nSensors measurements
    #space contains the place where each measurement happens
    M = self.nMeasure
    N = len(time)  #self.nSensors

    F = np.zeros((self.Ntheta,  M*N ))
    for theta in range(0,self.Ntheta):
      for s in range(0,N):
        F[theta][s*M:(s+1)*M] = self.data[theta][time[s]][ self.nMeasure*space[s] : self.nMeasure*(space[s]+1) ] 

    retval = 0.0
    T1,T2 = np.meshgrid(time ,time )
    X1,X2 = np.meshgrid(space,space)
    block = self.sigma * self.sigma * np.exp( -self.distance(T1,X1,T2,X2) )      
    cov = np.kron(np.eye(self.nMeasure,dtype=int), block)

    if np.abs( np.linalg.det(cov) ) < 1e-7:
      cov[np.diag_indices_from(cov)] += 1e-3

    if np.abs( np.linalg.det(cov) ) < 1e-3:
      return 0.0

    for theta in range(0,self.Ntheta):
      mean = F[theta][:]  
      y = np.random.multivariate_normal(mean, cov, self.Ny)

      rv = sp.multivariate_normal(mean, cov,allow_singular=True)
      s1 = np.mean(rv.logpdf(y))
      
      Pdf_inner = np.empty((self.Ntheta,self.Ny))
      for k in range(0,self.Ntheta):
        rvInner = sp.multivariate_normal(F[k][:], cov,allow_singular=True)
        Pdf_inner[k,:] = rvInner.pdf(y)
      s2 = np.mean ( np.log( np.mean( Pdf_inner,axis=0) ) )

      retval += (s1-s2)

    retval /= self.Ntheta
    
    return retval

  ######################
  def AreaBombing(self):
  ######################
    t_locations = self.data.shape[1]
    x_locations = int(self.data.shape[2] / self.nMeasure)

    t = np.arange(0,t_locations)
    x = np.arange(0,x_locations)

    self.value = np.zeros( (  t_locations**self.nSensors , x_locations**self.nSensors ) )

    max_v = -1e50
    max_t = np.zeros(self.nSensors)
    max_s = np.zeros(self.nSensors) 
    min_v =  1e50
    min_t = np.zeros(self.nSensors)
    min_s = np.zeros(self.nSensors) 


    time  = np.zeros(self.nSensors)
    space = np.zeros(self.nSensors)
    
    #those loops are wrong because they do not allow to place two sensors
    # at different locations but at the same time
    for s_x in itertools.combinations_with_replacement( range(x_locations), self.nSensors ):
      for s_t in itertools.combinations( range(t_locations), self.nSensors ):
        time = s_t
        space = s_x
        t_index = self.index(s_t,t_locations)
        x_index = self.index(s_x,x_locations)
        self.value[t_index][x_index] = self.EvaluateUtility(space + time)
        print (space[0],time[0],self.value[t_index][x_index],flush=True)
        if self.value[t_index][x_index] > max_v:
          max_t = time
          max_s = space
          max_v = self.value[t_index][x_index]

        if self.value[t_index][x_index] < min_v:
          min_t = time
          min_s = space
          min_v = self.value[t_index][x_index]
      print(' ',flush=True)
      print(' ',flush=True)


    print ("Maximum utility:", max_v)
    print ("Optimal sensor locations")
    for n in range (self.nSensors):
      print ("Sensor",n,"in location ",max_s[n],"at time",max_t[n])

    return max_t,max_s,max_v
'''
