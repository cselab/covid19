import numpy as np
import scipy.stats as sp

class OSP:

  #########################################################################################################
  def __init__(self, path, nSensors = 1, Ntheta = 100, Ny = 100 ,korali = True, start_day = -1,days=-1):
  #########################################################################################################
    self.path       = path           # path to tensor_.npy
    self.nSensors   = nSensors       # how many sensors to place
    self.Ny         = Ny             # how many samples to draw when evaluating utility
    self.Ntheta     = Ntheta         # model parameters samples
    self.korali     = korali         # whether Korali will be used or not  
    self.l          = 7.0            # arbitrary choice for time correlation length 
    self.days       = days           # how many days (i.e. sensor locations)
    self.start_day  = start_day      
    self.sigma_mean = np.zeros(self.days)  

    #self.sigma_mean1 = np.zeros((self.days,26))

    self.sigma      = [0.025,0.075,0.125,0.175,0.225,0.275]
    self.wsigma     = (0.30-0.00)/len(self.sigma)
    self.wsigma    /= (0.30-0.00)

    #self.sigma      = [0.01,0.03,0.05,0.07,0.09]
    #self.wsigma     = (0.10-0.00)/len(self.sigma)
    #self.wsigma    /= (0.10-0.00)

    for i in range(self.days):
      temp = np.load(path+"/tensor_Ntheta={:05d}.npy".format(i))
      self.sigma_mean[i] = np.mean(temp.flatten())
      #for c in range(26):
      #    self.sigma_mean1[i,c] = np.mean(temp[:,:,c].flatten())
         
  
  #############################################################################################
  def distance(self,time1,space1,time2,space2):
  #############################################################################################
    s = space1 - space2
    s[ s != 0  ] = np.iinfo(np.int32).max #assume different locations in space are uncorrelated
    retval = np.abs(time1-time2)/self.l + s 
    return retval

  #############################################################################################
  def EvaluateUtility(self,argument):
  #############################################################################################
    space = []
    time  = []
    if self.korali:
      st = argument["Parameters"]
      assert( len(st)%2 == 0 )
      n = int ( len(st)/ 2 )
      for i in range(n*2):
        if i%2 == 0:
          space.append(int(st[-(i+1)]))
        else:
          time.append(int(st[-(i+1)]))
    else:
      n = int ( len(argument)/2 )
      for i in range(n):
        space.append(argument[i])
        time.append(argument[i+n])

    n = len(time)
    for i in range(n):
        if time[i] < self.start_day:
          argument["F(x)"] = 0.0
          return

    for i in range(n-1):
        if time[n-1] == time[i] and space[n-1] == space[i]:
           argument["F(x)"] = 0.0
           return 


           st = []
           for j in range(n-1):
               st.append(space[j])
           for j in range(n-1):
               st.append(time[j])
           self.korali = False
           retval = self.EvaluateUtility(st)
           self.korali = True
           if self.korali:
              argument["F(x)"] = retval
              return 
           else:
              return retval
    Ntheta   = self.Ntheta
    Ny       = self.Ny
    F_tensor = np.zeros( (Ntheta, Ntheta, n ))
    for s in range(n):
        temp = np.load(self.path+"/tensor_Ntheta={:05d}.npy".format(time[s]))
        F_tensor[:,:,s] = temp[:,:,space[s]]

    #Estimate covariance matrix as a function of the sensor locations (time and space)
    rv_list = []
    covariances = []
    T1,T2 = np.meshgrid(time ,time )
    X1,X2 = np.meshgrid(space,space)
    block = np.exp(-self.distance(T1,X1,T2,X2)) 
    cov   = np.kron(np.eye(1), block)
    sigma_mean = np.zeros(n)
    for i in range(n):
      sigma_mean[i] = self.sigma_mean[time[i]]
      #sigma_mean[i] = self.sigma_mean1[time[i],space[i]]
    sig  = sigma_mean * np.eye(n)
    aux  = sig*cov*sig
    for sigma_var in self.sigma:
      aux1 = (sigma_var*sigma_var)*aux
      rv_list.append(sp.multivariate_normal(np.zeros(n), aux1, allow_singular=True))
      covariances.append(aux1)

    #compute utility
    retval = 0.0
    for jjj in range(len(self.sigma)):
      aux  = covariances[jjj] 
      rv   = rv_list[jjj]
      for theta in range(Ntheta):


        mean = F_tensor[theta,theta,:]
        y    = np.random.multivariate_normal(mean=mean, cov=aux, size=Ny)  
        s1   = np.mean(rv.logpdf(y-mean))

        #this is a faster way to avoid a second for loop over Ntheta
        evidence = np.zeros((Ntheta,Ny))
        s2 = 0.0
        m1,m2 = np.meshgrid(y[:,0],F_tensor[:,theta,0] )
        new_shape1 = m1.shape[0]*m1.shape[1]
        new_shape2 = m2.shape[0]*m2.shape[1]
        m1 = m1.reshape((new_shape1,1))
        m2 = m2.reshape((new_shape2,1))
        for ns in range(1,n):
             m1_tmp,m2_tmp = np.meshgrid(y[:,ns],F_tensor[:,theta,ns] )
             m1_tmp = m1_tmp.reshape(m1_tmp.shape[0]*m1_tmp.shape[1],1)
             m2_tmp = m2_tmp.reshape(m2_tmp.shape[0]*m2_tmp.shape[1],1)
             m1=np.concatenate( (m1,m1_tmp), axis= 1 )
             m2=np.concatenate( (m2,m2_tmp), axis= 1 )
        for i_sigma in range(len(self.sigma)):
          evidence += (rv_list[i_sigma].pdf(m1-m2)).reshape((Ntheta,Ny))
        s2 += np.mean ( np.log( self.wsigma*np.mean( evidence,axis=0) ) )

        retval += (s1-s2)
    retval *= self.wsigma/self.Ntheta 
    if self.korali:
      print ("Retval=",retval,"space=",space,"time=",time,flush=True)
      argument["F(x)"] = retval
    else:
      return retval
