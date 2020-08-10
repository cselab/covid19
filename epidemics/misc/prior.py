#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

import scipy.stats
import scipy.integrate as integrate

# Assuming an incubation period distribution of mean 5.2 days from a separate study of early COVID-19 cases1, 
# we inferred that infectiousness started from 2.3 days (95% CI, 0.8–3.0 days) before 
# symptom onset and peaked at 0.7 days (95% CI, −0.2–2.0 days) before symptom onset (Fig. 1c). 

# The mean incubation period was 5.2 days (95% confidence interval [CI], 4.1 to 7.0), 
# with the 95th percentile of the distribution at 12.5 days.

if __name__ == '__main__':
        ip_mean = 5.2
        ip_sdev = 2.8
        ip_var  = ip_sdev**2
    
        # stats gamma, a (shape, k) scale (scale, theta)
        # mean     = k*theta
        # variance = k*theta^2

        # k*theta = 5.2
        # k*theta**2 = 2.8

        th_ip = ip_var/ip_mean
        k_ip  = ip_mean/th_ip
        
        x0 = np.linspace(0,50,1000)
        y0 = scipy.stats.gamma.pdf(x0, a=k_ip, scale=th_ip) 

        pct = scipy.stats.gamma.cdf(12.5, a=k_ip, scale=th_ip)
        cup  = scipy.stats.gamma.cdf(7.0, a=k_ip, scale=th_ip)
        clow = scipy.stats.gamma.cdf(4.1, a=k_ip, scale=th_ip)
    
        print("INCUBATION PERIOD")
        print("Percentile {} should be 95%?".format(pct))
        print("Conf intervall {} should be 95%?".format(cup-clow))
        
        plt.title('Estimated incubaition period PDF')
        plt.plot(x0,y0)
        plt.show()



        ib_mean=2.3
        ib_cup=3.0

        print("INFECTIOUSNESS BEFORE")
        print("SEARCH PARAM")
        for k in np.linspace(32.60, 32.65, 20):
            th  = ib_mean/k
            cup = scipy.stats.gamma.cdf(ib_cup, a=k, scale=th)
            print("k {} - th {}: conf. {}".format(k, th, cup))


        k_ib  = 32.62105263157895
        th_ib = 0.07050661503710874

        pct_ib = scipy.stats.gamma.cdf(3.0, a=k_ib, scale=th_ib)
        print("Percentile {} should be 95%?".format(pct_ib))
    
        pct_max = scipy.stats.gamma.cdf(ip_mean, a=k_ib, scale=th_ib)
        print("Percentile {} at {} (~max)".format(pct_max,ip_mean))
        
        x1 = np.linspace(0,10,1000)
        y1 = scipy.stats.gamma.pdf(x1, a=k_ib, scale=th_ib) 
        plt.title('Estimated preasymtomatic period PDF')
        plt.plot(x1,y1)
        plt.show()


 
        print("SIMULATE LATENCY PERIOD")
        rvs_ip = scipy.stats.gamma.rvs(a=k_ip, scale=th_ip, size=100000)
        rvs_ib = scipy.stats.gamma.rvs(a=k_ib, scale=th_ib, size=100000)
        rvs_lat = rvs_ip-rvs_ib
        plt.title('Simulated incubation/preasymptomatic/latency periods')
        plt.hist(rvs_ip, density=False, histtype='stepfilled', bins=100, alpha=0.2, label='INC')
        plt.hist(rvs_ib, density=False, histtype='stepfilled', bins=100, alpha=0.2, label='PRE')
        plt.hist(rvs_lat, density=False, histtype='stepfilled', bins=100, alpha=0.2, label='LAT')
        plt.legend(['INC','PRE','LAT'])
        plt.show()
        
        lat_mu = np.mean(rvs_lat)
        lat_pct = np.percentile(rvs_lat, 95)

        print("Simulated mean latency period {}".format(lat_mu))
        print("Simulated 95-percentile latency period {}".format(lat_pct ))


        print("SEARCH PARAM LATENCY PERIOD")
        for k in np.linspace(1.166, 1.169, 20):
            th  = lat_mu/k
            cup = scipy.stats.gamma.cdf(lat_pct, a=k, scale=th)
            print("k {} - th {}: conf. {}".format(k, th, cup))
        k_lat  = 1.1671052631578946
        th_lat = 2.4931393061857263
 
        pct_lat_check = scipy.stats.gamma.cdf(lat_pct, a=k_lat, scale=th_lat)
        print("Percentile {} should be 95%?".format(pct_lat_check))
 
        x2 = np.linspace(0,10,1000)
        y2 = scipy.stats.gamma.pdf(x2, a=k_lat, scale=th_lat) 
        plt.hist(rvs_lat, density=True, histtype='stepfilled', bins=100, alpha=0.2, label='LAT')
        plt.title('Estimated latency period PDF')
        plt.plot(x2,y2)
        plt.show()

        print("\n\n")
        print("X"*20)
        print("INCUBATION PERIOD")
        print("Shape {} - Scale {}".format(k_ip, th_ip))
        print("X"*20)
        print("PREASYMPTOMATIC PERIOD")
        print("Shape {} - Scale {}".format(k_ib, th_ib))
        print("X"*20)
        print("LATENCY PERIOD")
        print("Shape {} - Scale {}".format(k_lat, th_lat))
        print("X"*20)


        th_rec = 3.970
        k_rec  = 4.337
        
        x3 = np.linspace(0,50,1000)
        y3 = scipy.stats.gamma.pdf(x3, a=k_rec, scale=th_rec) 

        pct = scipy.stats.gamma.cdf(14.0, a=k_rec, scale=th_rec)
    
        mean, var = scipy.stats.gamma.stats(a = k_rec, scale=th_rec)
        print("RECOVERY PERIOD")
        print("Percentile {} should be 50%?".format(pct))
        print("Mean {} Var {}".format(mean, var))

        lam_rec = np.log(2)/14
        
        x4 = np.linspace(0,50,1000)
        y4 = scipy.stats.expon.pdf(x4,scale=1./lam_rec) 

        pct21 = scipy.stats.expon.cdf(21.0, scale=1./lam_rec)
        pct42 = scipy.stats.expon.cdf(42.0, scale=1./lam_rec)
        pct50 = scipy.stats.expon.cdf(50.0, scale=1./lam_rec)
    
        print("RECOVERY PERIOD")
        print("Percentile {} / {} at 3 and 6 weeks".format(pct21, pct42))
        print("Percentile {} at 50 days".format(pct50))
        
        plt.title('Estimated recovery period PDF')
        plt.plot(x3,y3)
        plt.plot(x4,y4)
        plt.show()


