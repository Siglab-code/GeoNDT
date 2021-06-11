==============
Usage examples
==============

Bender Element Test
==============

The GeoNDT can efficiently study the three-dimensional wave propagation within soil specimens 
in the Bender element (BE) test.  In this example,  the soil sample is assumed to be 7.0 cm in diameter and 14 cm 
in height. The density of the dry sand is 1,800 kg/m$^3$. The P-wave and S-wave velocity of the 
soil sample is assumed as 380 m/s and 240 m/s, respectively. The code example for BE modeling is shown as follows: 

    >>> import numpy as np  
    >>> from geondt import one_phase_dynamic   
    >>> tlin = 200; tmin = 3*10**(-4); tmax = 150*10**(-5)   
    >>> vp = np.array([380.0])
    >>> vs = np.array([240.0])
    >>> rho  = np.array([1800])
    >>> E = rho1*vs**2*(3*vp**2-4*vs**2)/(vp**2-vs**2)
    >>> mu = (vp**2-2*vs**2)/(2*(vp**2-vs**2))
    >>> H = [140/1000]
    >>> rmax = 70/2/1000
    >>> BE = one_phase_dynamic(tmin, tmax, tlin, E, mu, rho, H, rmax/100, rmax, 1, 3, 0, 50)
    >>> yt = BE.run_f() 



Ultrasonic pulse velocity
======================
Ultrasoni cpulse velocity (UPV) test is commonly used for anomaly detection and strength evaluationof construction materials 
(e.g., concrete and steel). For the demonsration purpose, the example code for the  pile  integrity  test is given as follows: 


    >>> import numpy as np  
    >>> from geondt import one_phase_dynamic    
    >>> def load(s):  
        ''' Define external load with 50 kHz in the Laplace domain (for pile-soil interaction example)'''
        fn1 = -100*10**3*(np.exp(-np.complex(s)/(10*10**3)))*np.pi/(10*10.0**9 *np.pi**2 + np.complex(s)**2)  
        fn2 = 100*10**3*(np.exp(-np.complex(s)/(12.5*10**3)))*np.pi/(10*10.0**9 *np.pi**2 + np.complex(s)**2)  
        fn = fn1 + fn2      
        return fn  
    >>> tlin = 500; tmin =2*10**(-5); tmax = 100*10**(-5)      
    >>> vp = np.array([3500.0**2, 3500.0**2, 100.0**2])
    >>> vs = np.array([2000.0**2, 2000.0**2, (100/1.5)**2])
    >>> rho  = np.array([2400, 2400, 1800])
    >>> E = rho1*vs*(3*vp-4*vs)/(vp-vs)
    >>> mu = (vp**2-2*vs**2)/(2*(vp**2-vs**2))
    >>> H = [0.4,0.1,10]
    >>> rmax = 0.5
    >>> Pile = one_phase_dynamic(tmin, tmax, tlin, E, mu, rho, H, rmax/100, rmax, 4, 2, 0, 50) 
    >>> Pile.load_i = load
    >>> yt = Pile.run_f() 


Falling Weight Deflectometer
======================

Falling  weight  deflectometer  (FWD)  is  anotherin-situ testing method used to evaluate the mechanical
properties of pavement structures. The FWD test for a three-layer pavement system can be studied by the GeoNDT: 

    >>> import numpy as np  
    >>> from geondt import one_phase_dynamic  
    >>> tmin = 0.001; tmax=0.05;  tlin=200; rmax = 20 
    >>> E =  [1000.0*(10.0**6), 200.0*(10.0**6), 100.0*(10.0**6)] # Young ’s modulus  in each  layer
    >>> mu = [0.35, 0.35, 0.35]  # Poisson ’s ratio in each  layer
    >>> rho1 = [2300.0, 2000.0,1500.0] # density  in each  layer
    >>> H1 = [0.15,0.25,5.0]  # thickness  in each  layer
    >>> FWD= one_phase_dynamic(tmin, tmax, tlin, E, mu, rho1, H1, 0, rmax, 2, 2, 0)   
    >>> Displacement1 = FWD.run_i()  

Multichannel Analysis of Surface Waves
====================== 

Multichannel Analysis of Surface Waves (MASW) is one of the most popular techniquesfor the in-situ 
evaluation of the shear wave velocity in different soil layers.  The code example for MASW modeling
for a three-layer system is shown as follows: 

    >>> fmin = 5; fmax = 50; flin = 50 
    >>> omega = np.linspace(fmin,fmax,flin)*np.pi*2
    >>> E =  [100.0*(10.0**6), 200.0*(10.0**6), 500.0*(10.0**6)]
    >>> mu = [0.2, 0.3, 0.35] 
    >>> rho = [2300.0, 2000.0,2500.0]
    >>> H = [2,3,5.0]
    >>> MASW = one_phase_dispersion(fmin,fmax,flin, E, mu, rho, H)  
    >>> yt = MASW.run()