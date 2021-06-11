"""
This module is used to perform dynamic analysis for finite and infinite domain in one-phase materials. 
"""
from dataclasses import dataclass    
from math import erfc  
from joblib import Parallel, delayed  
import json  
import numpy as np 
from geondt import PoroSEM
from scipy.special import jv,jn_zeros
from scipy import special,optimize  
from .inverselaplace import ilt

@dataclass
class one_phase_dynamic:
    ''' 
    Object for wave propagation modeling of one-phase material in finite and infinite space.

    :param tmin: A mandatory (floating) parameter that defines the minimum time (s)
    :param tmax: A mandatory (floating) parameter that defines the maximum time (s)
    :param tlin: A mandatory (integer) parameter that defines the number of points within the maximum and minimum time
    :param E:    A mandatory (floating) list that defines the Young's modulus (Pa) for all layers  
    :param mu:   A mandatory (floating) list that defines the Poisson's ratio for all layers     
    :param rho:  A mandatory (floating) list that defines the density (kg/m^3) for all layers
    :param H:    A mandatory (floating) list that defines the thickness (m) for all layers
    :param r:    A mandatory (floating) parameter that defines the radial distance (m) of the location of interest 
    :param rmax: A mandatory (floating) parameter that defines the maximum radius (m) of the domain 
    :param node: A mandatory (integer) parameter that spesfies the ID of the output node (for instance, the 'node' for surface vertical displacement is 2)
    :param loc:  An optioanl (integer) parameter that defines the ID of the input load (for instance, the 'loc' for surface vertical load is 2)
    :param st:   An optional (integer) parameter (either 0 or 1) that defined the output type, where 0 represnts displacement and 1 represents stress output
    :param m_num: An optional (integer) parameter that defines the number of modes in the Fourier-Bessel series
    :param lap_num: An optional (integer) parameter that defines the number of iteration in the inverse laplace transform
    :param a1: An optional (integer) parameter that defines the radius (m) of the contact area between external load and domain  
    :param ncore: An optional (integer) parameter that defines the number of cores (-1 is used in default to use all available cores)
    :return: dynamic response in terms of displacement or stress
    ''' 
    tmin: float 
    tmax: float 
    tlin: int 
    E: float
    mu: float 
    rho: float
    H: float
    r: float
    rmax: float  
    node: int 
    loc: int = 2 
    st: int = 0  
    m_num: int = 300 
    lap_num: int = 35 
    ncore: int = -1
    a1: float = 0.15 


    def model_f(self, s):  
        ''' Define one-phase model for finite domain'''
        km = jn_zeros(0,self.m_num)/self.rmax
        fm=np.ones(self.m_num,dtype = float) # point load
        fn = self.load_f(s)  
        yt =  PoroSEM.one_phase_finite(np.complex(s), self.E, self.mu, self.rho, self.H, km, fm, self.r, np.complex(fn), self.node, self.loc, self.st, self.m_num, len(self.E))  
        return yt 
    
    def model_i(self, s): 
        ''' Define one-phase model for infinite domain'''
        km = jn_zeros(0,self.m_num)/self.rmax
        fm=np.ones(self.m_num,dtype = float) # point load
        for i in range(self.m_num): 
            fm[i] = 2*self.a1*jv(1,km[i]*self.a1)/(self.rmax**2*km[i]*jv(1,self.rmax*km[i])**2) #*(m_num + 1-i)/(m_num+1)*1.05*1.02 
           
        fn = self.load_i(s)   
        yt =  PoroSEM.one_phase_infinite(np.complex(s), self.E, self.mu, self.rho, self.H, km, fm, self.r, np.complex(fn), self.node, self.loc, self.st, self.m_num, len(self.E))  
        return yt
 
    def load_f(self, s):  
        ''' Define external load in the Laplace domain (for BE example)'''
        # fn1 = -100*10**3*(np.exp(-np.complex(s)/(10*10**3)))*np.pi/(10*10.0**9 *np.pi**2 + np.complex(s)**2)  
        # fn2 = 100*10**3*(np.exp(-np.complex(s)/(12.5*10**3)))*np.pi/(10*10.0**9 *np.pi**2 + np.complex(s)**2)  
        # fn = fn1 + fn2   
        fn1 = -20*10.0**3*(np.exp(-np.complex(s)/(2000)))*np.pi/(400*10.0**6*np.pi**2 + np.complex(s)**2)  
        fn2 = 20*10.0**3*(np.exp(-np.complex(s)/(2500)))*np.pi/(400*10.0**6*np.pi**2 + np.complex(s)**2)  
        fn = fn1 + fn2      
        return fn      

    def load_i(self, s):  
        ''' Define external load in the Laplace domain (for FWD example)'''
        fn = 3200.0*(1-np.exp(-np.complex(s)/40.0))*np.pi**2/(6400.0*np.pi**2*np.complex(s)+np.complex(s)**3)*50000   
        return fn     

    def inv_f(self, i):  
        ''' Define function to be used in the Parallel computing for finite domain'''
        t = np.linspace(self.tmin, self.tmax, self.tlin)  
        y = ilt(self.model_f,t[i], self.lap_num, "cme") 
        # y = invertlaplace(aimf,t[i],degree = 20)  #  option 2 
        return y        
    
    def run_f(self):  
        ''' Run the model using the  Parallel computing for finite domain''' 
        t = np.linspace(self.tmin, self.tmax, self.tlin)  
        yt = Parallel(n_jobs=self.ncore)(delayed(self.inv_f)(i) for i in range(len(t))) 
        return np.real(yt)  

    def inv_i(self, i):  
        ''' Define function to be used in the Parallel computing for infinite domain'''
        t = np.linspace(self.tmin, self.tmax, self.tlin)  
        y = ilt(self.model_i,t[i], self.lap_num, "cme") 
        # y = invertlaplace(aimf,t[i],degree = 20)  #  option 2 
        return y        
    
    def run_i(self):  
        ''' Run the model using the  Parallel computing for half-space domain'''
        t = np.linspace(self.tmin, self.tmax, self.tlin)  
        yt = Parallel(n_jobs=self.ncore)(delayed(self.inv_i)(i) for i in range(len(t))) 
        return np.real(yt)  

