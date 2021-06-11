"""
This module is used to perform dispersion analysis for infinite domain in one-phase materials. 
"""
from dataclasses import dataclass    
from math import erfc  
from joblib import Parallel, delayed  
import json  
import numpy as np 
from geondt import PoroSEM
from scipy.special import jv,jn_zeros
from scipy import special,optimize 

@dataclass
class one_phase_dispersion:
    ''' 
    Object for the dispersion analysis of one-phase material with the last layer as infinite half space.

    :param f1: A mandatory (floating) parameter that defines the lower bound of frequency (Hz)
    :param f2: A mandatory (floating) parameter that defines the upper bound frequency (Hz)
    :param flin: A mandatory (integer) parameter that defines the number of points within the lower and upper bound of frequency
    :param E:    A mandatory (floating) list that defines the Young's modulus (Pa) for all layers  
    :param mu:   A mandatory (floating) list that defines the Poisson's ratio for all layers     
    :param rho:  A mandatory (floating) list that defines the density (kg/m^3) for all layers
    :param H:    A mandatory (floating) list that defines the thickness (m) for all layers
    :param ncore: An optional (integer) parameter that defines the number of cores (-1 is used in default to use all available cores)
    :return:  Dispersion relations
    ''' 
    f1: float 
    f2: float  
    flin: int  
    E: float
    mu: float 
    rho: float
    H: float  
    ncore: int = -1

    def model_d(self, i): 
        ''' Define dispersion modelling and root searching algorithm'''
        def fun(k): 
            matrix =PoroSEM.one_phase_dispersion(i, self.E, self.mu, self.rho, self.H, k, len(self.E))
            matrix1 = np.array(matrix)
            matrix2 = np.real(matrix1)/10**11
            sign, logdet = np.linalg.slogdet(matrix2)
            return sign * np.exp(logdet)
        c_test = 10
        incre = i / c_test
        root = 0.00001
        for j in range(10**6):
            past = incre
            val1 = fun(incre)
            incre =  incre - 0.01
            val2 = fun(incre)
            if (np.real(val1) * np.real(val2) <= 0):
                root =  optimize.brentq(fun,incre, past)      
                break 
        return (i/root)     #give one value at a frequency 

    def final(self, n): 
        ''' Define function to be used for the Parallel computing'''
        omega = np.linspace(self.f1, self.f2,self.flin)*np.pi*2
        y =  self.model_d(omega[n])
        return y    
 
    def run(self):  
        ''' Run the model using the Parallel computing for infinite domain'''
        omega = np.linspace(self.f1, self.f2,self.flin)*np.pi*2  
        yt = Parallel(n_jobs=self.ncore)(delayed(self.final)(i) for i in range(len(omega))) 
        return np.real(yt)  
 
