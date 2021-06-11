==============
References
==============

Considering the computational efficiency and user-friendly environment, the GeoNDT software is developed using a dual layer/hybrid Python and Fortran environment 
to benefit from the strengths of the two languages, as follows: a) Fortran is a compiled language; it is closer to the material architecture it is executed on; 
it benefits from established mathematical libraries and thus compensates for the lower computational performance of interpreted languages such as Python under 
CPU intensive tasks; b)  Python is a user-friendly language; it has a wide support online; it has a rich set of high-quality scientific computational libraries 
and frameworks; and it offers improved code reusability, faster and cost-effective development. The detailed discription of the functions in GeoNDT can be found 
in following reference: 


one_phase_dynamic
====================== 

one_phase_dynamic function is used to study the wave propagation and dynamic response in one-phase material. The inputs of this function 
is listed as follows:   

       :parameter tmin: A mandatory (floating) parameter that defines the minimum time (s)
       :parameter tmax: A mandatory (floating) parameter that defines the maximum time (s)
       :parameter tlin: A mandatory (integer) parameter that defines the number of points within the maximum and minimum time
       :parameter E:    A mandatory (floating) list that defines the Young's modulus (Pa) for all layers  
       :parameter mu:   A mandatory (floating) list that defines the Poisson's ratio for all layers     
       :parameter rho:  A mandatory (floating) list that defines the density (kg/m^3) for all layers
       :parameter H:    A mandatory (floating) list that defines the thickness (m) for all layers
       :parameter r:    A mandatory (floating) parameter that defines the radius (m) of the location of interest 
       :parameter rmax: A mandatory (floating) parameter that defines the maximum radius (m) of the domain 
       :parameter node: A mandatory (integer) parameter that spesfies the ID of the output node (for instance, the 'node' for surface vertical displacement is 2)
       :parameter loc:  An optioanl (integer) parameter that defines the ID of the input load (for instance, the 'loc' for surface vertical load is 2)
       :parameter st:   An optional (integer) parameter (either 0 or 1) that defined the output type, where 0 represnts displacement and 1 represents stress output
       :parameter m_num: An optional (integer) parameter that defines the number of modes in the Fourier-Bessel series
       :parameter lap_num: An optional (integer) parameter that defines the number of iteration in the inverse Laplace transform
       :parameter a1: An optional (integer) parameter that defines (m) the radius on the contact area between external load and domain  
       :parameter ncore: An optional (integer) parameter that defines the number of cores (-1 is used in default to use all available cores)

 
one_phase_dispersion
====================== 

one_phase_dynamic function is used to perform dispersion analysis of one-phase material 
with the last layer as infinite half space.  The inputs of this function 
is listed as follows:   

       :parameter f1: A mandatory (floating) parameter that defines the lower bound of frequency (Hz)
       :parameter f2: A mandatory (floating) parameter that defines the upper bound frequency (Hz)
       :parameter flin: A mandatory (integer) parameter that defines the number of points within the lower and upper bound of frequency
       :parameter E:    A mandatory (floating) list that defines the Young's modulus (Pa) for all layers  
       :parameter mu:   A mandatory (floating) list that defines the Poisson's ratio for all layers     
       :parameter rho:  A mandatory (floating) list that defines the density (kg/m^3) for all layers
       :parameter H:    A mandatory (floating) list that defines the thickness (m) for all layers
       :parameter ncore: An optional (integer) parameter that defines the number of cores (-1 is used in default to use all available cores)


two_phase_dynamic
======================  

Object for wave propagation modeling of two-phase material in finite and infinite space.  The inputs of this function 
is listed as follows:   

    :parameter tmin: A mandatory (floating) parameter that defines the minimum time (s)
    :parameter tmax: A mandatory (floating) parameter that defines the maximum time (s)
    :parameter tlin: A mandatory (integer) parameter that defines the number of points within the maximum and mimumun time
    :parameter E:    A mandatory (floating) list that defines the Young's modulus (Pa) of bulk soil for all layers 
    :parameter mu:   A mandatory (floating) list that defines the Poisson's ratio of bulk soil for all layers  
    :parameter rho:  A mandatory (floating) list that defines the density (kg/m^3) of bulk soil for all layers
    :parameter H:    A mandatory (floating) list that defines the thickness (m) for all layers 
    :parameter kh:    A mandatory (floating) list that defines the permeability coefcient (m^2) of bulk soil for all layers 
    :parameter porosity: A mandatory (floating) list that defines the porosity of bulk soil for all layers 
    :parameter r:    A mandatory (floating) parameter that defines the radius (m) of the location of interest 
    :parameter rmax: A mandatory (floating) parameter that defines the maximum radius (m) of the domain 
    :parameter node: A mandatory (integer) parameter that spesfies the ID of the output node (for instance, the 'node' for surface vertical displacement is 2)
    :parameter loc:  An optioanl (integer) parameter that defines the ID of the input load (for instance, the 'loc' for surface vertical load is 2)
    :parameter st:   An optional (integer) parameter (either 0 or 1) that defined the output type, where 0 represnts displacement and 1 represents stress output
    :parameter vis: An optional (floating) parameter that defines the viscosity (Pa-s) of water
    :parameter kf: An optional (floating) parameter that defines the bulk modulus (Pa) of the fluid
    :parameter m_num: An optional (integer) parameter that defines the number of modes in the Fourier-Bessel series
    :parameter lap_num: An optional (integer) parameter that defines the number of iteration in the inverse Laplace transform
    :parameter a1: An optional (integer) parameter that defines the radius (m) of the contact area between external load and domain  


two_phase_dispersion
====================== 

Object for one-phase material in infinite space dispersion modeling.  The inputs of this function 
is listed as follows:   

    :parameter f1: A mandatory (floating) parameter that defines the lower bound of frequency (Hz)
    :parameter f2: A mandatory (floating) parameter that defines the upper bound frequency (Hz)
    :parameter flin: A mandatory (integer) parameter that defines the number of points within the lower and upper bpund of frequency
    :parameter E:    A mandatory (floating) list that defines the Young's modulus (Pa) of bulk soil for all layers 
    :parameter mu:   A mandatory (floating) list that defines the Poisson's ratio of bulk soil for all layers  
    :parameter rho:  A mandatory (floating) list that defines the density (kg/m^3) of bulk soil for all layers
    :parameter H:    A mandatory (floating) list that defines the thickness (m) for all layers 
    :parameter kh:    A mandatory (floating) list that defines the permeability coefcient (m^2) of bulk soil for all layers 
    :parameter porosity: A mandatory (floating) list that defines the porosity of bulk soil for all layers 
    :parameter vis: An optional (floating) parameter that defines the viscosity (Pa-s) of water
    :parameter kf: An optional (floating) parameter that defines the bulk modulus (Pa) of the fluid
    :parameter ncore: An optional (integer) parameter that defines the number of cores (-1 is used in default to use all available cores) 


three_phase_dynamic
======================  

Object for wave propagation modeling of three-phase material in finite and infinite space.  The inputs of this function 
is listed as follows:   

    :parameter tmin: A mandatory (floating) parameter that defines the minimum time (s)
    :parameter tmax: A mandatory (floating) parameter that defines the maximum time (s)
    :parameter tlin: A mandatory (integer) parameter that defines the number of points within the maximum and mimumun time
    :parameter kks:  A mandatory (floating) list that defines the bulk modulus of soild skeleton (Pa) for all layers  
    :parameter muus: A mandatory (floating) list that defines the shear modulus  of soild skeleton (Pa) for all layers    
    :parameter rho:  A mandatory (floating) list that defines the density of soild skeleton (kg/m^3) for all layers   
    :parameter phiw:  A mandatory (floating) list that defines the volumetric water content for all layers  
    :parameter phii:  A mandatory (floating) list that defines the volumetric ice content for all layers  
    :parameter H:    A mandatory (floating) list that defines the thickness (m) for all layers  
    :parameter r:    A mandatory (floating) parameter that defines the radius of the location of interest 
    :parameter rmax: A mandatory (floating) parameter that defines the maximum radius of the domain 
    :parameter node: A mandatory (integer) parameter that spesfies the ID of the output node (for surface vertical displacement, it is 2)
    :parameter loc:  A mandatory (integer) parameter that defines the ID of the input load
    :parameter st:   A mandatory (integer) parameter (either 0 or 1) that defined the output type where 0 represnts displacement and 1 represents stress
    :parameter KKi:  An optional (floating) parameter that defines the bulk modulus of ice (Pa)  
    :parameter muui:  An optional (floating) parameter that defines the shear modulus of ice (Pa) 
    :parameter Kww:  An optional (floating) parameter that defines the bulk modulus of water (Pa)         
    :parameter kappas:  An optional (floating) parameter related to the peameability of solid skeleton
    :parameter kappai:  An optional (floating) parameter  to the peameability of ice 
    :parameter b013:  An optional (floating) parameter that defines friction coefficient between the solid skeletal frame and ice matrix      
    :parameter m_num: A optional (integer) parameter that defines the number of modes 
    :parameter lap_num: A optional (integer) parameter that defines the number of iteration in the inverse Laplace transform
    :parameter a1: A optional (integer) parameter that defines the radius (m) of the contact area between external load and domain  
    :parameter ncore: A optional (integer) parameter that defines the number of cores (-1 is used in default to use all available cores) 



three_phase_dispersion
====================== 
Object for one-phase material in infinite space dispersion modeling. The inputs of this function 
is listed as follows:   

    :parameter f1: A mandatory (floating) parameter that defines the lower bound of frequency (Hz)
    :parameter f2: A mandatory (floating) parameter that defines the upper bound frequency (Hz)
    :parameter flin: A mandatory (integer) parameter that defines the number of points within the lower and upper bpund of frequency
    :parameter kks:  A mandatory (floating) list that defines the bulk modulus of soild skeleton (Pa) for all layers  
    :parameter muus: A mandatory (floating) list that defines the shear modulus  of soild skeleton (Pa) for all layers    
    :parameter rho:  A mandatory (floating) list that defines the density of soild skeleton (kg/m^3) for all layers   
    :parameter phiw:  A mandatory (floating) list that defines the volumetric water content for all layers  
    :parameter phii:  A mandatory (floating) list that defines the volumetric ice content for all layers  
    :parameter H:    A mandatory (floating) list that defines the thickness (m) for all layers  
    :parameter KKi:  An optional (floating) parameter that defines the bulk modulus of ice (Pa)  
    :parameter muui:  An optional (floating) parameter that defines the shear modulus of ice (Pa) 
    :parameter Kww:  An optional (floating) parameter that defines the bulk modulus of water (Pa)         
    :parameter kappas:  An optional (floating) parameter related to the peameability of solid skeleton
    :parameter kappai:  An optional (floating) parameter  to the peameability of ice 
    :parameter b013:  An optional (floating) parameter that defines friction coefficient between the solid skeletal frame and ice matrix  
    :parameter ncore: An optional (integer) parameter that defines the number of cores (-1 is used in default to use all available cores) 



functions inside dynamic modules 
====================== 

model_f deefine the dynamic model for finite domain. For instance::

    def model_f(self, s):  
        ''' Define one-phase model for finite domain'''
        km = jn_zeros(0,self.m_num)/self.rmax
        fm=np.ones(self.m_num,dtype = float) # point load
        fn = self.external_load(s)  
        yt =  PoroSEM.one_phase_finite(np.complex(s), self.E, self.mu, self.rho, self.H, km, fm, self.r, np.complex(fn), self.node, self.loc, self.st, self.m_num, len(self.E))  
        return yt 

model_i deefine the dynamic model for infinite domain. For instance::  
    
    def model_i(self, s): 
        ''' Define one-phase model for infinite domain'''
        km = jn_zeros(0,self.m_num)/self.rmax
        fm=np.ones(self.m_num,dtype = float) # point load
        for i in range(self.m_num): 
            fm[i] = 2*self.a1*jv(1,km[i]*self.a1)/(self.rmax**2*km[i]*jv(1,self.rmax*km[i])**2) #*(m_num + 1-i)/(m_num+1)*1.05*1.02 
        fn = self.external_load2(s)   
        yt =  PoroSEM.one_phase_infinite(np.complex(s), self.E, self.mu, self.rho, self.H, km, fm, self.r, np.complex(fn), self.node, self.loc, self.st, self.m_num, len(self.E))  
        return yt 

In above function, the external load is required to be defined in the Lapalce domain.  For example, 
the 10 kHz load in the BE test can be defined as::

    def external_load(self, s):  
        ''' Define external load in the Lapalce domain (for BE example)''' 
        fn1 = -20*10.0**3*(np.exp(-np.complex(s)/(2000)))*np.pi/(400*10.0**6*np.pi**2 + np.complex(s)**2)  
        fn2 = 20*10.0**3*(np.exp(-np.complex(s)/(2500)))*np.pi/(400*10.0**6*np.pi**2 + np.complex(s)**2)  
        fn = fn1 + fn2      
        return fn  

To run the dynamic response in the finite domain, users should use 'run_f' function, defined as:: 

    def run_f(self):  
        ''' Run the model using the  Parallel computing for finite domain''' 
        t = np.linspace(self.tmin, self.tmax, self.tlin)  
        yt = Parallel(n_jobs=self.ncore)(delayed(self.inv_f)(i) for i in range(len(t))) 
        return np.real(yt)  

To run the dynamic response in the infinite domain, users should use 'run_i' function, defined as::  

    def run_i(self):  
        ''' Run the model using the  Parallel computing for infinite domain'''
        t = np.linspace(self.tmin, self.tmax, self.tlin)  
        yt = Parallel(n_jobs=self.ncore)(delayed(self.inv_i)(i) for i in range(len(t))) 
        return np.real(yt)   


For dispersion analysis, users can call 'run' function directly, defined as follows:: 

    def run(self):  
        ''' Run the model using the Parallel computing for infinite domain'''
        omega = np.linspace(self.f1, self.f2,self.flin)*np.pi*2  
        yt = Parallel(n_jobs=self.ncore)(delayed(self.final)(i) for i in range(len(omega))) 
        return np.real(yt)  