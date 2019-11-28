import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

#%%
np.set_printoptions(precision=3)

#%% CG single, first way

class tia_power_consumption:
    def __init__(self):
        self.Id_1 = 0
        self.Id_2 = 0
        self.Id_3 = 0
        self.R = 0
        self.I_ref = 0
        
    def set_R(R):
        self.R = R
        
    def set_I_ref(I_ref):
        self.I_ref
        
    def set_Ids(Id_1, Id_2, Id_3):
        self.Id_1 = Id_1
        self.Id_2 = Id_2
        self.Id_3 = Id_3
        self.upd_power()
    
    def upd_power(self):
        self.power1 = 2*5*(self.Id_1 + self.Id_2 + self.Id_3) + 5*self.I_ref + 50/(4*self.R)
        self.power2 = 2*5*(self.Id_1 + self.Id_2 + self.Id_3) 
        self.power3 = 5*self.I_ref
        self.power4 = 50/(4*self.R_set)
        
    def get_power(self):
        return self.power1
    
    def print_power(self, n):
        if n >= 0:
            print('power consumption: {self.power1}')
        if n == 0:
            print('power consumption Ids:   {self.power2}')
            print('power consumption I_ref: {self.power3}')
            print('power consumption R:     {self.power4}')

def unit_test_tia_power_consumption():
    tia_power_module_test = tia_power_consumption()
    tia_power_module_test.set_Ids(1.2, 5E-5)
    
class tia:
    
    unCox = 50E-6
    upCox = 25E-6
    # lookup table for l=1um nmos
    wl_1u    = [2,     5,     10,    20,    50,    100,  200,  500,  1000, 2000]
    cdtot_1u = [5.60,  9.50,  16.0,  29.0,  68.0,  134,  263,  653,  1303, 2603]
    cgtot_1u = [3.45,  8.60,  17.2,  34.5,  86.3,  173,  345,  863,  1727, 3453]
    cstot_1u = [5.62,  9.50,  16.1,  29.2,  68.5,  133,  265,  658,  1313, 2623]
    cbtot_1u = [10.6,  17.6,  29.2,  52.3,  122,   238,  469,  1164, 2322, 4640]
    cgs_1u   = [1.02,  2.55,  5.10,  10.2,  25.5,  50,   102,  255,  510,  1020]
    cgd_1u   = [1.00,  2.50,  5.00,  10.0,  25.0,  51,   100,  250,  500,  1000]
    cdtot_lookup_1u = interpolate.interp1d(wl_1u, cdtot_1u)
    cgtot_lookup_1u = interpolate.interp1d(wl_1u, cgtot_1u)
    cstot_lookup_1u = interpolate.interp1d(wl_1u, cstot_1u)
    cbtot_lookup_1u = interpolate.interp1d(wl_1u, cbtot_1u)
    cgs_lookup_1u   = interpolate.interp1d(wl_1u, cgs_1u  )
    cgd_lookup_1u   = interpolate.interp1d(wl_1u, cgd_1u  )
    
    # lookup table for l=1um
    w_2u     = [2,     5,     10,    20,    50,    100,   200,   500]
    wl_2u    = [1,     2.5,   5,     10,    25,    50,    100,   250]
    cdtot_2u = [5.60,  9.50,  16.0,  29.0,  68.0,  133,   263,   653]
    cgtot_2u = [4.91,  12.3,  24.5,  49.1,  123,   245,   490,   1227]
    cstot_2u = [5.64,  9.60,  16.2,  29.4,  68.9,  135,   267,   663]
    cbtot_2u = [12.1,  21.2,  36.3,  66.7,  158,   310,   613,   1523]
    cgs_2u   = [1.04,  2.60,  5.20,  10.4,  26,    52,    103,   260]
    cgd_2u   = [1.00,  2.50,  5.00,  10.0,  25.0,  50,    100,   250]
    cdtot_lookup_2u = interpolate.interp1d(wl_2u, cdtot_2u)
    cgtot_lookup_2u = interpolate.interp1d(wl_2u, cgtot_2u)
    cstot_lookup_2u = interpolate.interp1d(wl_2u, cstot_2u)
    cbtot_lookup_2u = interpolate.interp1d(wl_2u, cbtot_2u)
    cgs_lookup_2u   = interpolate.interp1d(wl_2u, cgs_2u  )
    cgd_lookup_2u   = interpolate.interp1d(wl_2u, cgd_2u  )
    
    R_LDM = 10000 # ohms
    C_LDM = 500   # fF
    C_in  = 100   # fF
    
    def __init__(self):
        # configuration variables
        self.Vov_1 = 0.2
        self.Vov_2 = 0.2
        self.A_cs  = 3
        self.Vov_3 = 0.2
        self.Vov_B = 0.4
        self.Vov_L = 0.6
        self.Id_1 = 3E-5
        self.Id_2 = 6E-5
        self.Id_3 = 3E-5
        
        
        self.I_ref= 2E-5
        self.power1 = -1
        self.power2 = -1
        self.power3 = -1
        self.power4 = -1
        
        self.gm1 = 0
        self.gm2 = 0
        self.gm3 = 0
        self.R_set = 20000
        self.ro2= 0
        self.R1 = 0
        self.R2 = 0
        self.R3 = 0
        self.R4 = 0
        self.A1 = 0
        self.A2 = 0
        
        self.WL_1  = 100
        self.WL_2  = 100
        self.WL_3  = 100
        
        self.WL_L1 = 4
        self.WL_L2 = 4
        self.WL_B1 = 2
        self.WL_B2 = 2
        self.WL_B3 = 2
        
        self.WL_i1 = 2
        self.WL_i2 = 2
        self.WL_i3 = 2
        
        self.tran_res = 0
        
        self.C1 = 0
        self.C2 = 0
        self.C3 = 0
        self.C4 = 0
        
        self.b11 = 0
        self.b12 = 0
        self.b13 = 0
        self.b14 = 0
        self.b1  = 0
        self.w3dB = 0
        self.f3dB_ZVTC = 0
        
        self.n = 50
        self.f_sweep = np.logspace( 5, 9, self.n )
        
        self.mag_v_out  = np.zeros(self.n)
        self.mag_in     = np.zeros(self.n)
        self.mag_cg_out = np.zeros(self.n)
        self.mag_cs_out = np.zeros(self.n)
        
        self.mag_v_out_lookup  = np.zeros(self.n)
        self.mag_in_lookup     = np.zeros(self.n)
        self.mag_cg_out_lookup = np.zeros(self.n)
        self.mag_cs_out_lookup = np.zeros(self.n)
        
        self.f3dB_v_out  = -1
        self.f3dB_in     = -1
        self.f3dB_cg_out = -1
        self.f3dB_cs_out = -1
        
        self.FOM = 0
        
    def print_all(self):
        print('Vov_1: %1.2f  v' %self.Vov_1)
        print('Vov_2: %1.2f  v' %self.Vov_2)
        print('Vov_3: %1.2f  v' %self.Vov_3)
        print('Vov_B: %1.2f  v' %self.Vov_B)
        print('Vov_L: %1.2f  v' %self.Vov_L)
        print()
        
        print('Id_1:  %3.1f uA' %(1E+6*self.Id_1))
        print('Id_2:  %3.1f uA' %(1E+6*self.Id_2))
        print('Id_3:  %3.1f uA' %(1E+6*self.Id_3))
        print('I_ref: %3.1f uA' %(1E+6*self.I_ref))
        print()
        
        print('gm1:   %3.2f   mA/V' %(1E+3*self.gm1))
        print('gm2:   %3.2f   mA/V' %(1E+3*self.gm2))
        print('gm3:   %3.2f   mA/V' %(1E+3*self.gm3))
        print()
        
        print('R1:    %2.2f   kohms' %(1E-3*self.R1))
        print('R_set: %2.2f   kohms' %(1E-3*self.R_set))
        print('ro2:   %3.2f   kohms' %(1E-3*self.ro2))
        print('R2:    %2.2f   kohms' %(1E-3*self.R2))
        
        print('R3:    %2.2f   kohms' %(1E-3*self.R3))
        print('R4:    %2.2f   kohms' %(1E-3*self.R4))
        print()
        
        print('WL_L1: %3.1f,  w=%3.1f' %(self.WL_L1, 2*self.WL_L1))
        print('WL_1:  %3.1f' %(self.WL_1) )
        print('WL_B1: %3.1f,  w=%3.1f' %(self.WL_B1, 2*self.WL_B1))
        print()
        
        print('WL_L2: %3.1f,  w=%3.1f' %(self.WL_L2, 2*self.WL_L2))
        print('WL_2:  %3.1f' %(self.WL_2) )
        print('WL_B2: %3.1f,  w=%3.1f' %(self.WL_B2, 2*self.WL_B2))
        print()
        
        print('WL_3:  %3.1f' %(self.WL_3) )
        print('WL_B3: %3.1f,  w=%3.1f' %(self.WL_B3, 2*self.WL_B3))
        print()
        
        print('WL_i1: %3.1f,  w=%3.1f' %(self.WL_i1, 2*self.WL_i1))
        print('WL_i2: %3.1f,  w=%3.1f' %(self.WL_i2, 2*self.WL_i2))
        print('WL_i3: %3.1f,  w=%3.1f' %(self.WL_i3, 2*self.WL_i3))
        print()
        
        print('A1: %2.2f' %(self.A1) )
        print('check1 A1: %2.2f' %( np.sqrt((self.WL_2) / self.WL_L2) ) )
        print('check2 A1: %2.2f' %( self.Vov_L / (self.Vov_2) ) )
        
        print('R1:    %2.2f   kohms' %(1E-3*self.R1))
        print('R_set: %2.2f   kohms' %(1E-3*self.R_set))
        print('ro2:   %3.2f   kohms' %(1E-3*self.ro2))
        print('R2:    %2.2f   kohms' %(1E-3*self.R2))
        
        print('R3:    %2.2f   kohms' %(1E-3*self.R3))
        print('R4:    %2.2f   kohms' %(1E-3*self.R4))
        print()
        
        print('A2: %2.2f' %(self.A2))
        print('C1: %4.1f fF' %(1E+15*self.C1) )
        print('C2: %4.1f fF' %(1E+15*self.C2) )
        print('C3: %4.1f fF' %(1E+15*self.C3) )
        print('C4: %4.1f fF' %(1E+15*self.C4) )
        
        print('b11: %3.2f ns' %(1E+9*self.b11) )
        print('b12: %3.2f ns' %(1E+9*self.b12) )
        print('b13: %3.2f ns' %(1E+9*self.b13) )
        print('b14: %3.2f ns' %(1E+9*self.b14) )
        print('b1:  %3.2f ns' %(1E+9*self.b1 ) )
        
        print('f3dB_in:    %3.2f MHz' %(self.f3dB_in/1E+6) )
        print('f3dB_cg_out: %3.2f MHz' %(self.f3dB_cg_out/1E+6) )
        print('f3dB_cs_out: %3.2f MHz' %(self.f3dB_cs_out/1E+6) )
        print('f3dB_v_out: %3.2f MHz' %(self.f3dB_v_out/1E+6) )
        print('f3dB_ZVTC:  %3.2f MHz' %(self.f3dB_ZVTC/1E+6) )
        print()
        print('tran_res_gain: %2.2f kohms, %2.2f dB' %( self.gain/1000, 20*np.log10(self.gain) ) )
        print()
        print('power consumption total: %2.2f mW' %( 1000*self.power1 ) )
        print('power consumption Id: %2.2f mW' %( 1000*self.power2 ) )
        print('power consumption ref: %2.2f mW' %( 1000*self.power3 ) )
        print('power consumption resistors: %2.2f mW' %( 1000*self.power4 ) )
        
        print('FOM: %5.3f' %self.FOM)
        
    def upd(self):
        
        # update mosfets
        
        
        # calc gm's and R's
        self.gm1 = 2*self.Id_1 / self.Vov_1
        self.gm2 = 2*self.Id_2 / self.Vov_2
        self.gm3 = 2*self.Id_3 / self.Vov_3
        self.gmL2= 2*self.Id_2 / self.Vov_L
        self.ro2 = 1/(0.1*self.Id_2)
        self.R1 = (1/self.gm1)
        self.R2 = self.parallel(self.R_set, self.ro2)
        self.R3 = 1 / self.gmL2
        self.R4 = self.parallel( 0.84/self.gm3, self.R_LDM )
        
        # calc W/L = ( Id ) / ( 1/2*unCox*Vov^2 )
        self.WL_1  = (2*self.Id_1) / (self.unCox*self.Vov_1**2)
        self.WL_2  = (2*self.Id_2) / (self.unCox*self.Vov_2**2)
        self.WL_3  = (2*self.Id_3) / (self.unCox*self.Vov_3**2)
        
        self.WL_L1 = (2*self.Id_1) / (self.upCox*self.Vov_B**2)
        self.WL_L2 = (2*self.Id_2) / (self.unCox*(self.Vov_2*self.A_cs)**2)
        
        self.WL_B1 = (2*self.Id_1) / (self.unCox*self.Vov_B**2)
        self.WL_B2 = (2*self.Id_2) / (self.unCox*self.Vov_B**2)
        self.WL_B3 = (2*self.Id_3) / (self.unCox*self.Vov_B**2)
        
        self.WL_i1 = (2*self.I_ref) / (self.unCox*self.Vov_B**2)
        self.WL_i2 = (2*self.I_ref) / (self.unCox*self.Vov_B**2)
        self.WL_i3 = (2*self.I_ref) / (self.upCox*self.Vov_B**2)
        
        if self.WL_1 <= 2:
            print('WL_1 = %2.2f, below range' %self.WL_1)
        if self.WL_2 <= 2:
            print('WL_2 = %2.2f, below range' %self.WL_2)
        if self.WL_3 <= 2:
            print('WL_3 = %2.2f, below range' %self.WL_3)
            
        if self.WL_L1 <= 1:
            print('WL_L1 = %2.2f, below range' %self.WL_L1)
        if self.WL_L2 <= 1:
            print('WL_L2 = %2.2f, below range' %self.WL_L2)
       
        if self.WL_B1 <= 1:
            print('WL_B1 = %2.2f, below range' %self.WL_B1)
        if self.WL_B2 <= 1:
            print('WL_B2 = %2.2f, below range' %self.WL_B2)
        if self.WL_B3 <= 1:
            print('WL_B3 = %2.2f, below range' %self.WL_B3)
        
        
        self.A1 = ( self.WL_2 * self.Vov_2 ) / ( self.WL_L2 * self.Vov_L )
        self.A2 = -self.gm3 * self.parallel( 0.84/self.gm3, self.R_LDM )
        
        # calc C's
        self.C1 = 1E-15*( self.C_in 
                        + self.cstot_lookup_1u(self.WL_1) 
                        + self.cdtot_lookup_2u(self.WL_B1) )
        
        self.C2 = 1E-15*( self.cdtot_lookup_2u(self.WL_L1) 
                        + self.cdtot_lookup_1u(self.WL_1) 
                        + self.cgs_lookup_1u(self.WL_2) 
                        + (1 + self.A1)*self.cgd_lookup_1u(self.WL_2) )
        
        self.C3 = 1E-15*( (1 + 1/self.A1)*self.cgd_lookup_1u(self.WL_2)
                        + self.cdtot_lookup_1u(self.WL_2) - self.cgd_lookup_1u(self.WL_2)
                        + self.cstot_lookup_2u(self.WL_L2)
                        + self.cgd_lookup_1u(self.WL_3) 
                        + (1 + self.A2)*self.cgs_lookup_1u(self.WL_3) )
        
        self.C4 = 1E-15*( (1 + 1/self.A2)*self.cgs_lookup_1u(self.WL_3)
                        + self.cstot_lookup_1u(self.WL_3) - self.cgs_lookup_1u(self.WL_3) 
                        - self.cdtot_lookup_2u(self.WL_B3) 
                        + self.C_LDM )
        
        # calc b's
        self.b11 = self.R1 * self.C1
        self.b12 = self.R2 * self.C2
        self.b13 = self.R3 * self.C3
        self.b14 = self.R4 * self.C4
        self.b1  = self.b11 + self.b12 + self.b13 + self.b14
        self.f3dB_ZVTC = 1/(2*np.pi*self.b1)
        
        for i in range( np.size(self.f_sweep) ):
#            print(f'i: {i}')
            s = 2j*np.pi*self.f_sweep[i]
#            print(f'freq: {self.f_sweep[i]/1000000} MHz')
            
            self.mag_in[i]     = np.abs(  (1/self.b11) / (s + 1/self.b11) )
            self.mag_cg_out[i] = np.abs( ((1/self.b11)*(1/self.b12)) / ((s + 1/self.b11)*(s + 1/self.b12)) )
            self.mag_cs_out[i] = np.abs( ((1/self.b11)*(1/self.b12)*(1/self.b13)) / ((s + 1/self.b11)*(s + 1/self.b12)*(s + 1/self.b13)) )
            self.mag_v_out[i]  = np.abs( ((1/self.b11)*(1/self.b12)*(1/self.b13)*(1/self.b14)) / ((s + 1/self.b11)*(s + 1/self.b12)*(s + 1/self.b13)*(s + 1/self.b14)) )
#            print(f'mag_in:     {self.mag_in[i]}')
#            print(f'mag_cg_out: {self.mag_cg_out[i]}')
#            print(f'mag_cs_out: {self.mag_cs_out[i]}')
#            print(f'mag_v_out:  {self.mag_v_out[i]}')
        
        
#        plt.plot( (self.f_sweep ), self.mag_in ) 
#        plt.plot( (self.f_sweep ), self.mag_cg_out ) 
#        plt.plot( (self.f_sweep ), self.mag_cs_out ) 
#        plt.plot( (self.f_sweep ), self.mag_v_out ) 
                
        self.mag_in_lookup     = interpolate.interp1d( self.mag_in,     self.f_sweep )
        self.mag_cg_out_lookup = interpolate.interp1d( self.mag_cg_out, self.f_sweep )
        self.mag_cs_out_lookup = interpolate.interp1d( self.mag_cs_out, self.f_sweep )
        self.mag_v_out_lookup  = interpolate.interp1d( self.mag_v_out,  self.f_sweep )
        
        self.f3dB_in     = self.mag_in_lookup(0.71)
        self.f3dB_cg_out = self.mag_cg_out_lookup(0.71)
        self.f3dB_cs_out = self.mag_cs_out_lookup(0.71)
        self.f3dB_v_out  = self.mag_v_out_lookup(0.71)
        
        self.gain = -self.R2*self.A1*self.A2
        self.FOM = (self.gain*1E-3) * (self.f3dB_v_out*1E-6) / (self.power1*1E+3)
        
     
    def set_Vov_1(self, Vov_1):
        self.Vov_1 = Vov_1
    def set_Vov_2(self, Vov_2):
        self.Vov_2 = Vov_2
    def set_Vov_3(self, Vov_3):
        self.Vov_3 = Vov_3
    def set_Vov_L(self, Vov_L):
        self.Vov_L = Vov_L
    def set_Vov_B(self, Vov_B):
        self.Vov_B = Vov_B
        
    def set_Id_1(self, Id_1):
        self.Id_1 = Id_1
    def set_Id_2(self, Id_2):
        self.Id_2 = Id_2
    def set_Id_3(self, Id_3):
        self.Id_3 = Id_3
        
    def set_R_set(self, R_set):
        self.R_set = R_set
        
    def get_FOM(self):
        return self.FOM
    def get_gain(self):
        return self.gain
    def get_f3dB(self):
        return self.f3dB_v_out
    def get_power(self):
        return self.power1
        
        
    def parallel(self, R1, R2):
        return (R1*R2)/(R1+R2)
  
        
def unit_test_tia():
    tia1 = tia()
    
    tia1.set_Vov_1( 0.25 )
    tia1.set_Vov_2( 0.4 )
    tia1.set_Vov_3( 0.3 )
    tia1.set_Vov_L( 0.8 )
    tia1.set_Vov_B( 0.8 )
    tia1.set_Id_1( 3E-5 )
    tia1.set_Id_2( 2E-5 )
    tia1.set_Id_3( 5E-5 )
    tia1.set_R_set( 10E+3 )
    
    tia1.upd()
    tia1.print_all()

unit_test_tia()

def sweep_params():
    tia2 = tia()
    
    Vov_1 = np.linspace(0.25, 0.5, n)
    Vov_2 = np.linspace(0.25, 0.5, n)
    A_cs  = np.linspace(1, 3, 3)
    Vov_3 = 0.3
    Vov_L = 0.4
    Vov_B = 0.4
    Id_1  = np.linspace(1E-5, 1E-4, n)
    Id_2  = np.linspace(1E-5, 1E-4, n)
    Id_3  = np.linspace(1E-5, 1E-4, n)
    R_set = np.linspace(1E+3, 5E+4, n)
    
    """ sweep parameters with for loops
    put results in the pre-defined results matrix
    """
    print(f'total iterations: {n**3}')
    print('size Vov: {shape(Vov_1)}')
#    results = np.zeros(( n**9, 14))
#    count=-1
#    for ind_Vov_1 in range(np.size(Vov_1)):
#        for ind_Vov_2 in range(np.size(Vov_2)):
#            for ind_Vov_3 in range(np.size(Vov_3)):
#                
#                for ind_Vov_L in range(np.size(Vov_L)):
#                    for ind_Vov_B in range(np.size(Vov_B)):
#                        
#                        for ind_Id_1 in range(np.size(Id_1)):
#                            for ind_Id_2 in range(np.size(Id_2)):
#                                for ind_Id_3 in range(np.size(Id_3)):
#                                    
#                                    for ind_R_set in range(np.size(R_set)):
#                            
#                
#                                        count+=1
#                                        print(f'count: {count}')
#                                        tia2.set_Vov_1( Vov_1[ind_Vov_1] )
#                                        tia2.set_Vov_2( Vov_2[ind_Vov_2] )
#                                        tia2.set_Vov_3( Vov_3[ind_Vov_3] )
#                                        tia2.set_Vov_L( Vov_L[ind_Vov_L] )
#                                        tia2.set_Vov_B( Vov_B[ind_Vov_B] )
#                                        tia2.set_Id_1(   Id_1[ind_Id_1]  )
#                                        tia2.set_Id_2(   Id_2[ind_Id_2]  )
#                                        tia2.set_Id_3(   Id_3[ind_Id_3]  )
#                                        tia2.set_R_set(  R_set[ind_R_set])
#                                        tia2.upd()
#                
#                                        results[count][0] = tia2.get_FOM()     # FOM
#                                        results[count][1] = tia2.get_gain()    # gain
#                                        results[count][2] = tia2.get_f3dB()    # f3dB
#                                        results[count][3] = tia2.get_power()   # power
#                                        results[count][4]  = Vov_1[ind_Vov_1]  # Vov_1
#                                        results[count][5]  = Vov_2[ind_Vov_2]  # Vov_2
#                                        results[count][6]  = Vov_3[ind_Vov_3]  # Vov_3
#                                        results[count][7]  = Vov_L[ind_Vov_L]  # Vov_L
#                                        results[count][8]  = Vov_B[ind_Vov_B]  # Vov_B
#                                        results[count][9]  = Id_1[ind_Id_1]    # Id_1
#                                        results[count][10] = Id_2[ind_Id_2]    # Id_2
#                                        results[count][11] = Id_3[ind_Id_3]    # Id_3
#                                        results[count][12] = R_set[ind_R_set]  # R_set
#                
#                
#    print(f'results: {results}')
#    print(f'results: {results[10,:]}')
#    
#    """ sort results matrix
#    """
#    FOMs = results[:,0]
#    print(f'FOMs: {FOMs}')
#    argsort_FOMs = np.argsort(FOMs)
#    print(f'argsort_FOMs: {argsort_FOMs}')
#    
#    n_show = 3
#    top_results = np.zeros((n_show, 12))
#    for i in range(n_show):
#        print(f'{i} top result:')
#        top_results[i, :] = results[argsort_FOMs[-(i+1)], :]
#        print('\tFOM: %1.2f,  gain: %2.2f kohms,  f3dB: %3.2f MHz,  power: %2.2f mW' %(top_results[i, 0], top_results[i, 1]*1E-3, top_results[i, 2]*1E-6, top_results[i, 3]*1E+3 ))
##    print(f'top_results: {top_results}')
        
    
    
    
sweep_params()
   
        
#def calc_CG(Vov, WL, RL, Cin, Cout):
#    # gain and power stuff
#    gm = unCox * WL * Vov
#    Av = -gm * RL
#    Id = 0.5 * unCox * Vov**2
#    power = 5 * Id
#    print(f'gm: {gm}')
#    print(f'Av: {Av}')
#    print('Id: %3.2fuA' %(Id*1000000) )
#    print('Power: %3.2fmW' %(power*1000) )
#    
#    # BW stuff
#    Cg = f_cgtot(WL)*1E-15
#    Rin = 1/gm
#    b1 = Rin * Cg
#    w3dB = 1/b1
#    print(f'w3dB: {w3dB}')
#    
#calc_CG(0.2, 20, 10000, 100E-15, 100E-15)
    

#%% CG single, second way
# units of mv, kohms, uA

#gm_on_Id = np.linspace(10, 30, 5)
#print(f'gm_on_Id: {gm_on_Id}')
#
#Id = np.linspace(0.5, 2, 5)
#print(f'Id: {Id}')
#
#RL = np.linspace(5, 20, 5)
#print(f'RL: {RL}')


#class CG(tech_params):
#    
#    def __init__(self, Vov, WL, RL):
#        self.Vov = Vov
#        self.WL  = WL
#        self.RL  = RL
#        
#        self.gm = self.unCox * self.WL * self.Vov
#        self.Id = 0.5 * self.unCox * self.Vov**2
#        self.power = 5 * self.Id
#        
#        # BW stuff
#        self.Cg = self.f_cgtot(self.WL)*1E-15
#        self.Rin = 1/self.gm
#        self.b1 = self.Rin * self.Cg
#        self.w3dB = 1/self.b1
#        
#    def printout(self):
#        print(f'Vov: {self.Vov}v')
#        print(f'WL: {self.WL}')
#        print(f'RL: {(self.RL)/1000}kohms')
#        
#        print(f'gm: {self.gm}')
#        print('Id: %3.2fuA' %(self.Id*1000000) )
#        print('Power: %3.2fmW' %(self.power*1000) )
#        print(f'w3dB: {self.w3dB}')