import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

#%%


#%% constants
unCox = 0.00005



#%% CG single, first way

Vov1 = np.linspace(0.2, 0.3, 5)
print(f'Vov1: {Vov1}')
WL_1 = np.linspace(10, 100, 5)
print(f'WL_1: {WL_1}')
RL = np.linspace(5000, 20000, 5)
print(f'RL: {RL}')
print()

# input:  Vov, WL, RL
# output: gain, BW, power
class tech_params:
    unCox = 50E-6
    
    w     = [2,     5,     10,    20,    50,    100,  200,  500,  1000, 2000]
    cdtot = [5.60,  9.50,  16.0,  29.0,  68.0,  134,  263,  653,  1303, 2603]
    cgtot = [3.45,  8.60,  17.2,  34.5,  86.3,  173,  345,  863,  1727, 3453]
    cstot = [5.62,  9.50,  16.1,  29.2,  68.5,  133,  265,  658,  1313, 2623]
    cbtot = [10.6,  17.6,  29.2,  52.3,  122,   238,  469,  1164, 2322, 4640]
    cgs   = [1.02,  2.55,  5.10,  10.2,  25.5,  50,   102,  255,  510,  1020]
    cgd   = [1.00,  2.50,  5.00,  10.0,  25.0,  51,   100,  250,  500,  1000]
    f_cdtot = interpolate.interp1d(w, cdtot)
    f_cgtot = interpolate.interp1d(w, cgtot)
    f_cstot = interpolate.interp1d(w, cstot)
    f_cbtot = interpolate.interp1d(w, cbtot)
    f_cgs   = interpolate.interp1d(w, cgs  )
    f_cgd   = interpolate.interp1d(w, cgd  )
    
class amp(tech_params):
    def __init__(self, )
    
class CG(tech_params):
    
    def __init__(self, Vov, WL, RL):
        self.Vov = Vov
        self.WL  = WL
        self.RL  = RL
        
        self.gm = self.unCox * self.WL * self.Vov
        self.Id = 0.5 * self.unCox * self.Vov**2
        self.power = 5 * self.Id
        
        # BW stuff
        self.Cg = self.f_cgtot(self.WL)*1E-15
        self.Rin = 1/self.gm
        self.b1 = self.Rin * self.Cg
        self.w3dB = 1/self.b1
        
    def printout(self):
        print(f'Vov: {self.Vov}v')
        print(f'WL: {self.WL}')
        print(f'RL: {(self.RL)/1000}kohms')
        
        print(f'gm: {self.gm}')
        print('Id: %3.2fuA' %(self.Id*1000000) )
        print('Power: %3.2fmW' %(self.power*1000) )
        print(f'w3dB: {self.w3dB}')
        
CG_stage = CG(0.2, 20, 1000)
CG_stage.printout()
        
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









""" 
calculate the bandwidth

independent variables:
    1. Vov1
    2. Vov2
    3. Vov3
    4. W/L_1
    5. W/L_2
    6. W/L_3
    7. R1
    8. R2
    
independent variables:
    1. Id1
    2. Id2
    3. Id3
    4. gm1/Id1
    5. gm2/Id2
    6. gm3/Id3
    7. R1
    8. R2

"""
# give this function any six independent variables and it will calculate 
# the remaining 2 from the gain
#def calc_BW()