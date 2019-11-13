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
    def __init__(self):
        self.Vov1 = 0.2
        self.Vov2 = 0.2
        self.Vov3 = 0.2
        self.Id_1 = 0.0001
        self.Id_2 = 0.0001
        self.Id_3 = 0.0001
        self.R2   = 10000
        self.R3   = 10000
        
        self.gm1 = 0
        self.gm2 = 0
        self.gm3 = 0
        self.R1 = 0
        self.R4 = 0
        self.A  = 0
        
        self.WL_1 = 1
        self.WL_2 = 1
        self.WL_3 = 1
        
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
        
    def print_all(self):
        print(f'Vov1: {self.Vov1}')
        print(f'Vov2: {self.Vov2}')
        print(f'Vov3: {self.Vov3}')
        print(f'Id_1: {self.Id_1}')
        print(f'Id_2: {self.Id_2}')
        print(f'Id_3: {self.Id_3}')
        print(f'R2: {self.R2}')
        print(f'R3: {self.R3}')
        
        print(f'gm1: {self.gm1}')
        print(f'gm2: {self.gm2}')
        print(f'gm3: {self.gm3}')
        
        print(f'WL_1: {self.WL_1}')
        print(f'WL_2: {self.WL_2}')
        print(f'WL_3: {self.WL_3}')
        
        print(f'A: {self.A}')
        print('C1: %4.1ffF' %(1E+15*self.C1) )
        print('C2: %4.1ffF' %(1E+15*self.C2) )
        print('C3: %4.1ffF' %(1E+15*self.C3) )
        print('C4: %4.1ffF' %(1E+15*self.C4) )
        
        print(f'b11: {self.b11}')
        print(f'b12: {self.b12}')
        print(f'b13: {self.b13}')
        print(f'b14: {self.b14}')
        print(f'b1:  {self.b1}')
        print('w3dB: %4.1f MHz' %(self.w3dB/1000000) )
        print('tran_res_gain: %5.0f kohms' %( self.gm2 * self.R2 * self.R3 ) )
        print('power consumption: %2.2f mW' %( 5000*(self.Id_1 + self.Id_2 + self.Id_3)))
        
    def upd(self):
        self.tran_res = -0.84 * self.R1 * self.gm2 * self.R2
        
        # calc gm's and R's
        self.gm1 = 2*self.Id_1 / self.Vov1
        self.gm2 = 2*self.Id_2 / self.Vov2
        self.gm3 = 2*self.Id_3 / self.Vov3
        self.R1 = 1/self.gm1
        self.R4 = 1/self.gm3
        self.A = self.gm2 * self.R2
        
        # calc W/L's
        self.WL_1 = (2*self.Id_1) / (self.unCox*self.Vov1**2)
        self.WL_2 = (2*self.Id_2) / (self.unCox*self.Vov2**2)
        self.WL_3 = (2*self.Id_3) / (self.unCox*self.Vov3**2)
        
        # calc C's
        self.C1 = 1E-15*( 100 + self.f_cstot(self.WL_1) )
        self.C2 = 1E-15*( self.f_cdtot(self.WL_1) + self.f_cgs(self.WL_2) + (1+self.A)*self.f_cgd(self.WL_2) )
        self.C3 = 1E-15*( self.f_cdtot(self.WL_2) - self.f_cgd(self.WL_2) + (1 + 1/self.A)*self.f_cgd(self.WL_2) + self.f_cgd(self.WL_3) + 0.14*self.f_cgs(self.WL_3) )
        self.C4 = 1E-15*( -0.2*self.f_cgs(self.WL_3) + self.f_cstot(self.WL_3) - self.f_cgs(self.WL_3) + 250 )
        
        # calc b's
        self.b11 = self.R1 * self.C1
        self.b12 = self.R2 * self.C2
        self.b13 = self.R3 * self.C3
        self.b14 = self.R4 * self.C4
        self.b1  = self.b11 + self.b12 + self.b13 + self.b14
        self.w3dB = 1/self.b1
        
     
        
    def set_Vov1(self, Vov1):
        self.Vov1 = Vov1
    def set_Vov2(self, Vov2):
        self.Vov2 = Vov2
    def set_Vov3(self, Vov3):
        self.Vov3 = Vov3
        
    def set_Id_1(self, Id_1):
        self.Id_1 = Id_1
    def set_Id_2(self, Id_2):
        self.Id_2 = Id_2
    def set_Id_3(self, Id_3):
        self.Id_3 = Id_3
        
    def set_R1(self, R1):
        self.R1 = R1
    def set_R2(self, R2):
        self.R2 = R2
  
        
def unit_test_amp():
    amp1 = amp()
    amp1.upd()
    amp1.print_all()

unit_test_amp()
        
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