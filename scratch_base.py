import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

#%% power consumption
class tia_power_consumption:
    def __init__(self):
        self.Id_1 = 0
        self.Id_2 = 0
        self.Id_3 = 0
        self.R = 10000
        self.I_ref = 0
        self.power1 = 0
        self.power2 = 0
        self.power3 = 0
        self.power4 = 0
        
    def set_R(self, R):
        self.R = R
        
    def set_I_ref(self, I_ref):
        self.I_ref
        
    def set_Ids(self, Id_1, Id_2, Id_3):
        self.Id_1 = Id_1
        self.Id_2 = Id_2
        self.Id_3 = Id_3
        self.upd_power()
        
    def set_params(self, ratio_1_to_3, fraction_2 ):
        self.Id_1 = ( ratio_1_to_3 / (1 + ratio_1_to_3) )*(150E-6*(1-fraction_2))
        self.Id_2 = 150E-6*fraction_2
        self.Id_3 = ( 1 / (1 + ratio_1_to_3) )*(150E-6*(1-fraction_2))
        self.upd_power()
    
    def upd_power(self):
        self.power1 = 2*5*(self.Id_1 + self.Id_2 + self.Id_3) + 5*self.I_ref + 50/(4*self.R)
        self.power2 = 2*5*(self.Id_1 + self.Id_2 + self.Id_3) 
        self.power3 = 5*self.I_ref
        self.power4 = 50/(4*self.R)
        
    def get_power(self):
        return [self.power1, self.power2, self.power3, self.power4]
    
    def print_power(self, n):
        if n >= 0:
            print('total power consumption: %3.3f mw' %(self.power1/1000) )
        if n == 1:
            print('power consumption Ids:   %3.3f mw' %(self.power2/1000) )
            print('power consumption I_ref: %3.3f mw' %(self.power3/1000) )
            print('power consumption R:     %3.3f mw' %(self.power4/1000) )

def unit_test_tia_power_consumption(test_type):
    if test_type >= 0:
        tia_power_module_test = tia_power_consumption()
        tia_power_module_test.set_params(1.2, 5E-5)
        tia_power_module_test.print_power(1)
    

#unit_test_tia_power_consumption(1)

#%% MOSFET

class mosfet:
    
    unCox = 50E-6
    upCox = 25E-6
    lambd = 0.1
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
    
    # lookup table for l=2um
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
    
    def __init__(self, mosfet_type, debug):  # 0=nmos l=1, 1=nmos l=2, 2=pmos l=1, 3=pmos l=2
        # set variables
        self.Vov = 0.2
        self.Id  = 1E-5
        self.type= mosfet_type
        self.debug = debug
        self.type_dict = {0:'nmos l=1u', 1:'nmos l=2u', 2:'pmos l=1u', 3:'pmos l=2u',}
        # derived variables
        self.WL  = -1   
        self.L   = -1 # um
        self.W   = -1 # um
        self.gm  = -1 # A/V
        self.gmp = -1 # gm + gmb = 1.2*gm
        self.ro  = -1
        self.cgs = -1
        self.cgd = -1
        self.cgb = -1
        self.csb = -1
        self.cdb = -1
        
    def print_basic(self):
        print(f'type: {self.type_dict[self.type]}')
        print('input variables: \n  Vov= %2.2fV, \n  Id =%3.0fuA' %(self.Vov, 1E+6*self.Id) )
        print('calc parameters: \n  WL = %2.2f,  \n  W  = %1.1fum,  \n  L  = %1.0fum,  \n  gm = %1.2fmA/V' 
              %(self.WL, self.W, self.L, 1000*self.gm) )
    
    def print_caps(self):
        print('cgs: %3.1f fF' %self.cgs)
        print('cgd: %3.1f fF' %self.cgd)
        print('cgb: %3.1f fF' %self.cgb)
        print('csb: %3.1f fF' %self.csb)
        print('cdb: %3.1f fF\n' %self.cdb)
        
    def print_all(self):
        self.print_basic()
        self.print_caps()
    
    def set_WL(self, WL):
        self.WL = WL
        self.upd_caps()
        
    def upd_caps(self):
        if self.type == 0 or self.type == 2: 
            self.cgs = self.cgs_lookup_1u(self.WL)
            self.cgd = self.cgd_lookup_1u(self.WL)
            self.cgb = self.cgtot_lookup_1u(self.WL) - self.cgs - self.cgd
            self.csb = self.cstot_lookup_1u(self.WL) - self.cgs
            self.cdb = self.cdtot_lookup_1u(self.WL) - self.cgd
        if self.type == 1 or self.type == 3: 
            self.cgs = self.cgs_lookup_2u(self.WL)
            self.cgd = self.cgd_lookup_2u(self.WL)
            self.cgb = self.cgtot_lookup_2u(self.WL) - self.cgs - self.cgd
            self.csb = self.cstot_lookup_2u(self.WL) - self.cgs
            self.cdb = self.cdtot_lookup_2u(self.WL) - self.cgd
        
    def set_params(self, Vov, Id):
        self.Vov = Vov
        self.Id  = Id
        # WL = 2*Id / unCox*Vov^2
        if self.type == 0 or self.type == 2:
            self.WL = 2*self.Id / ( self.unCox*self.Vov**2 )
            self.W  = self.WL*(self.type+1)
        if self.type == 1 or self.type == 3:
            self.WL = 2*self.Id / ( self.upCox*self.Vov**2 )
            self.W  = self.WL*(self.type-1)
        # gm = 2*Id/Vov
        self.gm = 2*self.Id/self.Vov
        self.gmp= 1.2*self.gm
        self.ro = 2/(self.lambd*self.Id)
        self.upd_caps()
        if self.debug:
            self.print_all()
        
    def get_gm(self):
        return self.gm
        
    def get_params(self):
        return [self.Vov, self.Id, self.WL, self.W, self.gm, self.ro]
        
    def get_caps(self):
        return [self.cgs, self.cgd, self.cgb, self.csb, self.cdb]
    
    def sweep_caps(self):
        WL_sweep_1u  = [2,     5,     10,    20,    50,    100,  200,  500,  1000, 2000]
        WL_sweep_2u  = [1,     2.5,   5,     10,    25,    50,    100,   250]
        for self.type in range(4):
            if self.type == 0 or self.type == 2:
                sweep = np.zeros((10, 6))
                sweep[:,0] = WL_sweep_1u
            if self.type == 1 or self.type ==3:
                sweep = np.zeros((8,6))
                sweep[:,0] = WL_sweep_2u
                
            if self.debug: print(f'sweep: \n{sweep}')
            for i in range(sweep[:,0].shape[0]):  # iterate rows
                self.set_WL(sweep[i,0])
                sweep[i,0] = sweep[i,0]
                sweep[i, 1:6]   = self.get_caps()
                
            
            
            fig = plt.figure(self.type+1)
            fig.clf()
            fig.suptitle(self.type_dict[self.type])
            ax = fig.subplots(1,1)
            ax.set_xlim([0,100])
            ax.set_ylim([0,100])
            ax.set_xlabel('W/L')
            ax.set_ylabel('C [fF]')
            
            names_list = ['cgs', 'cgd', 'cgb', 'csb', 'cdb']
            trend_list = []
            trend_max = 100
            for i in range(5):
                trend_list.append( np.polyfit(sweep[:,0], sweep[:,i+1], 1) )
                ax.plot(sweep[:,0], sweep[:,i+1], label='%s=%3.2f %3.2f' %(names_list[i], trend_list[i][0], trend_list[i][1]))
                ax.plot([0, trend_max], [trend_list[i][1], trend_max*trend_list[i][0]+trend_list[i][1]], 'k--', linewidth=0.5)
            ax.legend()
            ax.grid()
    
        
def unit_test_mosfet(test_type):
    nmos_1 = mosfet(0, 0)
    nmos_1.print_all()
    
#unit_test_mosfet(1)

#%%
class ee:
    @staticmethod
    def parallel(R1, R2):
        return (R1*R2)/(R1+R2)
    
    @staticmethod
    def print_R(name, val):
        if val > 1E+19:
            print('%s = inf' %name)
        if val <= 1E+19 and val > 1000:
            print('%s = %2.3f kohms' %(name, val/1000))
        if val <= 1000:
            print('%s = %3.0f ohms' %(name, val))
            
    @staticmethod
    def print_C(name, val):
        if val > 1E+6:
            print('%s = %3.0f uF' %(name, val/1E+6))
        if val <= 1E+6:
            print('%s = %3.0f fF' %(name, val))
            
    @staticmethod
    def print_A(name, units, val):
        if val > 1E+3:
            print('%s = %3.2f k%s' %(name, val/1E+3, units))
        if val <= 1E+3:
            print('%s = %3.2f %s' %(name, val, units))

#%% CG
class CG(mosfet):
    
    def __init__(self):
        self.M1   = mosfet(0, 0)
        self.RL   = -1
        self.TI   = -1
        self.Rin  = -1
        self.Cin  = -1
        self.Rout = -1
        self.Cout = -1
        
    def _set(self, Vov, Id, RL):
        self.RL    = RL
        self.M1.set_params(Vov, Id)
        
        self.Rin  = 1/self.M1.gmp
        self.Rout = ee.parallel(2*self.M1.ro, self.RL)
        self.Cin  = self.M1.cgs + self.M1.csb
        self.Cout = self.M1.cgd + self.M1.cdb
        self.TI   = self.Rout
        self._print()
        
        return [self.TI, self.Rin, self.Cin, self.Rout, self.Cout]
        
    def _print(self):
        print('CG stage parameters:')
        ee.print_A('   A', 'V/A', self.TI)
        ee.print_R('   Rin ', self.Rin)
        ee.print_C('   Cin ', self.Cin)
        ee.print_R('   Rout', self.Rout)
        ee.print_C('   Cout', self.Cout)
    
    
cg_stage = CG()
ret = cg_stage._set(0.2, 1E-4, 10000)
print()
        

#%% CS
class CS(mosfet):
    
    def __init__(self):
        self.M2   = mosfet(0, 0)
        self.A1   = -1
        self.Rin  = -1
        self.Cin  = -1
        self.Rout = -1
        self.Cout = -1
#        self.M.print_all()
        
    def _set(self, Vov, Id, gain):  # A = -gain
        self.A1 = -gain
        self.M2.set_params(Vov, Id)
        
        self.Rin  = 1E+21
        self.Rout = ee.parallel(2*self.M2.ro, 1/(self.A1*self.M2.gmp))
        self.Cin  = self.M2.cgs + self.M2.cgb + (1+self.A1)*self.M2.cgd
        self.Cout = (1+1/self.A1)*self.M2.cgd + self.M2.cdb
        self._print()
        
        return [self.A1, self.Rin, self.Cin, self.Rout, self.Cout]
        
    def _print(self):
        print('CS stage parameters:')
        ee.print_A('   A', 'V/V', self.A1)
        ee.print_R('   Rin ', self.Rin)
        ee.print_C('   Cin ', self.Cin)
        ee.print_R('   Rout', self.Rout)
        ee.print_C('   Cout', self.Cout)
    
    
cs_stage = CS()
ret = cs_stage._set(0.2, 1E-5, -4)
print()
    

#%% CD
class CD(mosfet):
    
    def __init__(self):
        self.M3    = mosfet(0, 0)
        self.A2    = -0.84
        self.Rin  = -1
        self.Cin  = -1
        self.Rout = -1
        self.Cout = -1
#        self.M.print_all()
        
    def _set(self, Vov, Id): 
        self.M3.set_params(Vov, Id)
        
        self.Rin  = 1E+21
        self.Rout = ee.parallel(2*self.M3.ro, 0.84/2*self.M3.Vov/self.M3.Id)
        self.Cin  = self.M3.cgs + self.M3.cgb + (1+self.A2)*self.M3.cgd
        self.Cout = (1+1/self.A2)*self.M3.cgd + self.M3.cdb
        self._print()
        
        return [self.A2, self.Rin, self.Cin, self.Rout, self.Cout]
        
    def _print(self):
        print('CS stage parameters:')
        ee.print_A('   A', 'V/V', self.A2)
        ee.print_R('   Rin ', self.Rin)
        ee.print_C('   Cin ', self.Cin)
        ee.print_R('   Rout', self.Rout)
        ee.print_C('   Cout', self.Cout)
    
    
cd_stage = CD()
ret = cd_stage._set(0.2, 1E-5)
print()














