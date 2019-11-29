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
        # check WL limits
        if self.WL <= 2 or self.WL > 2000:
            print(f'out of range: \n   Vov = {self.Vov} \n   Id = {self.Id} \n   WL={self.WL}')
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
        
    def _print(self):
        print('CG stage parameters:')
        ee.print_A('   TI', 'ohms', self.TI)
        ee.print_R('   Rin ', self.Rin)
        ee.print_R('   Rout', self.Rout)
        ee.print_C('   Cin ', self.Cin)
        ee.print_C('   Cout', self.Cout)
        
    def _set(self, Vov, Id, RL):
        self.RL    = RL
        self.M1.set_params(Vov, Id)
        self.Rin  = 1/self.M1.gmp
        self.Rout = ee.parallel(2*self.M1.ro, self.RL)
        self.Cin  = self.M1.cgs + self.M1.csb
        self.Cout = self.M1.cgd + self.M1.cdb
        self.TI   = self.Rout
        
    def get_Rin(self):
        return self.Rin
    def get_Rout(self):
        return self.Rout
    def get_Cin(self):
        return self.Cin
    def get_Cout(self):
        return self.Cout
    def get_TI(self):
        return self.TI
        
    
def unit_test_CG():
    print('--- CG unit test --------------')
    cg = CG()
    cg._set(0.2, 1E-4, 10000)
    cg._print()


unit_test_CG()

#%% CG_lookup.  abstract CG stage to lookup3(Vov, Id, RL)
class CG_LK(CG):
    _Vov = np.linspace(0.15, 0.3,  5)
    _Id  = np.linspace(1E-5, 1E-4, 5)
    _RL  = np.linspace(1E+4, 4E+4, 5)
#    print(f'_Vov: {_Vov}')
#    print(f'_Id:  {_Id}')
#    print(f'_RL:  {_RL}')
    print()

    def __init__(self):
        super().__init__()
        self.n_Vov = self._Vov.shape[0]
        self.n_Id  = self._Id.shape[0]
        self.n_RL= self._RL.shape[0]
#        print(f'n_Vov: {self.n_Vov}')
#        print(f'n_Id:  {self.n_Id}')
#        print(f'n_RL:  {self.n_RL}')
        
        self.Rin_data  = np.zeros((self.n_Vov, self.n_Id, self.n_RL))
        self.Rout_data = np.zeros((self.n_Vov, self.n_Id, self.n_RL))
        self.Cin_data  = np.zeros((self.n_Vov, self.n_Id, self.n_RL))
        self.Cout_data = np.zeros((self.n_Vov, self.n_Id, self.n_RL))
        self.TI_data   = np.zeros((self.n_Vov, self.n_Id, self.n_RL))
        
        for i in range(self.n_Vov):
            for j in range(self.n_Id):
                for k in range(self.n_RL):
#                    print(f'i: {i}, j: {j}, k: {k}')
                    super()._set(self._Vov[i], self._Id[j], self._RL[k])
                    self.Rin_data [i,j,k] = super().get_Rin()
                    self.Rout_data[i,j,k] = super().get_Rout()
                    self.Cin_data [i,j,k] = super().get_Cin()
                    self.Cout_data[i,j,k] = super().get_Cout()
                    self.TI_data  [i,j,k] = super().get_TI()
        self.Rin_f  = RegularGridInterpolator((self._Vov, self._Id, self._RL), self.Rin_data)
        self.Rout_f = RegularGridInterpolator((self._Vov, self._Id, self._RL), self.Rout_data)
        self.Cin_f  = RegularGridInterpolator((self._Vov, self._Id, self._RL), self.Cin_data)
        self.Cout_f = RegularGridInterpolator((self._Vov, self._Id, self._RL), self.Cout_data)
        self.TI_f   = RegularGridInterpolator((self._Vov, self._Id, self._RL), self.TI_data)
                
    def get_Rin(self, Vov, Id, RL):
        return self.Rin_f([Vov, Id, RL])
    def get_Rout(self, Vov, Id, RL):
        return self.Rout_f([Vov, Id, RL])
    def get_Cin(self, Vov, Id, RL):
        return self.Cin_f([Vov, Id, RL])
    def get_Cout(self, Vov, Id, RL):
        return self.Cout_f([Vov, Id, RL])
    def get_TI(self, Vov, Id, RL):
        return self.TI_f([Vov, Id, RL])

def CG_LK_unit_test(Vov, Id, RL):
    print('--- CG_LK_unit_test ----------------')
    cg = CG_LK()
    ee.print_A('   TI', 'ohms', cg.get_TI  (Vov, Id, RL))
    ee.print_R('   Rin ', cg.get_Rin (Vov, Id, RL))
    ee.print_R('   Rout', cg.get_Rout(Vov, Id, RL))
    ee.print_C('   Cin ', cg.get_Cin (Vov, Id, RL))
    ee.print_C('   Cout', cg.get_Cout(Vov, Id, RL))
    print('------------------------------------\n')
    
CG_LK_unit_test(0.2, 1E-4, 10000)


#%% CS
class CS(mosfet):
    
    def __init__(self):
        self.M2   = mosfet(0, 0)
        self.A1   = -1
        self.Rin  = -1
        self.Cin  = -1
        self.Rout = -1
        self.Cout = -1

    def _print(self):
        print('CS stage parameters:')
        ee.print_A('   |A| ', 'V/V', self.A1)
        ee.print_R('   Rin ', self.Rin)
        ee.print_R('   Rout', self.Rout)
        ee.print_C('   Cin ', self.Cin)
        ee.print_C('   Cout', self.Cout)
        
    def _set(self, Vov, Id, A1):  # A = -gain
        self.A1 = A1
        self.M2.set_params(Vov, Id)
        
        self.Rin  = 1E+21
        self.Rout = ee.parallel(2*self.M2.ro, 1/(self.A1*self.M2.gmp))
        self.Cin  = self.M2.cgs + self.M2.cgb + (1+self.A1)*self.M2.cgd
        self.Cout = (1+1/self.A1)*self.M2.cgd + self.M2.cdb
#        self._print()
        
    def get_Rin(self):
        return self.Rin
    def get_Rout(self):
        return self.Rout
    def get_Cin(self):
        return self.Cin
    def get_Cout(self):
        return self.Cout
    def get_A1(self):
        return self.A1

def unit_test_CS():
    print('--- CS_unit_test --------------------')
    cs = CS()
    cs._set(0.2, 1E-4, 3)
    cs._print()
    print()

unit_test_CS()


#%%
class CS_LK(CS):
    _Vov = np.linspace(0.15, 0.3,  5)
    _Id  = np.linspace(1E-5, 1E-4, 5)
    _A1  = np.linspace(1, 5, 5)
#    print(f'_Vov: {_Vov}')
#    print(f'_Id:  {_Id}')
#    print(f'_gain:{_gain}')
    
    def __init__(self):
        super().__init__()
        self.n_Vov = self._Vov.shape[0]
        self.n_Id  = self._Id.shape[0]
        self.n_A1= self._A1.shape[0]
        print(f'n_Vov: {self.n_Vov}')
        print(f'n_Id:  {self.n_Id}')
        print(f'n_A1:  {self.n_A1}')
        
        self.Rin_data  = np.zeros((self.n_Vov, self.n_Id, self.n_A1))
        self.Rout_data = np.zeros((self.n_Vov, self.n_Id, self.n_A1))
        self.Cin_data  = np.zeros((self.n_Vov, self.n_Id, self.n_A1))
        self.Cout_data = np.zeros((self.n_Vov, self.n_Id, self.n_A1))
        self.A1_data   = np.zeros((self.n_Vov, self.n_Id, self.n_A1))
        
        for i in range(self.n_Vov):
            for j in range(self.n_Id):
                for k in range(self.n_A1):
#                    print(f'i: {i}, j: {j}, k: {k}')
                    super()._set(self._Vov[i], self._Id[j], self._A1[k])
                    self.Rin_data [i,j,k] = super().get_Rin()
                    self.Rout_data[i,j,k] = super().get_Rout()
                    self.Cin_data [i,j,k] = super().get_Cin()
                    self.Cout_data[i,j,k] = super().get_Cout()
                    self.A1_data  [i,j,k] = super().get_A1()
        self.Rin_f  = RegularGridInterpolator((self._Vov, self._Id, self._A1), self.Rin_data)
        self.Rout_f = RegularGridInterpolator((self._Vov, self._Id, self._A1), self.Rout_data)
        self.Cin_f  = RegularGridInterpolator((self._Vov, self._Id, self._A1), self.Cin_data)
        self.Cout_f = RegularGridInterpolator((self._Vov, self._Id, self._A1), self.Cout_data)
        self.A1_f   = RegularGridInterpolator((self._Vov, self._Id, self._A1), self.A1_data)
                
    def get_Rin(self, Vov, Id, A1):
        return self.Rin_f ([Vov, Id, A1])
    def get_Rout(self, Vov, Id, A1):
        return self.Rout_f([Vov, Id, A1])
    def get_Cin(self, Vov, Id, A1):
        return self.Cin_f ([Vov, Id, A1])
    def get_Cout(self, Vov, Id, A1):
        return self.Cout_f([Vov, Id, A1])
    def get_A1(self, Vov, Id, A1):
        return self.A1_f([Vov, Id, A1])
    
def CS_LK_unit_test():
    print('--- CS_LK_unit_test ----------------')
    cs = CS_LK()
    ee.print_A('   |A| ', 'V/V', cs.get_A1(0.2, 1E-4, 3))
    ee.print_R('   Rin ', cs.get_Rin(0.2,  1E-4, 3))
    ee.print_R('   Rout', cs.get_Rout(0.2, 1E-4, 3))
    ee.print_C('   Cin ', cs.get_Cin(0.2,  1E-4, 3))
    ee.print_C('   Cout', cs.get_Cout(0.2, 1E-4, 3))
    print('------------------------------------\n')
    
CS_LK_unit_test()
    

#%% CD
class CD(mosfet):
    
    def __init__(self):
        self.M3    = mosfet(0, 0)
        self.A2    = -0.84
        self.Rin  = -1
        self.Cin  = -1
        self.Rout = -1
        self.Cout = -1
        
    def _print(self):
        print('CD stage parameters:')
        ee.print_A('    |A| ', 'V/V', -self.A2)
        ee.print_R('    Rin ', self.Rin)
        ee.print_R('    Rout', self.Rout)
        ee.print_C('    Cin ', self.Cin)
        ee.print_C('    Cout', self.Cout)
        
    def _set(self, Vov, Id): 
        self.M3.set_params(Vov, Id)
        self.Rin  = 1E+21
        self.Rout = ee.parallel(2*self.M3.ro, 1/self.M3.gmp)
        self.Cin  = self.M3.cgs + self.M3.cgb + (1+self.A2)*self.M3.cgd
        self.Cout = (1+1/self.A2)*self.M3.cgd + self.M3.cdb
#        self._print()
    
    def get_A2(self):
        return self.A2
    def get_Rin(self):
        return self.Rin
    def get_Rout(self):
        return self.Rout
    def get_Cin(self):
        return self.Cin
    def get_Cout(self):
        return self.Cout
  
def unit_test_CD():
    print('--- CD_unit_test ---------------------')
    cd = CD()
    cd._set(0.2, 1E-4)
    cd._print()
    print('')

unit_test_CD()


#%%
class CD_LK(CD):
    _Vov = np.linspace(0.15, 0.3,  5)
    _Id  = np.linspace(1E-5, 1E-4, 5)
#    print(f'_Vov: {_Vov}')
#    print(f'_Id:  {_Id}')
    
    def __init__(self):
        super().__init__()
        self.n_Vov = self._Vov.shape[0]
        self.n_Id  = self._Id.shape[0]
        self.A2_data   = np.zeros((self.n_Vov, self.n_Id))
        self.Rin_data  = np.zeros((self.n_Vov, self.n_Id))
        self.Rout_data = np.zeros((self.n_Vov, self.n_Id))
        self.Cin_data  = np.zeros((self.n_Vov, self.n_Id))
        self.Cout_data = np.zeros((self.n_Vov, self.n_Id))
        print(f'n_Vov: {self.n_Vov}')
        print(f'n_Id:  {self.n_Id}')
        for i in range(self.n_Vov):
            for j in range(self.n_Id):
#                print(f'i: {i}, j: {j}')
                super()._set(self._Vov[i], self._Id[j])
                self.A2_data  [i,j] = super().get_A2()
                self.Rin_data[i,j]  = super().get_Rin()
                self.Rout_data[i,j] = super().get_Rout()
                self.Cin_data[i,j]  = super().get_Cin()
                self.Cout_data[i,j] = super().get_Cout()
        self.A2_f   = RegularGridInterpolator((self._Vov, self._Id), self.A2_data)
        self.Rin_f  = RegularGridInterpolator((self._Vov, self._Id), self.Rin_data)
        self.Rout_f = RegularGridInterpolator((self._Vov, self._Id), self.Rout_data)
        self.Cin_f  = RegularGridInterpolator((self._Vov, self._Id), self.Cin_data)
        self.Cout_f = RegularGridInterpolator((self._Vov, self._Id), self.Cout_data)
    
    def get_A2(self, Vov, Id):
        return self.A2_f([Vov, Id])
    def get_Rin(self, Vov, Id):
        return self.Rin_f([Vov, Id])
    def get_Rout(self, Vov, Id):
        return self.Rout_f([Vov, Id])
    def get_Cin(self, Vov, Id):
        return self.Cin_f([Vov, Id])
    def get_Cout(self, Vov, Id):
        return self.Cout_f([Vov, Id])
        
def CD_LK_unit_test():
    print('--- CD_LK_unit_test ----------------')
    cd = CD_LK()
    ee.print_A('   |A| ', 'V/V', -cd.get_A2(0.2, 1E-4))
    ee.print_R('   Rin ', cd.get_Rin(0.2, 1E-4))
    ee.print_R('   Rout', cd.get_Rout(0.2, 1E-4))
    ee.print_C('   Cin ', cd.get_Cin(0.2, 1E-4))
    ee.print_C('   Cout', cd.get_Cout(0.2, 1E-4))
    print('------------------------------------')
    
CD_LK_unit_test()

##%% TIA
#class TIA(CG, CS, CD):
#    def __init__(self):
#        self.CG = CG()
#        self.CS = CS()
#        self.CD = CD()
#        self.Cin  = 100
#        self.Rout = 1E+4
#        self.Cout = 500
        












