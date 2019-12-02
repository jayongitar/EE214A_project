import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import time


#%%
class ee:
    @staticmethod
    def parallel(R1, R2):
        if R1 > 1E+20 and R2 > 1E+20:
            return 1E+21
        if R1 > 1E+20 and R2 < 1E+20:
            return R2
        if R2 > 1E+20 and R1 < 1E+20:
            return R1
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
    def print_F(name, val):
        print('%s = %3.2f MHz' %(name, val*1E-6))
            
    @staticmethod
    def print_V(name, val):
        print('%s = %1.3f V' %(name, val))
            
    @staticmethod
    def print_I(name, val):
        if val < 1E-3:
            print('%s = %3.1f uA' %(name, val*1E+6))
        if val >= 1E-3:
            print('%s = %3.2f mA' %(name, val*1E+3))
            
    @staticmethod
    def print_A(name, units, val):
        if val > 1E+3:
            print('%s = %3.2f k%s' %(name, val/1E+3, units))
        if val <= 1E+3:
            print('%s = %3.2f %s' %(name, val, units))

#%% MOSFET

class mosfet:
    
    unCox = 50E-6
    upCox = 25E-6
    lambd = 0.1
    
    coef_n1 = np.array([
            [2.033, 0.000],
            [0.505, 0.000],
            [0.800, 3.000],
            [0.553, 2.188]
            ])
    
    coef_n2 = np.array([
            [3.567, 0.000],
            [0.509, 0.000],
            [0.800, 3.000],
            [0.553, 2.195]
            ])
    
    coef_n4 = np.array([
            [6.633, 0.000],
            [0.506, 0.000],
            [0.8,   3.000],
            [0.678, 2.609]
            ])
    
    coef_p1 = np.array([
            [2.033, 0.000],
            [0.502, 0.000],
            [1.250, 2.100],
            [1.033, 1.826]
            ])
    
    coef_p2 = np.array([
            [3.567, 0.000],
            [0.509, 0.000],
            [1.250, 2.100],
            [0.816, 1.541]
            ])
    
    coef_p4 = np.array([
            [6.633, 0.000],
            [0.506, 0.000],
            [1.250, 2.100],
            [1.033, 1.826]
            ])
    
    
    def __init__(self, name, mosfet_type, debug):  # 0=nmos l=1, 1=nmos l=2, 2=pmos l=1, 3=pmos l=2
        # set variables
        self.name = name
        self.Vov = 0.2
        self.Id  = 2E-5
        self.type= mosfet_type
        self.debug = debug
        self.type_dict = {1:'nmos l=1u', 2:'nmos l=2u', 4:'nmos l=4u', 5:'pmos l=1u', 6:'pmos l=2u', 8:'pmos l=4u'}
        # derived variables
        self.WL  = -1   
        self.W   = -1 # um
        self.L   = -1 # um
        self.gm  = -1 # A/V
        self.gmp = -1 # gm + gmb = 1.2*gm
        self.ro  = -1
        self.cgs = -1
        self.cgd = -1
        self.csb = -1
        self.cdb = -1
        
    def _print(self):
        print(f'mosfet {self.name}')
        print(f' type: {self.type_dict[self.type]}')
        print(' input variables:')
        print('  Vov= %2.2fV' %self.Vov)
        print('  Id =%3.0fuA' %(1E+6*self.Id))
        print(' calc parameters:')
        print('  W  = %1.1fum' %self.W)
        print('  L  = %1.0fum' %self.L)
        print('  WL = %2.1f'   %self.WL)
        print('  gm = %1.2fmA/V, %2.2f kohms' %(1E+3*self.gm, 1E-3*(1/self.gm)))
        ee.print_R('  ro', self.ro)
    
    def _print_caps(self):
        print('caps:')
        print('  cgs: %3.1f fF' %self.cgs)
        print('  cgd: %3.1f fF' %self.cgd)
        print('  csb: %3.1f fF' %self.csb)
        print('  cdb: %3.1f fF\n' %self.cdb)
        
    def _print_all(self):
        print('-'*30)
        self._print()
        self._print_caps()
        
    def _set(self, Vov, Id):
        self.Vov = Vov
        self.Id  = Id
        # WL = 2*Id / unCox*Vov^2
        # 0=nmos l=1, 1=nmos l=2, 2=pmos l=1, 3=pmos l=2
        if self.type == 1 or self.type == 2 or self.type == 4:
            self.WL = 2*self.Id / ( self.unCox*self.Vov**2 )
            self.L  = self.type
            self.W  = self.WL*(self.L)
        if self.type == 5 or self.type ==6 or self.type == 8:
            self.WL = 2*self.Id / ( self.upCox*self.Vov**2 )
            self.L  = self.type-4
            self.W  = self.WL*(self.L)
            
        # check WL limits
        if self.W < 0.1 or self.W > 1000:
            print(f'mosfet {self.name}: out of range parameter: \n   Vov = {self.Vov} \n   Id = {self.Id} \n   WL={self.WL}')
            return -1
            
        # gm = 2*Id/Vov
        self.gm = 2*self.Id/self.Vov
        self.gmp= 1.2*self.gm
        self.ro = 1/(self.lambd*self.Id)
        self.upd_caps()
        if self.debug:
            self._print_all()
        return 0
             
    def upd_caps(self):
#        print(self.coef_n1.shape)
        try:
            if self.type == 1: 
                self.cgs = self.coef_n1[0, 0]*self.W + self.coef_n1[0, 1]
                self.cgd = self.coef_n1[1, 0]*self.W + self.coef_n1[1, 1]
                self.csb = self.coef_n1[2, 0]*self.W + self.coef_n1[2, 1]
                self.cdb = self.coef_n1[3, 0]*self.W + self.coef_n1[3, 1]
            if self.type == 2:
                self.cgs = self.coef_n2[0, 0]*self.W + self.coef_n2[0, 1]
                self.cgd = self.coef_n2[1, 0]*self.W + self.coef_n2[1, 1]
                self.csb = self.coef_n2[2, 0]*self.W + self.coef_n2[2, 1]
                self.cdb = self.coef_n2[3, 0]*self.W + self.coef_n2[3, 1]
            if self.type == 4:
                self.cgs = self.coef_n4[0, 0]*self.W + self.coef_n4[0, 1]
                self.cgd = self.coef_n4[1, 0]*self.W + self.coef_n4[1, 1]
                self.csb = self.coef_n4[2, 0]*self.W + self.coef_n4[2, 1]
                self.cdb = self.coef_n4[3, 0]*self.W + self.coef_n4[3, 1]
                
            if self.type == 5:
                self.cgs = self.coef_p1[0, 0]*self.W + self.coef_p1[0, 1]
                self.cgd = self.coef_p1[1, 0]*self.W + self.coef_p1[1, 1]
                self.csb = self.coef_p1[2, 0]*self.W + self.coef_p1[2, 1]
                self.cdb = self.coef_p1[3, 0]*self.W + self.coef_p1[3, 1]
            if self.type == 6: 
                self.cgs = self.coef_p2[0, 0]*self.W + self.coef_p2[0, 1]
                self.cgd = self.coef_p2[1, 0]*self.W + self.coef_p2[1, 1]
                self.csb = self.coef_p2[2, 0]*self.W + self.coef_p2[2, 1]
                self.cdb = self.coef_p2[3, 0]*self.W + self.coef_p2[3, 1]
            if self.type == 8: 
                self.cgs = self.coef_p4[0, 0]*self.W + self.coef_p4[0, 1]
                self.cgd = self.coef_p4[1, 0]*self.W + self.coef_p4[1, 1]
                self.csb = self.coef_p4[2, 0]*self.W + self.coef_p4[2, 1]
                self.cdb = self.coef_p4[3, 0]*self.W + self.coef_p4[3, 1]
        except:
            print('error in mosfet.upd_caps()')
            
    def _get_W(self):
        return self.W
    def _get_L(self):
        return self.L
    def _get(self):
        return [np.round(1E+0*self.Vov, 2), 
                np.round(1E+6*self.Id,  2), 
                np.round(1E+0*self.type,2), 
                np.round(1E+0*self.WL,  2), 
                np.round(1E+0*self.W,   2), 
                np.round(1E+0*self.L,   2), 
                np.round(1E+3*self.gm,  2), 
                np.round(1E+3*self.gmp, 2), 
                np.round(1E-3*self.ro,  2), 
                np.round(1E+0*self.cgs, 2), 
                np.round(1E+0*self.cgd, 2), 
                np.round(1E+0*self.csb, 2), 
                np.round(1E+0*self.cdb, 2)]
        
def unit_test_mosfet(test_type):
    nmos_1 = mosfet(0, 0)  # (type, debug)
    nmos_1._set(0.2, 1E-5)
    nmos_1.print_all()
    
#unit_test_mosfet(1)

#%% CG
class CG():
    def __init__(self):
        print_mosfet = 0
        self.M1L  = mosfet('M1L', 6, print_mosfet)
        self.M1   = mosfet('M1',  1, print_mosfet)
        self.M1B  = mosfet('M1B', 2, print_mosfet)
        self.Id_1 = -1
        self.R_LCG= -1
        self.TI   = -1
        self.Rin  = -1
        self.Cin  = -1
        self.Rout = -1
        self.Cout = -1
        self.C_IN = 100
        
    def _print(self):
        print('CG stage parameters:')
        ee.print_A('   TI', 'ohms', self.TI)
        ee.print_R('   Rin ', self.Rin)
        ee.print_R('   Rout', self.Rout)
        ee.print_C('   Cin ', self.Cin)
        ee.print_C('   Cout', self.Cout)
        
    def _set(self, Vov_1, Vov_N, Vov_P, Id_1, R_LCG):
        self.Id_1 = Id_1
        self.R_LCG = R_LCG
        r1 = self.M1._set (Vov_1, Id_1)
        r2 = self.M1L._set(Vov_P,  Id_1)
        r3 = self.M1B._set(Vov_N,  Id_1)
        if r1 == -1 or r2 == -1 or r3 == -1:
            print('failed to set a mosfet parameters in CG.')
            return -1  
#        self.Rin_approx = 1/self.M1.gmp
#        self.Cin_approx = self.M1.cgs + self.M1.csb
#        self.Cout_approx = self.M1.cgd + self.M1.cdb
        self.Rin  = (self.M1.ro + 2*self.R_LCG) / (self.M1.gmp*(self.M1.ro + self.R_LCG))
        self.Rout = ee.parallel(2*self.M1.ro, self.R_LCG)
        self.Cin  = self.M1.cgs + self.M1.csb + self.M1B.cgd + self.M1B.cdb + self.C_IN
        self.Cout = self.M1.cgd + self.M1.cdb + self.M1L.cgd + self.M1L.cdb
        self.TI   = self.Rout
        return 0
        
    def get_TI(self):
        return self.TI
    def get_Rin(self):
        return self.Rin
    def get_Rout(self):
        return self.Rout
    def get_Cin(self):
        return self.Cin
    def get_Cout(self):
        return self.Cout
        
def unit_test_CG():
    print('--- CG unit test --------------')
    cg = CG()
    r = cg._set(0.2, 0.4, 0.4, 1E-5, 10000)
#    print(f'r: {r}')
    cg._print()

#unit_test_CG()

#%% CS
class CS(mosfet):
    def __init__(self):
        print_mosfet = 0
        self.M2   = mosfet('M2',  1, print_mosfet)
        self.M2L  = mosfet('M2L', 4, print_mosfet)
        self.M2B  = mosfet('M2B', 2, print_mosfet)
        self.Id_2 = -1
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
        
    def _set(self, Vov_2, Vov_N, Id_2, A1):  # A = -gain
        self.Id_2 = Id_2
        self.A1 = A1
        r1 = self.M2._set(Vov_2, Id_2)
        r2 = self.M2B._set(Vov_N, Id_2)
        r3 = self.M2L._set(1.2*Vov_2*A1, Id_2)
        if r1 == -1 or r2 == -1 or r3 == -1:
            return -1      
        if self.M2L.Vov > 0.67:
            return -2
        
        self.Rin  = 1E+21
        self.Rout = ee.parallel(0.5*self.M2.ro, 1/(self.M2L.gmp))
        self.Cin  = self.M2.cgs + (1+self.A1)*self.M2.cgd
        self.Cout = (1+1/self.A1)*self.M2.cgd + self.M2.cdb + self.M2L.csb + self.M2L.cgs
        return 0
    
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
    r = cs._set(0.2, 0.5, 4E-5, 5)
    print(f'r: {r}')
    cs._print()
    print()

#unit_test_CS()

#%% CD
class CD(mosfet):
    
    def __init__(self):
        print_mosfet = 0
        self.M3   = mosfet('M3',  1, print_mosfet)
        self.M3B  = mosfet('M3B', 2, print_mosfet)
        self.Id_3 = -1
        self.A2   = -1
        self.Rin  = -1
        self.Cin  = -1
        self.Rout = -1
        self.Rout_check = -1
        self.Cout = -1
        self.R_OUT= 10000
        self.C_OUT= 500
        
    def _print(self):
        print('CD stage parameters:')
        ee.print_A('    |A| ', 'V/V', -self.A2)
        ee.print_R('    Rin ', self.Rin)
        ee.print_R('    Rout', self.Rout)
        ee.print_R('    Rout_check', self.Rout_check)
        ee.print_C('    Cin ', self.Cin)
        ee.print_C('    Cout', self.Cout)
        
    def _set(self, Vov_3, Vov_N, Id_3): 
        self.Id_3 = Id_3
        r1 = self.M3._set(Vov_3, Id_3)
        r2 = self.M3B._set(Vov_N, Id_3)
        if r1 == -1 or r2 == -1:
            return -1  
        
        self.Rin  = 1E+21
        self.Rout = ee.parallel( 1 / self.M3.gmp, ee.parallel(0.5*self.M3.ro, self.R_OUT))
        self.Rout_check= ee.parallel(1 / self.M3.gmp, self.R_OUT)
        self.Cin  = self.M3.cgd + (1+self.A2)*self.M3.cgs
        self.Cout = (1+1/self.A2)*self.M3.cgs + self.M3.csb + self.M3B.cgd + self.M3B.cdb + self.C_OUT
        self.A2   = -self.M3.gm*self.Rout
#        self._print()
        return 0
    
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
    cd._set(0.2, 0.6, 1E-4)
    cd._print()
    print('')

#unit_test_CD()

#%% Power Control Module
class PCM():
    I_ref = 1E-5
    Vsup   = 5.
    p_total = 1.8E-3
    p_I_ref = I_ref*Vsup
    
    def __init__(self):
        self.V1        = -1
        self.R_LCG     = -1
        self.Ru        = -1
        self.Rd        = -1
        self.p_Res     = -1
        self.ratio_1   = -1  # ratio Id_1 to Id_3
        self.ratio_2   = -1  # ratio Id_2 to total
        self.Id_half   = -1
        self.p_Id_half = -1
        self.Id_1      = -1
        self.Id_2      = -1
        self.Id_3      = -1
    
    def _print(self):
        print('--- PCM ------------------------')
        print('V1:        %1.3f V'    %self.V1)
        print('R_LCG:     %2.3f kohms' %(1E-3*self.R_LCG))
        print('Ru:        %2.3f kohms' %(1E-3*self.Ru))
        print('Rd:        %2.3f kohms' %(1E-3*self.Rd))
        print('p_Res:     %1.3f mW' %(1E+3*self.p_Res))
        print('ratio_1:   %1.2f' %self.ratio_1)
        print('ratio_2:   %1.2f' %self.ratio_2)
        print('Id_half:   %3.1f uA' %(1E+6*self.Id_half))
        print('p_Id_half: %3.1f mW' %(1E+3*self.p_Id_half))
        print('Id_1:      %3.1f uA' %(1E+6*self.Id_1))
        print('Id_2:      %3.1f uA' %(1E+6*self.Id_2))
        print('Id_3:      %3.1f uA' %(1E+6*self.Id_3))
        print('--------------------------------')
        
    def _set(self, R_LCG, V1, ratio_1, ratio_2):
        self.R_LCG   = R_LCG
        self.V1      = V1
        self.ratio_1 = ratio_1  # ratio Id_1 to Id_3
        self.ratio_2 = ratio_2  # ratio Id_2 to total
        self._upd()
        
    def _upd(self):
        self.Ru   = self.R_LCG*(5/(2.5+self.V1))
        self.Rd   = self.R_LCG*(5/(2.5-self.V1))
        self.p_Res     = (2.5)**2 / self.Rd + (2.5)**2 / self.Rd
        self.p_Id_half = 0.5*(self.p_total - self.p_I_ref) - self.p_Res
        self.Id_half  = self.p_Id_half/self.Vsup
        
        self.Id_1 = self.Id_half * (self.ratio_1*(1-self.ratio_2)) / (1+self.ratio_1)
        self.Id_2 = self.Id_half * self.ratio_2
        self.Id_3 = self.Id_half * (1-self.ratio_2) / (1+self.ratio_1)
#        self._print()
    
    def get_R_LCG(self):
        return self.R_LCG
    def get_V1(self):
        return self.V1
    def get_ratio_1(self):
        return self.ratio_1
    def get_ratio_2(self):
        return self.ratio_2
    def get_Id_1(self):
        return self.Id_1
    def get_Id_2(self):
        return self.Id_2
    def get_Id_3(self):
        return self.Id_3
        
def PCM_unit_test():
    pcm = PCM()
    pcm._set(1.5E+4, 0, 1, 0.2)
    
#PCM_unit_test()

#%% Current Mirror
class CM:
    def __init__(self):
        print_mosfet = 0
        self.Mi1 = mosfet('Mi1', 2, print_mosfet)
        self.Mi2 = mosfet('Mi2', 2, print_mosfet)
        self.Mi3 = mosfet('Mi2', 6, print_mosfet)
        self.Vov_N     = -1
        self.Vov_P     = -1
        self.Id_mirror = -1
        self.mirror_results = []
        
    def _set(self, Vov_N, Vov_P, Id_mirror):
        self.mirror_results = []
#        try:
        self.Vov_N     = Vov_N
        self.Vov_P     = Vov_P
        self.Id_mirror = Id_mirror
        
        # calculate Mi1, Mi2, Mi3 sizing
        self.Mi1._set(self.Vov_N, Id_mirror)
        self.Mi2._set(self.Vov_N, Id_mirror)
        self.Mi3._set(self.Vov_P, Id_mirror)
        
        self.mirror_results.append([self.Mi1._get_W(), self.Mi1._get_L()])
        self.mirror_results.append([self.Mi2._get_W(), self.Mi2._get_L()])
        self.mirror_results.append([self.Mi3._get_W(), self.Mi3._get_L()])

    def _print(self):
        print(f'CM._print() {self.mirror_results}')
        
    def _get(self):
        return self.mirror_results
    
        
def CM_unit_test():
    cm = CM()
    cm._set(0.5, 0.5, 1E-5)
    cm._print()
    results = cm._get()
    print(f'results: {results}')
    

CM_unit_test()
    
    
    
    
    
    
    
    
    
    
    
    