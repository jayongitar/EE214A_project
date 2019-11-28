import numpy as np
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
    
    def __init__(self, mosfet_type):  # 0=nmos l=1, 1=nmos l=2, 2=pmos l=1, 3=pmos l=2
        # set variables
        self.Vov = 0.2
        self.Id  = 1E-5
        self.type= mosfet_type
        self.type_dict = {0:'nmos l=1u', 1:'nmos l=2u', 2:'pmos l=1u', 3:'pmos l=2u',}
        # derived variables
        self.WL  = 3   
        self.L   = 1  # um
        self.W   = 2  # um
        self.gm  = 0.001 # A/V
        self.cgs = 0
        self.cgd = 0
        self.cgb = 0
        self.csb = 0
        self.cdb = 0
        
    def print_basic(self):
        print(f'type: {self.type_dict[self.type]}')
        print('input variables: Vov = %2.2fV, Id =%3.0fuA' %(self.Vov, 1E+6*self.Id) )
        print('calc parameters: \nWL  = %2.2f,  W = %1.1fum,  L  = %1.0fum,  gm = %1.2fmA/V' 
              %(self.WL, self.W, self.L, 1000*self.gm) )
    
    def print_caps(self):
        print('cgs: %3.1f fF' %self.cgs)
        print('cgd: %3.1f fF' %self.cgd)
        print('cgb: %3.1f fF' %self.cgb)
        print('csb: %3.1f fF' %self.csb)
        print('cdb: %3.1f fF' %self.cdb)
        
    def print_all(self):
        self.print_basic()
        self.print_caps()
        
    def set_params(self, Vov, Id):
        self.Vov = Vov
        self.Id  = Id*1E-6
        # WL = 2*Id / unCox*Vov^2
        if self.type == 0 or self.type == 1:
            self.WL = 2*self.Id / ( self.unCox*self.Vov**2 )
            self.W  = self.WL*(self.type+1)
        if self.type == 2 or self.type == 3:
            self.WL = 2*self.Id / ( self.upCox*self.Vov**2 )
            self.W  = self.WL*(self.type-1)
        # gm = 2*Id/Vov
        self.gm = 2*self.Id/self.Vov
        
        if self.type == 0 or self.type == 1: 
            self.cgs = self.cgs_lookup_1u(self.WL)
            self.cgd = self.cgd_lookup_1u(self.WL)
            self.cgb = self.cgtot_lookup_1u(self.WL) - self.cgs - self.cgd
            self.csb = self.cstot_lookup_1u(self.WL) - self.cgs
            self.cdb = self.cdtot_lookup_1u(self.WL) - self.cgd
        if self.type == 0 or self.type == 1: 
            self.cgs = self.cgs_lookup_2u(self.WL)
            self.cgd = self.cgd_lookup_2u(self.WL)
            self.cgb = self.cgtot_lookup_2u(self.WL) - self.cgs - self.cgd
            self.csb = self.cstot_lookup_2u(self.WL) - self.cgs
            self.cdb = self.cdtot_lookup_2u(self.WL) - self.cgd
        

def unit_test_mosfet(test_type):
    nmos_1 = mosfet(0)
    nmos_1.print_all()
    print()
    nmos_1.set_params(0.2, 10)
    nmos_1.print_all()
    print()
    nmos_1.set_params(0.4, 70)
    nmos_1.print_all()
    print()

unit_test_mosfet(1)






























