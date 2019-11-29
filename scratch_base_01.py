import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import time

#%% TIA

class TIA(CG_LK, CS_LK, CD_LK, PCM):
    
    def __init__(self):
        self.cg_LK = CG_LK()
        self.cs_LK = CS_LK()
        self.cd_LK = CD_LK()
        self.pcm = PCM()
        
        # inputs
        self.Vov_1     = -1
        self.Vov_2     = -1
        self.Vov_3     = -1
        
        # stage parameters
        self.Cin     = 100E-15
        
        self.mag_A0  = -1
        self.Rin_CG  = -1
        self.Rout_CG = -1
        self.Cin_CG  = -1
        self.Cout_CG = -1
        
        self.mag_A1  = -1
        self.Rin_CS  = -1
        self.Rout_CS = -1
        self.Cin_CS  = -1
        self.Cout_CS = -1
        
        self.mag_A2  = -1
        self.Rin_CD  = -1
        self.Rout_CD = -1
        self.Cin_CD  = -1
        self.Cout_CD = -1
        
        self.Rout    = 1E+4
        self.Cout    = 500E-15
        self.mag_A3  = -1
        
        # gain and pole location
        self.mag_A0 = -1
        self.C0     = -1
        self.R0     = -1
        self.p0     = -1
        self.mag_A1 = -1
        self.C1     = -1
        self.R1     = -1
        self.p1     = -1
        self.mag_A2 = -1
        self.C2     = -1
        self.R2     = -1
        self.p2     = -1
        self.mag_A3 = -1
        self.C3     = -1
        self.R3     = -1
        self.p3     = -1
        
        # outputs
        self.power = -1
        self.gain  = -1
        self.BW    = -1
        self.FOM   = -1
        
    def _print_io(self):
        print('\nTIA inputs:')
        print('Vov_1:   %3.3f V' %self.Vov_1)
        print('Vov_2:   %3.3f V' %self.Vov_2)
        print('Vov_3:   %3.3f V' %self.Vov_3)
        print('R_LCG:   %3.1f kohms' %(1E-3*self.pcm.R_LCG))
        print('V1:      %3.3f V'    %self.pcm.V1)
        print('ratio_1: %3.3f' %self.pcm.ratio_1)
        print('ratio_2: %3.3f' %self.pcm.ratio_2)
        
        print('TIA outputs:')
        print('power: %3.2f mw'  %(1E+3*self.power))
        print('gain:  %3.3f kohms' %(1E-3*self.gain))
        print('BW:    %3.1f MHz' %(1E-6*self.BW))
        print('FOM:   %3.3f V'   %self.FOM)
        
    def _print_stage_params(self):
        ee.print('Cin', self.Cin)
        ee.print('mag_A0', self.mag_A0)
        
        self.Rin_CG  = -1
        self.Rout_CG = -1
        self.Cin_CG  = -1
        self.Cout_CG = -1
        
        self.mag_A1  = -1
        self.Rin_CS  = -1
        self.Rout_CS = -1
        self.Cin_CS  = -1
        self.Cout_CS = -1
        
        self.mag_A2  = -1
        self.Rin_CD  = -1
        self.Rout_CD = -1
        self.Cin_CD  = -1
        self.Cout_CD = -1
        
        self.Rout    = 1E+4
        self.Cout    = 500E-15
        self.mag_A3  = -1
       
    def _set(self, Vov_1, Vov_2, Vov_3, R_LCG, V1, ratio_1, ratio_2):
        self.Vov_1 = Vov_1
        self.Vov_2 = Vov_2
        self.Vov_3 = Vov_3
        
        # pcm calculates drain currents from power consumption constrains and 2 input params
        self.pcm._set(R_LCG, V1, ratio_1, ratio_2)
        # now the gain of the CG and CD stages are fixed and the gain (A1) for the CS stage can
        # be calculated from the total gain requirement of 40kohms TI
        try: 
            self.mag_A0  = self.cg_LK.get_TI  (self.Vov_1, self.pcm.Id_1, self.pcm.R_LCG)
            self.Rin_CG  = self.cg_LK.get_Rin (self.Vov_1, self.pcm.Id_1, self.pcm.R_LCG)
            self.Rout_CG = self.cg_LK.get_Rout(self.Vov_1, self.pcm.Id_1, self.pcm.R_LCG)
            self.Cin_CG  = self.cg_LK.get_Cin (self.Vov_1, self.pcm.Id_1, self.pcm.R_LCG)
            self.Cout_CG = self.cg_LK.get_Cout(self.Vov_1, self.pcm.Id_1, self.pcm.R_LCG)
        except:
            print('Value in CG out of range...')
        try:
            self.mag_A2  = 0.84
            self.Rin_CD  = self.cd_LK.get_Rin (self.Vov_3, self.pcm.Id_1)
            self.Rout_CD = self.cd_LK.get_Rout(self.Vov_3, self.pcm.Id_1)
            self.Cin_CD  = self.cd_LK.get_Cin (self.Vov_3, self.pcm.Id_1)
            self.Cout_CD = self.cd_LK.get_Cout(self.Vov_3, self.pcm.Id_1)
            self.mag_A3  = self.Rout / (self.Rout + self.Rout_CD)
        except:
            print('Value in CD out of range...')
        try:
            self.mag_A1      = 40000 / (self.mag_A0*self.mag_A2*self.mag_A3)
            self.Rin_CS  = self.cs_LK.get_Rin (self.Vov_2, self.pcm.Id_2, self.A1)
            self.Rout_CS = self.cs_LK.get_Rout(self.Vov_2, self.pcm.Id_2, self.A1)
            self.Cin_CS  = self.cs_LK.get_Cin (self.Vov_2, self.pcm.Id_2, self.A1)
            self.Cout_CS = self.cs_LK.get_Cout(self.Vov_2, self.pcm.Id_2, self.A1)
        except:
            print('Value in CS out of range...')
        
        
def TIA_unit_test():
    tia = TIA()
    tia._set(0.2, 0.2, 0.2, 2E+4, 0, 1, 0.2)
    tia._print_io()


TIA_unit_test()



#%% Sweep

class SWEEP(TIA):
    
    _Vov_1     = np.linspace(0.2, 0.3, 2)
    _Vov_2     = np.linspace(0.2, 0.3, 2)
    _Vov_3     = np.linspace(0.2, 0.3, 2)
    _I_ratio_1 = np.linspace(0.7, 1.4, 2)
    _I_ratio_2 = np.linspace(0.8, 1.25,2)
    _A1        = np.linspace(0.2, 0.3, 2)
    _RL        = np.linspace(2E+4,3E+4,2)
    
    def __init__(self):
        self.CG_LK = CG_LK()
        self.CS_LK = CS_LK()
        self.CD_LK = CD_LK()
        
        # variables
        self.Vov_1 = -1
        self.Id_1  = -1
        self.RL    = -1
        
        self.Vov_2 = -1
        self.Id_2  = -1
        self.A1    = -1
        
        self.Vov_3 = -1
        self.Id_3  = -1 
        
        # parameters
        self.I_ratio_1 = -1
        self.I_ratio_2 = -1
        
    
    def iterate(self):
        start_ms = int(round(time.time() * 1000))
        count = -1
        max_count = self._Vov_1.shape[0]*self._Vov_2.shape[0]*self._Vov_3.shape[0]*self._I_ratio_1.shape[0]*self._I_ratio_2.shape[0]*self._A1.shape[0]*self._RL.shape[0]
        for ind_Vov_1 in range(self._Vov_1.shape[0]):
            for ind_Vov_2 in range(self._Vov_2.shape[0]):
                for ind_Vov_3 in range(self._Vov_3.shape[0]):
                    for ind_I_ratio_1 in range(self._I_ratio_1.shape[0]):
                        for ind_I_ratio_2 in range(self._I_ratio_2.shape[0]):
                            for ind_A1 in range(self._A1.shape[0]):
                                for ind_RL in range(self._RL.shape[0]):
                                    count+=1
                                    print(f'count: {count}')
                                
                                
                               
                                
                                

        print(f'max_count: {max_count}')
        end_ms = int(round(time.time() * 1000))
        print(f'total time: {end_ms-start_ms} ms')
    
def SWEEP_unit_test():
    sweep = SWEEP()
    sweep.iterate()

#SWEEP_unit_test()





