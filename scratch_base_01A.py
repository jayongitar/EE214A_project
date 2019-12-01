import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import time

#%% TIA

class TIA(CG, CS, CD, PCM):
    
    def __init__(self):
        self.cg = CG()
        self.cs = CS()
        self.cd = CD()
        self.pcm = PCM()
        # inputs
        self.Vov_1 = -1
        self.Vov_2 = -1
        self.Vov_3 = -1
        self.Vov_N = -1
        self.Vov_P = -1
        self.Vov_2L= -1
        # stage parameters: mag, Rin, Rout, Cin, Cout
        self.CG_list = [-1, -1, -1, -1, -1]
        self.CS_list = [-1, -1, -1, -1, -1]
        self.CD_list = [-1, -1, -1, -1, -1]
        # BW:  R, C, P [MHz]
        self.sys1_list = [-1, -1, -1]
        self.sys2_list = [-1, -1, -1]
        self.sys3_list = [-1, -1, -1]
        self.sys4_list = [-1, -1, -1]
        # magnitudes
        self.M_list = [-1, -1, -1]
        # outputs
        self.power = -1
        self.gain  = -1
        self.f3dB  = -1
        self.f3dB_list = [-1, -1, -1, -1]
        self.FOM   = -1
        
    def _print_input(self):
        print('-'*40)
        print('TIA inputs:')
        print('  Vov_1:   %3.3f V' %self.Vov_1)
        print('  Vov_2:   %3.3f V' %self.Vov_2)
        print('  Vov_3:   %3.3f V' %self.Vov_3)
        
        print('  Vov_N:   %3.3f V' %self.Vov_N)
        print('  Vov_P:   %3.3f V' %self.Vov_P)
        
        print('  R_LCG:   %3.1f kohms' %(1E-3*self.pcm.R_LCG))
        print('  V1:      %3.3f V'    %self.pcm.V1)
        print('  ratio_1: %3.3f' %self.pcm.ratio_1)
        print('  ratio_2: %3.3f' %self.pcm.ratio_2)
    
    def _print_output(self):
        print('TIA outputs:')
#        print('   power: %3.2f mw'  %(1E+3*self.power))
#        print('   gain:  %3.3f kohms' %(1E-3*self.gain))
        print('   f3dB:  %3.1f MHz' %(1E-6*self.f3dB))
        print('   FOM:   %3.3f'   %self.FOM)
        
    def _print_stage_params(self):
        print('-'*40)
        print('CG stage params:')
        ee.print_A('   mag_A0 ', 'A/V', self.CG_list[0])
        ee.print_R('   Rin_CG ', self.CG_list[1])
        ee.print_R('   Rout_CG', self.CG_list[2])
        ee.print_C('   Cin_CG ', self.CG_list[3])
        ee.print_C('   Cout_CG', self.CG_list[4])
        print('CS stage params:')
        ee.print_A('   mag_A1 ', 'V/V', self.CS_list[0])
        ee.print_R('   Rin_CS ', self.CS_list[1])
        ee.print_R('   Rout_CS', self.CS_list[2])
        ee.print_C('   Cin_CS' , self.CS_list[3])
        ee.print_C('   Cout_CS', self.CS_list[4])
        print('CD stage params:')
        ee.print_A('   mag_A2 ', 'V/V', self.CD_list[0])
        ee.print_R('   Rin_CD ', self.CD_list[1])
        ee.print_R('   Rout_CD', self.CD_list[2])
        ee.print_C('   Cin_CD ', self.CD_list[3])
        ee.print_C('   Cout_CD', self.CD_list[4])
        print('-'*40)
       
    def _print_poles(self):
        print('-'*40)
        print('TIA poles')
        ee.print_R('   R0', self.sys1_list[0])
        ee.print_C('   C0', self.sys1_list[1])
        ee.print_F('   p0', self.sys1_list[2])
        print()
        ee.print_R('   R1', self.sys2_list[0])
        ee.print_C('   C1', self.sys2_list[1])
        ee.print_F('   p1', self.sys2_list[2])
        print()
        ee.print_R('   R2', self.sys3_list[0])
        ee.print_C('   C2', self.sys3_list[1])
        ee.print_F('   p2', self.sys3_list[2])
        print()
#        ee.print_A('   mag_A3', 'V/A', self.mag_A3)
        ee.print_R('   R3', self.sys4_list[0])
        ee.print_C('   C3', self.sys4_list[1])
        ee.print_F('   p3', self.sys4_list[2])
        print('-'*40)
        
    # _in = (Vov_1, Vov_2, Vov_3, V_BN, V_BP, R_LCG, V1, ratio_1, ratio_2)
    def _set(self, _in):
        # unpack _in to local variables
        Vov_1 =   _in[0]
        Vov_2 =   _in[1]
        Vov_3 =   _in[2]
        V_BN  =   _in[3]
        V_BP  =   _in[4]
        R_LCG =   _in[5]
        V1    =   _in[6]
        ratio_1 = _in[7]  # ratio_1: Id_1 to Id_3
        ratio_2 = _in[8]  # ratio_2: Id_2 to total
        # pcm calculates drain currents from power consumption constrains and 2 input params
        self.pcm._set(R_LCG, V1, ratio_1, ratio_2)
        
        # cg._set(self, Vov_1, V_BN, V_BP, Id_1, R_LCG):
        r1 = self.cg._set(Vov_1, V_BN, V_BP, self.pcm.get_Id_1(), R_LCG)
        self.CG_list = [self.cg.get_TI(), self.cg.get_Rin(), self.cg.get_Rout(), self.cg.get_Cin(), self.cg.get_Cout()]   
        
        # cd._set(self, Vov_3, V_BN, Id_3):
        r3 = self.cd._set(Vov_3, V_BN, self.pcm.get_Id_3())
        self.CD_list = [-self.cd.get_A2(), self.cd.get_Rin(), self.cd.get_Rout(), self.cd.get_Cin(), self.cd.get_Cout()]
        
#        print(f'CG_mag: {self.CG_list[0]}')
#        print(f'CD_mag: {self.CD_list[0]}')
        A1 = 40000 / (self.CG_list[0] + self.CD_list[0])
#        print(f'A1: {A1}')
        # cs._set(self, Vov_2, V_BN, Id_2, A1):
        r2 = self.cs._set(Vov_2, V_BN, self.pcm.get_Id_2(), A1)
        self.CS_list = [A1, self.cs.get_Rin(), self.cs.get_Rout(), self.cs.get_Cin(), self.cs.get_Cout()]
        
        # stage parameters: 0:mag, 1:Rin, 2:Rout, 3:Cin, 4:Cout
        self.M_list = [self.CG_list[0], self.CS_list[0], self.CD_list[0]]
        # CG in
        self.sys1_list[0] = self.CG_list[1]  # Rin_CG
        self.sys1_list[1] = self.CG_list[3]  # Cin_CG
        self.sys1_list[2] = -1/(6.28*1E-15 * self.sys1_list[0] * self.sys1_list[1])
        # CS in, CG_out
        self.sys2_list[0] = ee.parallel(self.CS_list[1], self.CG_list[2])
        self.sys2_list[1] = self.CS_list[3]+ self.CG_list[4]
        self.sys2_list[2] = -1/(6.28*1E-15 * self.sys2_list[0] * self.sys2_list[1])
        # CD in, CS out
        self.sys3_list[0] = ee.parallel(self.CD_list[1], self.CS_list[2])
        self.sys3_list[1] = self.CD_list[3]+ self.CS_list[4]
        self.sys3_list[2] = -1/(6.28*1E-15 * self.sys3_list[0] * self.sys3_list[1])
        # CD out
        self.sys4_list[0] = self.CD_list[2]
        self.sys4_list[1] = self.CD_list[4]
        self.sys4_list[2] = -1/(6.28*1E-15 * self.sys4_list[0] * self.sys4_list[1])
        
        print(f'r1: {r1}')
        print(f'r2: {r2}')
        print(f'r3: {r3}')
        if r1 == -1 or r2 == -1 or r3 == -1:
            print('invalid _in to tia._set')
            return -1
        print(f'M: {self.M_list}')
        print('Done _setting')
        
        n   = 75
        mag_H0 = np.zeros(n)
        mag_H1 = np.zeros(n)
        mag_H2 = np.zeros(n)
        mag_H3 = np.zeros(n)
        mag_H  = np.zeros(n)
        _freq = np.logspace(6, 10, n)
        for i in range(n):
            mag_H0[i] = abs(self.sys1_list[2])  /np.sqrt(self.sys1_list[2]**2 + _freq[i]**2)
            mag_H1[i] = abs(self.sys2_list[2]) / np.sqrt(self.sys2_list[2]**2 + _freq[i]**2)
            mag_H2[i] = abs(self.sys3_list[2]) / np.sqrt(self.sys3_list[2]**2 + _freq[i]**2)
            mag_H3[i] = abs(self.sys4_list[2]) / np.sqrt(self.sys4_list[2]**2 + _freq[i]**2)
            mag_H [i] = mag_H0[i]*mag_H1[i]*mag_H2[i]*mag_H3[i]
            try:
                mag0_f = interpolate.interp1d(mag_H0, _freq)
                mag1_f = interpolate.interp1d(mag_H1, _freq)
                mag2_f = interpolate.interp1d(mag_H2, _freq)
                mag3_f = interpolate.interp1d(mag_H3, _freq)
                mag_f  = interpolate.interp1d(mag_H,  _freq)
                self.f3dB_0 = mag0_f(0.71)
                self.f3dB_1 = mag1_f(0.71)
                self.f3dB_2 = mag2_f(0.71)
                self.f3dB_3 = mag3_f(0.71)
                self.f3dB   = mag_f (0.71)
#                print('freq: %3.3f,   mag: %3.3f' %(1E-6*_freq[i], mag[i]))
            except:
                print('no cutoff freq found')
                return -2
#        print('f3dB_0: %3.3f MHz' %(1E-6*self.f3dB_0))
#        print('f3dB_1: %3.3f MHz' %(1E-6*self.f3dB_1))
#        print('f3dB_2: %3.3f MHz' %(1E-6*self.f3dB_2))
#        print('f3dB_3: %3.3f MHz' %(1E-6*self.f3dB_3))
#        print('f3dB:   %3.3f MHz' %(1E-6*self.f3dB))
        self.FOM = 40*(1E-6*self.f3dB)/2
        return 0
    
p = 0     
def TIA_unit_test():
    tia = TIA()
    # _in = (Vov_1, Vov_2, Vov_3, V_BN, V_BP, R_LCG, V1, ratio_1, ratio_2)
    _in = [0.2, 0.2, 0.2, 0.5, 0.5, 2E+4, 0, 1, 0.2]
    r = tia._set(_in)
    print(f'tia ret: {r}')
#    if p == 0:
#        tia._print_input()
#        tia._print_stage_params()
#        tia._print_poles()
#    else:
#        print('invalid input into TIA._set')


#TIA_unit_test()



#%% Sweep

class SWEEP(TIA):
    
    _Vov_1     = np.linspace(0.2, 0.4, 1)
    _Vov_2     = np.linspace(0.2, 0.4, 1)
    _Vov_3     = np.linspace(0.2, 0.4, 1)
    _V_BN      = np.linspace(0.3, 0.5, 1)
    _V_BP      = np.linspace(0.3, 0.5, 1)
    _R_LCG     = np.linspace(1.5E+4,2E+4,1)
    _V1        = np.linspace(-0.2, 0.2, 1)
    _I_ratio_1 = np.linspace(0.7, 1.4, 1)  # ratio_1:  Id_1 to Id_3
    _I_ratio_2 = np.linspace(0.2, 0.3, 1)  # ratio_2:  Id_2 to total
    
    _FOM_Vov_1 = np.zeros(_Vov_1.shape[0])
    _FOM_Vov_2 = np.zeros(_Vov_2.shape[0])
    _FOM_Vov_3 = np.zeros(_Vov_3.shape[0])
    _FOM_V_BN  = np.zeros(_V_BN.shape[0])
    _FOM_V_BP  = np.zeros(_V_BP.shape[0])
    _FOM_R_LCG = np.zeros(_R_LCG.shape[0])
    _FOM_V1    = np.zeros(_V1.shape[0])
    _FOM_I_ratio_1 = np.zeros(_I_ratio_1.shape[0])
    _FOM_I_ratio_2 = np.zeros(_I_ratio_2.shape[0])
    
   
    def __init__(self):
        self.tia = TIA()
        print('initializing sweep')
#        self._print_input()
        
    def _print_input(self):
        print('-'*40)
        print('Inputs')
        print(f'_Vov_1:     {self._Vov_1}')
        print(f'_Vov_2:     {self._Vov_2}')
        print(f'_Vov_3:     {self._Vov_3}')
        print(f'_V_BN:      {self._V_BN}')
        print(f'_V_BP:      {self._V_BP}')
        print(f'_R_LCG:     {self._R_LCG}')
        print(f'_V1:        {self._V1}')
        print(f'_I_ratio_1: {self._I_ratio_1}')
        print(f'_I_ratio_2: {self._I_ratio_2}')
        print('-'*40)
        
    def _print_iteration_input(self, i1, i2, i3, i4, i5, i6, i7):
        print('-'*40)
        print('Iteration Inputs')
        print(f'_Vov_1:     {self._Vov_1[i1]}')
        print(f'_Vov_2:     {self._Vov_2[i2]}')
        print(f'_Vov_3:     {self._Vov_3[i3]}')
        print(f'_R_LCG:     {self._R_LCG[i4]}')
        print(f'_V1:        {self._V1[i5]}')
        print(f'_I_ratio_1: {self._I_ratio_1[i6]}')
        print(f'_I_ratio_2: {self._I_ratio_2[i7]}')
    
    def iterate(self):
        start_ms = int(round(time.time() * 1000))
        count = -1
#        max_count = self._Vov_1.shape[0]*self._Vov_2.shape[0]*self._Vov_3.shape[0]*self._I_ratio_1.shape[0]*self._I_ratio_2.shape[0]*self._A1.shape[0]*self._RL.shape[0]
        for i1 in range(self._Vov_1.shape[0]):
            for i2 in range(self._Vov_2.shape[0]):
                for i3 in range(self._Vov_3.shape[0]):
                    
                    for i4 in range(self._V_BN.shape[0]):
                        for i5 in range(self._V_BP.shape[0]):
                    
                            for i6 in range(self._R_LCG.shape[0]):
                                for i7 in range(self._V1.shape[0]):
                                    for i8 in range(self._I_ratio_1.shape[0]):
                                        for i9 in range(self._I_ratio_2.shape[0]):
                                            count+=1
                                            print('*'*50)
#                                           self._print_iteration_input(i1, i2, i3, i4, i5, i6, i7)
                                            self.tia._print_input()
                                            self.tia.pcm._print()
                                            _in = (self._Vov_1[i1],
                                                   self._Vov_2[i2],
                                                   self._Vov_3[i3],
                                                   self._V_BN [i4],
                                                   self._V_BP [i5],
                                                   self._R_LCG[i6],
                                                   self._V1   [i7],
                                                   self._I_ratio_1[i8],
                                                   self._I_ratio_2[i9]
                                                   )
                                            self.tia._set(_in)
                                            self.tia._print_stage_params()
                                            self.tia._print_poles()
                                            self.tia._print_output()
#                                           self.tia._print_io()
                                    
                                            if self.tia.FOM > self._FOM_Vov_1[i1]:
                                                self._FOM_Vov_1[i1] = self.tia.FOM
                                            if self.tia.FOM > self._FOM_Vov_2[i2]:
                                                self._FOM_Vov_2[i2] = self.tia.FOM
                                            if self.tia.FOM > self._FOM_Vov_3[i3]:
                                                self._FOM_Vov_3[i3] = self.tia.FOM                                          
                                            if self.tia.FOM > self._FOM_R_LCG[i4]:
                                                self._FOM_R_LCG[i4] = self.tia.FOM
                                            if self.tia.FOM > self._FOM_V1[i5]:
                                                self._FOM_V1[i5] = self.tia.FOM
                                            if self.tia.FOM > self._FOM_I_ratio_1[i6]:
                                                self._FOM_I_ratio_1[i6] = self.tia.FOM
                                            if self.tia.FOM > self._FOM_I_ratio_2[i7]:
                                                self._FOM_I_ratio_2[i7] = self.tia.FOM         
                                

#        print(f'max_count: {max_count}')
        end_ms = int(round(time.time() * 1000))
        print(f'total time: {end_ms-start_ms} ms')
        
    def _plot(self):
#        mm_Vov_1 = abs(self._FOM_Vov_1[0] - self._FOM_Vov_1[-1])
#        mm_Vov_2 = abs(self._FOM_Vov_2[0] - self._FOM_Vov_2[-1])
#        mm_Vov_3 = abs(self._FOM_Vov_3[0] - self._FOM_Vov_3[-1])
#        mm_R_LCG = abs(self._FOM_R_LCG[0] - self._FOM_R_LCG[-1])
#        mm_V_1   = abs(self._FOM_V_1[0]   - self._FOM_v_1[-1])
#        mm_I_ratio_1 = abs(self._FOM_I_ratio_1[0] - self._FOM_I_ratio_1[-1])
#        mm_I_ratio_2 = abs(self._FOM_I_ratio_2[0] - self._FOM_I_ratio_2[-1])
#        mm = max([mm_Vov_1, mm_Vov_2, mm_Vov_3, mm_R_LCG, mm_V_1, mm_I_ratio_1, mm_I_ratio_2])
        
        fig1 = plt.figure(1)
        fig1.clf()
        fig1.suptitle('Vov_1')
        ax = fig1.add_subplot(111)
        ax.set_xlabel('Vov_1')
        ax.set_ylabel('FOM')
        ax.plot(self._Vov_1, self._FOM_Vov_1)
        ax.grid()
        
#        fig2 = plt.figure(2)
#        fig2.clf()
#        fig2.suptitle('Vov_2')
#        ax = fig2.add_subplot(111)
#        ax.set_xlabel('Vov_2')
#        ax.set_ylabel('FOM')
#        ax.plot(self._Vov_2, self._FOM_Vov_2)
#        ax.grid()
#        
#        fig3 = plt.figure(3)
#        fig3.clf()
#        fig3.suptitle('Vov_3')
#        ax = fig3.add_subplot(111)
#        ax.set_xlabel('Vov_3')
#        ax.set_ylabel('FOM')
#        ax.plot(self._Vov_3, self._FOM_Vov_3)
#        ax.grid()
#        
#        fig4 = plt.figure(4)
#        fig4.clf()
#        fig4.suptitle('R_LCG')
#        ax = fig4.add_subplot(111)
#        ax.set_xlabel('R_LCG')
#        ax.set_ylabel('FOM')
#        ax.plot(self._R_LCG, self._FOM_R_LCG)
#        ax.grid()
        
#        fig5 = plt.figure(5)
#        fig5.clf()
#        fig5.suptitle('V1')
#        ax = fig5.add_subplot(111)
#        ax.set_xlabel('V1')
#        ax.set_ylabel('FOM')
#        ax.plot(self._V1, self._FOM_V1)
#        ax.grid()        
#        
#        fig6 = plt.figure(6)
#        fig6.clf()
#        fig6.suptitle('I_ratio_1')
#        ax = fig6.add_subplot(111)
#        ax.set_xlabel('I_ratio_1')
#        ax.set_ylabel('FOM')
#        ax.plot(self._I_ratio_1, self._FOM_I_ratio_1)
#        ax.grid()
#        
#        fig7 = plt.figure(7)
#        fig7.clf()
#        fig7.suptitle('I_ratio_2')
#        ax = fig7.add_subplot(111)
#        ax.set_xlabel('_I_ratio_2')
#        ax.set_ylabel('FOM')
#        ax.plot(self._I_ratio_2, self._FOM_I_ratio_2)
#        ax.grid()
        
       
    
def SWEEP_unit_test():
    sweep = SWEEP()
    sweep.iterate()
    sweep._plot()

SWEEP_unit_test()





