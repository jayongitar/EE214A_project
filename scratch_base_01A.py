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
        self.Vov_L = -1
        self.Vov_B = -1
        # stage parameters: mag, Rin, Rout, Cin, Cout
        self.CG_list = [-1, -1, -1, -1, -1]
        self.CS_list = [-1, -1, -1, -1, -1]
        self.CD_list = [-1, -1, -1, -1, -1]
        # BW:  R, C, P [MHz]
        self.sys1_list = [-1, -1, -1]
        self.sys2_list = [-1, -1, -1]
        self.sys3_list = [-1, -1, -1]
        self.sys4_list = [-1, -1, -1]
        # outputs
        self.power = -1
        self.gain  = -1
        self.f3dB  = -1
        self.f3dB_list = [-1, -1, -1, -1]
        self.FOM   = -1
        
    def _print_input(self):
        print('-'*40)
        print('TIA inputs:')
        print('   Vov_1:   %3.3f V' %self.Vov_1)
        print('   Vov_2:   %3.3f V' %self.Vov_2)
        print('   Vov_3:   %3.3f V' %self.Vov_3)
        
        print('   Vov_N:   %3.3f V' %self.Vov_2)
        print('   Vov_3:   %3.3f V' %self.Vov_3)
        
        print('   R_LCG:   %3.1f kohms' %(1E-3*self.pcm.R_LCG))
        print('   V1:      %3.3f V'    %self.pcm.V1)
        print('   ratio_1: %3.3f' %self.pcm.ratio_1)
        print('   ratio_2: %3.3f' %self.pcm.ratio_2)
    
    def _print_output(self):
        print('TIA outputs:')
#        print('   power: %3.2f mw'  %(1E+3*self.power))
#        print('   gain:  %3.3f kohms' %(1E-3*self.gain))
        print('   f3dB:  %3.1f MHz' %(1E-6*self.f3dB))
        print('   FOM:   %3.3f'   %self.FOM)
        
    def _print_stage_params(self):
        print('-'*40)
        print('CG stage params:')
        ee.print_A('   mag_A0', 'A/V', self.mag_A0)
        ee.print_R('   Rin_CG', self.Rin_CG)
        ee.print_R('   Rout_CG', self.Rout_CG)
        ee.print_C('   Cin_CG', self.Cin_CG)
        ee.print_C('   Cout_CG', self.Cout_CG)
        print('CS stage params:')
        ee.print_A('   mag_A1', 'V/V', self.mag_A1)
        ee.print_R('   Rin_CS', self.Rin_CS)
        ee.print_R('   Rout_CS', self.Rout_CS)
        ee.print_C('   Cin_CS', self.Cin_CS)
        ee.print_C('   Cout_CS', self.Cout_CS)
        print('CD stage params:')
        ee.print_A('   mag_A23', 'V/V', self.mag_A2*self.mag_A3)
        ee.print_R('   Rin_CD', self.Rin_CD)
        ee.print_R('   Rout_CD', self.Rout_CD)
        ee.print_C('   Cin_CD', self.Cin_CD)
        ee.print_C('   Cout_CD', self.Cout_CD)
        print('-'*40)
       
    def _print_poles(self):
        print('-'*40)
        print('TIA poles')
        ee.print_A('   mag_A0', 'V/A', self.mag_A0)
        ee.print_R('   R0', self.R0)
        ee.print_C('   C0', self.C0)
        ee.print_F('   p0', self.p0)
        print()
        ee.print_A('   mag_A1', 'V/A', self.mag_A1)
        ee.print_R('   R1', self.R1)
        ee.print_C('   C1', self.C1)
        ee.print_F('   p1', self.p1)
        print()
        ee.print_A('   mag_A2', 'V/A', self.mag_A2)
        ee.print_R('   R2', self.R2)
        ee.print_C('   C2', self.C2)
        ee.print_F('   p2', self.p2)
        print()
        ee.print_A('   mag_A3', 'V/A', self.mag_A3)
        ee.print_R('   R3', self.R3)
        ee.print_C('   C3', self.C3)
        ee.print_F('   p3', self.p3)
        print('-'*40)
        
      
    def _set(self, Vov_1, Vov_2, Vov_3, R_LCG, V1, ratio_1, ratio_2):
        self.Vov_1 = Vov_1
        self.Vov_2 = Vov_2
        self.Vov_3 = Vov_3
        
        # pcm calculates drain currents from power consumption constrains and 2 input params
        self.pcm._set(R_LCG, V1, ratio_1, ratio_2)
        
        self.mag_A0  = self.cg_LK.get_TI  (self.Vov_1, self.pcm.Id_1,  self.pcm.R_LCG)
        self.Rin_CG  = self.cg_LK.get_Rin (self.Vov_1, self.pcm.Id_1,  self.pcm.R_LCG)
        self.Rout_CG = self.cg_LK.get_Rout(self.Vov_1, self.pcm.Id_1,  self.pcm.R_LCG)
        self.Cin_CG  = self.cg_LK.get_Cin (self.Vov_1, self.pcm.Id_1,  self.pcm.R_LCG)
        self.Cout_CG = self.cg_LK.get_Cout(self.Vov_1, self.pcm.Id_1,  self.pcm.R_LCG)
        
        self.mag_A2  = 0.84
        self.Rin_CD  = self.cd_LK.get_Rin (self.Vov_3, self.pcm.Id_3)
        self.Rout_CD = self.cd_LK.get_Rout(self.Vov_3, self.pcm.Id_3)
        self.Cin_CD  = self.cd_LK.get_Cin (self.Vov_3, self.pcm.Id_3)
        self.Cout_CD = self.cd_LK.get_Cout(self.Vov_3, self.pcm.Id_3)
        self.mag_A3  = self.Rout / (self.Rout + self.Rout_CD)
    
        self.mag_A1  = 40000 / (self.mag_A0*self.mag_A2*self.mag_A3)
        self.Rin_CS  = self.cs_LK.get_Rin (self.Vov_2, self.pcm.Id_2, self.mag_A1)
        self.Rout_CS = self.cs_LK.get_Rout(self.Vov_2, self.pcm.Id_2, self.mag_A1)
        self.Cin_CS  = self.cs_LK.get_Cin (self.Vov_2, self.pcm.Id_2, self.mag_A1)
        self.Cout_CS = self.cs_LK.get_Cout(self.Vov_2, self.pcm.Id_2, self.mag_A1)
        
        self.C0     = self.Cin_CG + self.Cin
        self.R0     = self.Rin_CG
        self.p0     = -1/(6.28*1E-15*self.C0*self.R0)
        
        self.C1     = self.Cout_CG + self.Cin_CS
        self.R1     = ee.parallel(self.Rout_CG, self.Rin_CS)
        self.p1     = -1/(6.28*1E-15*self.C1*self.R1)

        self.C2     = self.Cout_CS + self.Cin_CD
        self.R2     = ee.parallel(self.Rout_CS, self.Rin_CD)
        self.p2     = -1/(6.28*1E-15*self.C2*self.R2)

        self.C3     = self.Cout_CD
        self.R3     = ee.parallel(self.Rout_CD, self.Rout)
        self.p3     = -1/(6.28*1E-15*self.C3*self.R3)
        
        n   = 75
        mag0 = np.zeros(n)
        mag1 = np.zeros(n)
        mag2 = np.zeros(n)
        mag3 = np.zeros(n)
        mag  = np.zeros(n)
        _freq = np.logspace(6, 10, n)
        for i in range(n):
            mag0[i] = abs(self.p0)/np.sqrt(self.p0**2 + _freq[i]**2)
            mag1[i] = abs(self.p1)/np.sqrt(self.p1**2 + _freq[i]**2)
            mag2[i] = abs(self.p2)/np.sqrt(self.p2**2 + _freq[i]**2)
            mag3[i] = abs(self.p3)/np.sqrt(self.p3**2 + _freq[i]**2)
            mag[i]  = mag0[i]*mag1[i]*mag2[i]*mag3[i]
            try:
                mag0_f = interpolate.interp1d(mag0, _freq)
                mag1_f = interpolate.interp1d(mag1, _freq)
                mag2_f = interpolate.interp1d(mag2, _freq)
                mag3_f = interpolate.interp1d(mag3, _freq)
                mag_f  = interpolate.interp1d(mag,  _freq)
                self.f3dB_0 = mag0_f(0.71)
                self.f3dB_1 = mag1_f(0.71)
                self.f3dB_2 = mag2_f(0.71)
                self.f3dB_3 = mag3_f(0.71)
                self.f3dB   = mag_f(0.71)
#                print('freq: %3.3f,   mag: %3.3f' %(1E-6*_freq[i], mag[i]))
            except:
                print('no cutoff freq found')
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
    p = tia._set(0.2, 0.2, 0.2, 2E+4, 0, 1, 0.2)
    print(f'p: {p}')
    if p == 0:
        tia._print_io()
        tia._print_stage_params()
        tia._print_poles()
    else:
        print('invalid input into TIA._set')


#TIA_unit_test()



#%% Sweep

class SWEEP(TIA):
    
    _Vov_1     = np.linspace(0.2, 0.4, 4)
    _Vov_2     = np.linspace(0.2, 0.4, 4)
    _Vov_3     = np.linspace(0.2, 0.4, 4)
    _R_LCG     = np.linspace(1.5E+4,2E+4,4)
    _V1        = np.linspace(-0.2, 0.2, 1)
    _I_ratio_1 = np.linspace(0.7, 1.4, 1)  # ratio_1:  Id_1 to Id_3
    _I_ratio_2 = np.linspace(0.2, 0.3, 1)  # ratio_2:  Id_2 to total
    
    _FOM_Vov_1 = np.zeros(_Vov_1.shape[0])
    _FOM_Vov_2 = np.zeros(_Vov_2.shape[0])
    _FOM_Vov_3 = np.zeros(_Vov_3.shape[0])
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
        max_count = self._Vov_1.shape[0]*self._Vov_2.shape[0]*self._Vov_3.shape[0]*self._I_ratio_1.shape[0]*self._I_ratio_2.shape[0]*self._A1.shape[0]*self._RL.shape[0]
        for i1 in range(self._Vov_1.shape[0]):
            for i2 in range(self._Vov_2.shape[0]):
                for i3 in range(self._Vov_3.shape[0]):
                    for i4 in range(self._R_LCG.shape[0]):
                        for i5 in range(self._V1.shape[0]):
                            for i6 in range(self._I_ratio_1.shape[0]):
                                for i7 in range(self._I_ratio_2.shape[0]):
                                    count+=1
                                    print('*'*50)
#                                    self._print_iteration_input(i1, i2, i3, i4, i5, i6, i7)
                                    self.tia._print_input()
                                    self.tia.pcm._print()
                                    self.tia._set(self._Vov_1[i1],
                                                  self._Vov_2[i2],
                                                  self._Vov_3[i3],
                                                  self._R_LCG[i4],
                                                  self._V1[i5],
                                                  self._I_ratio_1[i6],
                                                  self._I_ratio_2[i7]
                                                  )
                                    self.tia._print_stage_params()
                                    self.tia._print_poles()
                                    self.tia._print_output()
#                                    self.tia._print_io()
                                    
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
                                

        print(f'max_count: {max_count}')
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
        
        fig2 = plt.figure(2)
        fig2.clf()
        fig2.suptitle('Vov_2')
        ax = fig2.add_subplot(111)
        ax.set_xlabel('Vov_2')
        ax.set_ylabel('FOM')
        ax.plot(self._Vov_2, self._FOM_Vov_2)
        ax.grid()
        
        fig3 = plt.figure(3)
        fig3.clf()
        fig3.suptitle('Vov_3')
        ax = fig3.add_subplot(111)
        ax.set_xlabel('Vov_3')
        ax.set_ylabel('FOM')
        ax.plot(self._Vov_3, self._FOM_Vov_3)
        ax.grid()
        
        fig4 = plt.figure(4)
        fig4.clf()
        fig4.suptitle('R_LCG')
        ax = fig4.add_subplot(111)
        ax.set_xlabel('R_LCG')
        ax.set_ylabel('FOM')
        ax.plot(self._R_LCG, self._FOM_R_LCG)
        ax.grid()
        
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





