import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import time
#np.set_printoptions(precision=2)
#np.set_printoptions(linewidth=1000)
#np.set_printoptions(threshold=np.inf)
#np.set_printoptions(suppress=False,
#   formatter={'float_kind':'{:2.2f}'.format})

np.set_printoptions(suppress=True,
   formatter={'float_kind':'{:7.2f}'.format}, linewidth=130)

#%% static inline

class out:
    
    @staticmethod
    def _print_design(design):
#        design = self.tia._get_design()
        print('              Vov_1   Vov_2   Vov_3   Vov_N   Vov_P   R_LCG    V1   ratio_1 ratio_2   Ru      Rd     FOM     f3dB [MHz]')
        print(f'design:    {design}')
        
    @staticmethod
    def _print_mosfets(mosfets):
#        mosfets = self.tia._get_mosfets()
        
        labels_list = ['Vov [V]   ', 
                       'Id [uA]   ', 
                       'type      ', 
                       'WL        ', 
                       'W [um]    ', 
                       'L[um]     ', 
                       'gm [mA/V] ', 
                       'gmp [mA/V]', 
                       'ro [kohms]', 
                       'cgs [fF]  ', 
                       'cgd [fF]  ', 
                       'csb [fF]  ', 
                       'cdb [fF]  ']
        
        print('mosfets:      M1L     M1      M1B     M2L     M2      M2B     M3      M3B     Mi1     Mi2     Mi3')
        for i in range(13):
            print(f'{labels_list[i]} {mosfets[i,:]}')

#%% TIA

class TIA(CG, CS, CD, PCM, CM):
    
    def __init__(self):
        self.cg = CG()
        self.cs = CS()
        self.cd = CD()
        self.pcm = PCM()
        self.cm = CM()
        # inputs
        self.Vov_1 = -1
        self.Vov_2 = -1
        self.Vov_3 = -1
        self.Vov_N = -1
        self.Vov_P  = -1
        # stage parameters: mag, Rin, Rout, Cin, Cout
        self.CG_list = [-1, -1, -1, -1, -1]
        self.CS_list = [-1, -1, -1, -1, -1]
        self.CD_list = [-1, -1, -1, -1, -1]
        # BW:  R, C, P [MHz]
        self.sys1_list = [-1, -1, -1]
        self.sys2_list = [-1, -1, -1]
        self.sys3_list = [-1, -1, -1]
        self.sys4_list = [-1, -1, -1]
        self.M_list = [-1, -1, -1]
        # outputs
        self.power = -1
        self.gain  = -1
        self.f3dB_list = [-1, -1, -1, -1]
        self.FOM   = -1
        self.performance_list = []
    
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
        
    # _in = (Vov_1, Vov_2, Vov_3, Vov_N, Vov_P, R_LCG, V1, ratio_1, ratio_2)
    def _set(self, _in):
        # unpack _in
        self.Vov_1       = _in[0]
        self.Vov_2       = _in[1]
        self.Vov_3       = _in[2]
        self.Vov_N       = _in[3]
        self.Vov_P       = _in[4]
        self.pcm._set(_in[5], _in[6], _in[7], _in[8])
        # _set(self, Vov_N, Vov_P, Id_mirror)
        self.cm._set(self.Vov_N, self.Vov_P, 1E-5)
        
        # cg._set(self, Vov_1, Vov_N, Vov_P, Id_1, R_LCG):
        r1 = self.cg._set(self.Vov_1, self.Vov_N, self.Vov_P, self.pcm.get_Id_1(), self.pcm.get_R_LCG())
        self.CG_list = [self.cg.get_TI(), self.cg.get_Rin(), self.cg.get_Rout(), self.cg.get_Cin(), self.cg.get_Cout()]   
        # cd._set(self, Vov_3, Vov_N, Id_3):
        r3 = self.cd._set(self.Vov_3, self.Vov_N, self.pcm.get_Id_3())
        self.CD_list = [-self.cd.get_A2(), self.cd.get_Rin(), self.cd.get_Rout(), self.cd.get_Cin(), self.cd.get_Cout()]
#        print(f'CG_mag: {self.CG_list[0]}')
#        print(f'CD_mag: {self.CD_list[0]}')
        A1 = 40000 / (self.CG_list[0]*self.CD_list[0])
#        print(f'A1: {A1}')
        # cs._set(self, Vov_2, Vov_N, Id_2, A1):
        r2 = self.cs._set(self.Vov_2, self.Vov_N, self.pcm.get_Id_2(), A1)
        self.CS_list = [-A1, self.cs.get_Rin(), self.cs.get_Rout(), self.cs.get_Cin(), self.cs.get_Cout()]
        
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
        
        
        if r1 == -1 or r2 == -1 or r3 == -1:
            print('invalid _in to tia._set')
            return -1
#        print(f'r1: {r1}')
#        print(f'r2: {r2}')
#        print(f'r3: {r3}')
#        print(f'M: {self.M_list}')
#        print('Done _setting')
        
        n   = 75
        mag1   = np.zeros(n)
        mag2   = np.zeros(n)
        mag3   = np.zeros(n)
        mag4   = np.zeros(n)
        mag_H1 = np.zeros(n)
        mag_H2 = np.zeros(n)
        mag_H3 = np.zeros(n)
        mag_H4 = np.zeros(n)
        
        _freq = np.logspace(6, 10, n)
        for i in range(n):
            mag1[i] = abs(self.sys1_list[2])  /np.sqrt(self.sys1_list[2]**2 + _freq[i]**2)
            mag2[i] = abs(self.sys2_list[2]) / np.sqrt(self.sys2_list[2]**2 + _freq[i]**2)
            mag3[i] = abs(self.sys3_list[2]) / np.sqrt(self.sys3_list[2]**2 + _freq[i]**2)
            mag4[i] = abs(self.sys4_list[2]) / np.sqrt(self.sys4_list[2]**2 + _freq[i]**2)
            mag_H1[i] = mag1[i]
            mag_H2[i] = mag1[i]*mag2[i]
            mag_H3[i] = mag1[i]*mag2[i]*mag3[i]
            mag_H4[i] = mag1[i]*mag2[i]*mag3[i]*mag4[i]
        
        try:
            mag1_f = interpolate.interp1d(mag_H1,  _freq)
            mag2_f = interpolate.interp1d(mag_H2, _freq)
            mag3_f = interpolate.interp1d(mag_H3, _freq)
            mag4_f = interpolate.interp1d(mag_H4, _freq)
            self.f3dB_list[0] = mag1_f(0.71)
            self.f3dB_list[1] = mag2_f(0.71)
            self.f3dB_list[2] = mag3_f(0.71)
            self.f3dB_list[3] = mag4_f(0.71)
#            print('f3dB1: %3.3f MHz' %(1E-6*self.f3dB_list[0]))
#            print('f3dB2: %3.3f MHz' %(1E-6*self.f3dB_list[1]))
#            print('f3dB3: %3.3f MHz' %(1E-6*self.f3dB_list[2]))
#            print('f3dB4: %3.3f MHz' %(1E-6*self.f3dB_list[3]))

        except:
            print('no cutoff freq found')
            return -2
        self.FOM = 40*(1E-6*self.f3dB_list[3])/2
        return 0
    
    def _get_design(self):
        design = np.zeros((13))
        inputs = self._get_input()
        performance = self._get_performance()
        design[0:11] = inputs
        design[11]   = performance[0]
        design[12]  = performance[3]
        return design  

    def _get_input(self):        
        Vov_1 = np.round(self.Vov_1, 3)
        Vov_2 = np.round(self.Vov_2, 3)
        Vov_3 = np.round(self.Vov_3, 3)
        Vov_N = np.round(self.Vov_N, 3)
        Vov_P = np.round(self.Vov_P, 3)
        R_LCG   = np.round(self.pcm.R_LCG/1000, 1)
        V1      = np.round(self.pcm.V1, 3)
        ratio_1 = np.round(self.pcm.ratio_1, 3)
        ratio_2 = np.round(self.pcm.ratio_2, 3)
        Ru      = np.round(self.pcm.Ru/1000, 2)
        Rd      = np.round(self.pcm.Rd/1000, 2)
    
        return_array = np.array((Vov_1, Vov_2, Vov_3, Vov_N, Vov_P, R_LCG, V1, ratio_1, ratio_2, Ru, Rd))
        return return_array
    
    def _print_performance(self):
        performance = self._get_performance()
        print('                FOM   f3dB1   f3dB2   f3dB3   f3dB4 [MHz]')
        print(f'tia perfm: {performance}')
    
    def _get_performance(self):        
        FOM   = np.round(self.FOM, 0)
        f3dB1 = np.round(self.f3dB_list[0]/1E+6, 2)
        f3dB2 = np.round(self.f3dB_list[1]/1E+6, 2)
        f3dB3 = np.round(self.f3dB_list[2]/1E+6, 2)
        f3dB4 = np.round(self.f3dB_list[3]/1E+6, 2)
        return_array = np.array((FOM, f3dB1, f3dB2, f3dB3, f3dB4))
        return return_array
    
    def _get_mosfets(self):
        M1L = self.cg.M1L._get()
        M1  = self.cg.M1._get()
        M1B = self.cg.M1B._get()
    
        M2L = self.cs.M2L._get()
        M2  = self.cs.M2._get()
        M2B = self.cs.M2B._get()
        
        M3  = self.cd.M3._get()
        M3B = self.cd.M3B._get()
        
        Mi1 = self.cm.Mi1._get()
        Mi2 = self.cm.Mi2._get()
        Mi3 = self.cm.Mi3._get()  
        
        mosfets = np.zeros((13, 11))
        for i in range(13):
            mosfets[i, 0] = M1L[i]
            mosfets[i, 1] = M1 [i]
            mosfets[i, 2] = M1B[i]
            
            mosfets[i, 3] = M2L[i]
            mosfets[i, 4] = M2 [i]
            mosfets[i, 5] = M2B[i]
            
            mosfets[i, 6] = M3 [i]
            mosfets[i, 7] = M3B[i]
            
            mosfets[i, 8 ] = Mi1[i]
            mosfets[i, 9 ] = Mi2[i]
            mosfets[i, 10] = Mi3[i]
        return mosfets
   
    
p = 0     
def TIA_unit_test():
    tia = TIA()
    
    _in = [0.2, 0.2, 0.2, 0.5, 0.5, 2E+4, 0.0, 1, 0.2]
    r = tia._set(_in)
    print('-'*100)
    print(f'tia ret: {r}')

    design = tia._get_design()
    out._print_design(design)
    print()
    
    mosfets = tia._get_mosfets()
    out._print_mosfets(mosfets)
    print()
    print('-'*100)
    
#TIA_unit_test()



#%% Sweep

class SWEEP(TIA):
    
    _Vov_1   = np.linspace(0.2,   0.4,   4)
    _Vov_2   = np.linspace(0.25,   0.25, 1)
    _Vov_3   = np.linspace(0.3,   0.4,   4)
    _Vov_N   = np.linspace(0.6,   0.7,   1)
    _Vov_P   = np.linspace(0.6,   0.7,   1)
    _R_LCG   = np.linspace(1.5E+4,2.0E+4,7)
    _V1      = np.linspace(0.6,   0.7,   1)
    _ratio_1 = np.linspace(0.3,   1.0,   7)  # ratio_1:  Id_1 to Id_3
    _ratio_2 = np.linspace(0.15,   0.3,  7)  # ratio_2:  Id_2 to total
    
    _FOM_Vov_1   = np.zeros(_Vov_1.shape[0])
    _FOM_Vov_2   = np.zeros(_Vov_2.shape[0])
    _FOM_Vov_3   = np.zeros(_Vov_3.shape[0])
    _FOM_Vov_N   = np.zeros(_Vov_N.shape[0])
    _FOM_Vov_P   = np.zeros(_Vov_P.shape[0])
    _FOM_R_LCG   = np.zeros(_R_LCG.shape[0])
    _FOM_V1      = np.zeros(_V1.shape[0])
    _FOM_ratio_1 = np.zeros(_ratio_1.shape[0])
    _FOM_ratio_2 = np.zeros(_ratio_2.shape[0])
    
    total_iterations =         \
            _Vov_1.shape[0]*   \
            _Vov_2.shape[0]*   \
            _Vov_3.shape[0]*   \
            _Vov_N.shape[0]*   \
            _Vov_P.shape[0]*   \
            _R_LCG.shape[0]*   \
            _V1.shape[0]*      \
            _ratio_1.shape[0]* \
            _ratio_2.shape[0]
                        
    print(f'total iterations: {total_iterations}')
                        
                        
    
   
    def __init__(self):
        self.tia = TIA()
        self.valid_point_list    = []
        self.invalid_point_list  = []
        self.results_list        = []
        self.sorted_results_list = []
        self.top_results_list    = []
        self.n_top_results = 3
        print('initializing sweep')
       
    def _print_input_vectors(self):
        print('-'*40)
        print('Input Vectors')
        print(f'  _Vov_1:   {self._Vov_1}')
        print(f'  _Vov_2:   {self._Vov_2}')
        print(f'  _Vov_3:   {self._Vov_3}')
        print(f'  _Vov_N:   {self._Vov_N}')
        print(f'  _Vov_P:   {self._Vov_P}')
        print(f'  _R_LCG:   {self._R_LCG}')
        print(f'  _V1:      {self._V1}')
        print(f'  _ratio_1: {self._ratio_1}')
        print(f'  _ratio_2: {self._ratio_2}')
        print('-'*40)
        
    def _set(self, i1, i2, i3, i4, i5, i6, i7, i8, i9):
        # inputs
        self.Vov_1    = self._Vov_1  [i1] # 6
        self.Vov_2    = self._Vov_2  [i2] # 7
        self.Vov_3    = self._Vov_3  [i3] # 8
        self.Vov_N    = self._Vov_N  [i4] # 5
        self.Vov_P    = self._Vov_P  [i5] # 9
        self.R_LCG    = self._R_LCG  [i6] # 10
        self.V1       = self._V1     [i7] # 11
        self.ratio_1  = self._ratio_1[i8] # 12
        self.ratio_2  = self._ratio_2[i9] # 13
        
        _in = [self.Vov_1,
               self.Vov_2,
               self.Vov_3,
               self.Vov_N,
               self.Vov_P,
               self.R_LCG,
               self.V1,
               self.ratio_1,
               self.ratio_2] 
        
#        print(f'_in: {_in}')
        set_ret = self.tia._set(_in)
        return set_ret
    
    
        
      
    def iterate(self):
        start_ms = int(round(time.time() * 1000))
        count = -1
        val_point_list = []
        self._print_input_vectors()
        for i1 in range(self._Vov_1.shape[0]):
            for i2 in range(self._Vov_2.shape[0]):
                for i3 in range(self._Vov_3.shape[0]):
                    for i4 in range(self._Vov_N.shape[0]):
                        for i5 in range(self._Vov_P.shape[0]):
                            for i6 in range(self._R_LCG.shape[0]):
                                for i7 in range(self._V1.shape[0]):
                                    for i8 in range(self._ratio_1.shape[0]):
                                        for i9 in range(self._ratio_2.shape[0]):
                                            count+=1
                                            set_ret = self._set(i1, i2, i3, i4, i5, i6, i7, i8, i9)
                                            print(f'count: {count}  set_ret: {set_ret}')
                                            if set_ret == 0:
                                                design = self.tia._get_design()
                                                mosfets = self.tia._get_mosfets()
                                                self.results_list.append([design, mosfets])
#                                                print('*'*100)
#                                                out._print_design(design)
#                                                out._print_mosfets(mosfets)   
#                                                print()
#                                                print('-'*100)
                                                
                                                if self.tia.FOM > self._FOM_Vov_1[i1]:
                                                    self._FOM_Vov_1[i1] = self.tia.FOM
                                                if self.tia.FOM > self._FOM_Vov_2[i2]:
                                                    self._FOM_Vov_2[i2] = self.tia.FOM
                                                if self.tia.FOM > self._FOM_Vov_3[i3]:
                                                    self._FOM_Vov_3[i3] = self.tia.FOM   
                                                
                                                if self.tia.FOM > self._FOM_Vov_N[i4]:
                                                    self._FOM_Vov_N[i4] = self.tia.FOM
                                                if self.tia.FOM > self._FOM_Vov_P[i5]:
                                                    self._FOM_Vov_P[i5] = self.tia.FOM
                                                
                                                if self.tia.FOM > self._FOM_R_LCG[i6]:
                                                    self._FOM_R_LCG[i6] = self.tia.FOM
                                                if self.tia.FOM > self._FOM_V1[i7]:
                                                    self._FOM_V1[i7] = self.tia.FOM
                                                if self.tia.FOM > self._FOM_ratio_1[i8]:
                                                    self._FOM_ratio_1[i8] = self.tia.FOM
                                                if self.tia.FOM > self._FOM_ratio_2[i9]:
                                                    self._FOM_ratio_2[i9] = self.tia.FOM    
#                                                print(f'point: {count} valid')
                                                self.valid_point_list.append('p:%d valid' %count)
                                                    
                                            if set_ret == -1:
#                                                print('input out of range')
                                                self.invalid_point_list.append('p:%d oor' %count)
                                            if set_ret == -2:
#                                                print('f3dB not found')
                                                self.invalid_point_list.append('p:%d nf3' %count)

        end_ms = int(round(time.time() * 1000))
        print(f'total time: {end_ms-start_ms} ms')
        print(f'{len(self.valid_point_list)} valid points')
        print(f'{len(self.invalid_point_list)} invalid points')
        
    def _sort(self):
        temp_results_list = self.results_list
        self.results_list = []
        
        _FOM = np.zeros(len(temp_results_list))
        for i in range(len(temp_results_list)):
            _FOM[i] = (temp_results_list[i][0][11])
            
#        print(f'_FOM: {_FOM}')
        print(f'FOM_max: {np.max(_FOM)}')
        argsort_FOM = np.argsort(-_FOM) # negative sign to sort from highes to lowest
#        print(f'args: {argsort_FOM}')
        
        self.top_results_list = []
        for i in range(len(temp_results_list)):
            ind    = argsort_FOM[i]
            result = temp_results_list[ind]
#            print(f'ind: {ind}')
#            print(f'result: {result}')
            self.results_list.append(result)
            if i < self.n_top_results:
                self.top_results_list.append(result)
        
    def _print(self):
        for ind, val in enumerate(self.results_list):
            print('*'*100)
            print(f'design: #{ind+1}')
            out._print_design(val[0])
            out._print_mosfets(val[1])
            
    def _print_top(self):
        for ind, val in enumerate(self.top_results_list):
            print('*'*100)
            print(f'design: #{ind+1}')
            out._print_design(val[0])
            out._print_mosfets(val[1])
        
    def _print_total_iterations(self):
        print(f'total iterations: {self.total_iterations}')
    
    def _plot_single(self, fig_num, name, x, y):
        fig = plt.figure(fig_num)
        fig.clf()
        fig.suptitle(name)
        ax = fig.add_subplot(111)
        ax.set_xlabel(name)
        ax.set_ylabel('FOM')
        ax.plot(x, y)
        ax.grid()

    # Vov_1, Vov_2, Vov_3, Vov_N, Vov_P, R_LCG, V1, ratio_1, ratio_2
    def _plot(self, f1, f2, f3, f4, f5, f6, f7, f8, f9):
        if (f1):
            self._plot_single(1, 'Vov_1', self._Vov_1, self._FOM_Vov_1)
        if (f2):
            self._plot_single(2, 'Vov_2', self._Vov_2, self._FOM_Vov_2)
        if (f3):
            self._plot_single(3, 'Vov_3', self._Vov_3, self._FOM_Vov_3)
        if (f4):
            self._plot_single(4, 'Vov_N', self._Vov_N, self._FOM_Vov_N)
        if (f5):
            self._plot_single(5, 'Vov_P', self._Vov_P, self._FOM_Vov_P)
        if (f6):
            self._plot_single(6, 'R_LCG', self._R_LCG, self._FOM_R_LCG)
        if (f7):
            self._plot_single(7, 'V1',    self._V1, self._FOM_V1)
        if (f8):
            self._plot_single(8, 'ratio_1', self._ratio_1, self._FOM_ratio_1)
        if (f9):
            self._plot_single(9, 'ratio_2', self._ratio_2, self._FOM_ratio_2)
        

def SWEEP_unit_test():
    sweep = SWEEP()
    sweep._print_total_iterations()
    
    sweep.iterate()
#    sweep._print()
    print()
    print()
    sweep._sort()

    print()
    print()
    sweep._print_top()
    sweep._plot(1, 1, 1, 1, 1, 1, 1, 1, 1)
#    sweep._plot()
    
SWEEP_unit_test()











