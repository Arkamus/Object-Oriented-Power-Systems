import matplotlib.pyplot as plt
import CoolProp.CoolProp as CP
import numpy as np
import pandas as pd
import math

class HeatExchanger:
    
    def __init__(self, T_hot_in, T_hot_out, T_cold_in, T_cold_out, m_hot, m_cold, p_hot_in, p_cold_in, hot_fluid, cold_fluid, flow_arrangement, **kwargs):
        
        """
        T_hot_in - temperature of hotter fluid at heat exchanger inlet, st. C
        T_hot_out - temperature of hotter fluid at heat exchanger outlet, st. C
        p_hot_in - pressure of hotter fluid at heat exchanger inlet, kPa(a)
        m_hot - mas flow of hotter fluid, kg/s
        
        T_cold_in - temperature of colder fluid at heat exchanger inlet, st. C
        T_cold_out - temperature of colder fluid at heat exchanger outlet, st. C
        p_cold_in - pressure of colder fluid at heat exchanger inlet, kPa(a)
        m_cold - mas flow of colder fluid, kg/s
        
        flow_arrangement - arangament of fluids flowing throug heat exchanger can be:
            * parallel - fluids are flowing in the same direction ← ← or → →
            * counter - fluids are flowing in counter directions → ←
            * cross - fluids flow directions are crosing each other ↑← or →↑ or ↓← or →↓ or other similar
        
        If any of fluid is a mixture of gases such as:
            * air
            * flue gases from combustion systems
            * exchaust gases from other proceses like melting metals
            
        The name schould be written as:
            'HEOS::CO2[z_CO2]&O2[z_O2]&N2[z_N2]&H2O[z_H2O]'
        
        Where:
            CO2, O2, N2, H2O and so on - name of mixture element
            z_CO2, z_O2, z_N2, z_H2O and so on - molar or volumetric content of mixture element as value form 0.0-1.0 with sum of each component not exceeding 1.0
            
        If contents needs to be variables then write name of fluid as:
            f'HEOS::CO2[{z_CO2}]&O2[{z_O2}]&N2[{z_N2}]&H2O[{z_H2O}]'
        
        Sources:
            https://en.wikipedia.org/wiki/Heat_exchanger
            
        !!!WARNING!!!
        No condensation in mixtures is allowed at this moment.
        
        """
        self.T_hot_in = T_hot_in
        self.T_hot_out = T_hot_out
        self.p_hot_in = p_hot_in
        self.m_hot = m_hot
        
        self.T_cold_in = T_cold_in
        self.T_cold_out = T_cold_out
        self.p_cold_in = p_cold_in
        self.m_cold = m_cold
        
        self.hot_fluid = hot_fluid
        self.cold_fluid = cold_fluid
        self.flow_arrangement = flow_arrangement
        
        if flow_arrangement == 'parallel':
            pass
        elif flow_arrangement == 'counter':
            pass
        elif flow_arrangement == 'cross':
            
            raise Exception("Flow arrangement can be only: parallel, counter or cross")
            
        else:
            
            raise Exception("Flow arrangement can be only: parallel, counter or cross")
            
        """
        Checking if one of fluids is a user defined mixture. If true then creates function T=f(Hmass). 
        This is due to limitations of CoolProp to have only limited avaliable input pairs parametes for user defined mixtures:
            * P/Q (pressure/quality)
            * T/Q (temperature/quality)
            * T/P (temperature/pressure)
        
        Sources:
            http://www.coolprop.org/fluid_properties/Mixtures.html#phase-envelope            
        """
         
        if 'HEOS' in self.hot_fluid:
            
            'Number of points for aproximation'
            num = 50
            
            p_data = np.full((1,num), self.p_hot_in).flatten()
            T_data = np.linspace(start = self.T_hot_in, stop = self.T_cold_out,num = num, endpoint=True)
             
            h_data = CP.PropsSI('Hmass','T|gas',T_data + 273.15,'P', p_data * 1000, self.hot_fluid)/1000
            
            self.T_hot_aprox = np.poly1d(np.polyfit(x = h_data, y = T_data , deg = 6))
            
    def LMTD(self):
        
        """
        Calculates Logarithmic mean temperature difference
        
        dTA - temperature difference between two streams (hot and cold) at one end of heat exchanger (A)
        dTB - temperature difference between two streams (hot and cold) at other end of heat exchanger (B)
        
        Source:
           https://en.wikipedia.org/wiki/Logarithmic_mean_temperature_difference 
        """
        if self.flow_arrangement == 'parallel':
        
            dTA = self.T_hot_in-self.T_cold_in
            dTB = self.T_hot_out-self.T_cold_out
        
        elif self.flow_arrangement == 'counter':
            
            dTA = self.T_hot_in-self.T_cold_out
            dTB = self.T_hot_out-self.T_cold_in
            
        elif self.flow_arrangement == 'cross':
            
            raise Exception("Cross flow arrangement not yet implemented, can be only: parallel or counter")
        
        self.LMTD = (dTA-dTB)/math.log(dTA/dTB)
        
    def HTA_LMTD(self, alfa):
        
        """
        Calculates heat transfer area using LMTD.
        
        Inputs:
            alafa - heat transfer coefficien
        """
            
    def calc_temp_dist(self, points = 4):
        
        """
        Preliminary ploting temperature distribution in function of transfered heat stream.
        Simplifyng assumption: no pressure loses in heat exchanger.
        
        ponits - number of characteristic points on graph 
        """
        self.points = points
        
        p_hot_out = self.p_hot_in
        p_cold_out = self.p_cold_in
        
        'Checking if hot fluid is a user defined mixture'
        if 'HEOS' in self.hot_fluid:
            T_param = 'T|gas'
            #Dew_point_hot_fluid = CP.PropsSI('T', 'Q', 1, 'P', self.p_hot_in, self.hot_fluid)-273.15
        else:
            T_param = 'T'
        
        'Calculating properties in key point'
        self.h_hot_in = CP.PropsSI('Hmass', T_param, self.T_hot_in+273.15, 'P', self.p_hot_in*1000, self.hot_fluid)/1000
        self.h_hot_out = CP.PropsSI('Hmass', T_param, self.T_hot_out+273.15, 'P', p_hot_out*1000, self.hot_fluid)/1000
        
        self.Q_hot = self.m_hot*(self.h_hot_in-self.h_hot_out)
        
        self.h_cold_in = CP.PropsSI('Hmass', 'T', self.T_cold_in+273.15, 'P', self.p_cold_in*1000, self.cold_fluid)/1000
        self.h_cold_out = CP.PropsSI('Hmass', 'T', self.T_cold_out+273.15, 'P', p_cold_out*1000, self.cold_fluid)/1000
        
        self.Q_cold = self.m_cold*(self.h_cold_out-self.h_cold_in)
        
        self.Q_step = self.Q_hot/(points-1)
        
        'Creating an array of intial steps to be evaluated and later used to make a plot'
        self.Q_intervals = np.linspace(start = 0,stop = self.Q_hot, num = self.points, endpoint = True)
        
        'Creating data frame for heat exchanger temperature distribution'
        self.temp_dist_df = pd.DataFrame(columns = ['Q', 'T_hot', 'T_cold', 'dT', 'h_hot', 'h_cold'])
        
        'Phase change is not allowed for incompressible fluids. Checking if hot fluid is incompressible'
        if 'INCOMP' in self.hot_fluid:
            
            phase_hot_in  = 'liquid'
            phase_hot_out = 'liquid'
            
        else:
         
            'Checking if ther is a phase change of hot fluid'
            phase_hot_in  = CP.PhaseSI(T_param, self.T_hot_in+273.15, 'P', self.p_hot_in*1000, self.hot_fluid)
            phase_hot_out = CP.PhaseSI(T_param, self.T_hot_out+273.15, 'P', self.p_hot_in*1000, self.hot_fluid)
            
        'Phase change is not allowed for incompressible fluids. Checking if cold fluid is incompressible'
        if 'INCOMP' in self.cold_fluid:
            
            phase_cold_in  = 'liquid'
            phase_cold_out = 'liquid'
            
        else:
            
            'Checking if there is a phase change of cold fluid'
            phase_cold_in  = CP.PhaseSI(T_param, self.T_cold_in+273.15, 'P', self.p_cold_in*1000, self.cold_fluid)
            phase_cold_out = CP.PhaseSI(T_param, self.T_cold_out+273.15, 'P', self.p_cold_in*1000, self.cold_fluid)
        
        'Calculating stauration points of cold fluid if it is changid phase in heat exchanger'
        if phase_cold_in != phase_cold_out:
            
            h_cold_evap0 = CP.PropsSI('Hmass', 'Q', 0, 'P', self.p_cold_in*1000, self.cold_fluid)/1000
            h_cold_evap1 = CP.PropsSI('Hmass', 'Q', 1, 'P', self.p_cold_in*1000, self.cold_fluid)/1000
            
            T_cold_evap = CP.PropsSI('T', 'Q', 1, 'P', self.p_cold_in*1000, self.cold_fluid)-273.15
            
            Q_cold_superheat = self.m_cold*(self.h_cold_out-h_cold_evap1)
            Q_cold_evaporation =  self.m_cold*(self.h_cold_out-h_cold_evap0)
        
            self.Q_intervals = np.append(self.Q_intervals, [Q_cold_superheat, Q_cold_evaporation])
            self.Q_intervals = np.sort(self.Q_intervals)
        
        'Calculating stauration points of hot fluid if it is changid phase in heat exchanger'
        if phase_hot_in != phase_hot_out:
            
            h_hot_evap0 = CP.PropsSI('Hmass', 'Q', 0, 'P', self.p_hot_in*1000, self.hot_fluid)/1000
            h_hot_evap1 = CP.PropsSI('Hmass', 'Q', 1, 'P', self.p_hot_in*1000, self.hot_fluid)/1000
            
            T_hot_evap = CP.PropsSI('T', 'Q', 1, 'P', self.p_hot_in*1000, self.hot_fluid)-273.15
            
            Q_hot_superheat = self.m_hot*(self.h_hot_in-h_hot_evap1)
            Q_hot_evaporation =  self.m_hot*(self.h_hot_in-h_hot_evap0)
        
            self.Q_intervals = np.append(self.Q_intervals, [Q_hot_superheat, Q_hot_evaporation])
            self.Q_intervals = np.sort(self.Q_intervals)
        
        'Calculating temperatures for each heat stream in frame'
        self.temp_dist_df['Q'] = self.Q_intervals
        
        counter = 0
        
        for Q in self.Q_intervals:
            
            h_hot = self.h_hot_in - (Q/self.m_hot)
            
            if 'HEOS' in self.hot_fluid:
                
                T_hot = self.T_hot_aprox(h_hot)
                
            else:
                
                T_hot = CP.PropsSI('T', 'Hmass', h_hot*1000, 'P', self.p_hot_in*1000, self.hot_fluid)-273.15
        
            h_cold = self.h_cold_out - Q/self.m_cold
            T_cold = CP.PropsSI('T', 'Hmass', h_cold*1000, 'P', self.p_cold_in*1000, self.cold_fluid)-273.15
            
            self.temp_dist_df.loc[counter, 'T_hot'] = T_hot
            self.temp_dist_df.loc[counter, 'T_cold'] = T_cold
            
            self.temp_dist_df.loc[counter, 'h_hot'] = h_hot
            self.temp_dist_df.loc[counter, 'h_cold'] = h_cold
            
            counter += 1
        
        'Realigning data acordning to flow arrangement'
        if self.flow_arrangement == 'counter': 
            
            self.temp_dist_df['T_cold'] = self.temp_dist_df['T_cold'].sort_values(ascending=False).values
            self.temp_dist_df['dT'] =  self.temp_dist_df['T_hot'] - self.temp_dist_df['T_cold']
            
        elif self.flow_arrangemnt == 'parallel':
            
            self.temp_dist_df['dT'] =  self.temp_dist_df['T_hot'] - self.temp_dist_df['T_cold']
            
        elif self.flow_arrangement == 'cross':
            pass
        
        self.dT_min = min(self.temp_dist_df['dT'].values)
        
    
    def plot_temp_dist_df(self):
        
        fig, axs = plt.subplots(nrows = 1, ncols = 1, constrained_layout = True )
        
        ((ax)) = axs
        
        QX = self.temp_dist_df.Q.values
        
        ax.plot(QX, self.temp_dist_df.T_hot.values, label = f'Hotter fluid: {self.hot_fluid}', c = 'red')
        ax.plot(QX, self.temp_dist_df.T_cold.values,label = f'Colder fluid: {self.cold_fluid}', c = 'blue' )
        ax.grid(axis ='both', linestyle = ':')
        ax.set_xlim(min(QX), max(QX))
        ax.set_xlabel('Heat transfered, kWt')
        ax.set_ylabel('Temperaure, $^\circ$C')
        ax.legend()
        
        return fig, ax
