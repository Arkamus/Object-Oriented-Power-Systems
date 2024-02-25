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
            
    def log_mean_temp(self, T_hot_in, T_hot_out, T_cold_in, T_cold_out, flow_arrangement):
        
        """
        Calculates logarithmic mean temperature difference method - source 1. 
        If there is phase change in heat exchanger, then LMTD should be calculated fore each zone:
            * preheating/ precooling
            * evaporation/condensation
            * superheating/subcoling
        
        
        Sources:
           1. https://en.wikipedia.org/wiki/Logarithmic_mean_temperature_difference 
        """
        if flow_arrangement == 'parallel':
        
            dTA = T_hot_in-T_cold_in
            dTB = T_hot_out-T_cold_out
        
        elif flow_arrangement == 'counter':
            
            dTA = T_hot_in-T_cold_out
            dTB = T_hot_out-T_cold_in
            
        elif flow_arrangement == 'cross':
            
            raise Exception("Cross flow arrangement not yet implemented, can be only: parallel or counter")
        
        LMTD = (dTA-dTB)/math.log(dTA/dTB)
        
        return LMTD
    
    def HTA_LMTD(self):
        
        """
        Calculation of heat transfer area using LMTD method
        U - overall heat transfer coefficient, W/m^2*K
        
        Sources:
        1. https://www.politesi.polimi.it/retrieve/a81cb05a-c4e7-616b-e053-1605fe0a889a/Marco%20Astolfi%20-%20PhD-last.pdf
        
        """
        
        'Checking if there is phase change in heat exchanger - if it is then heat exchanger will be calculated as three separated heat exchanger - for three zones of heat exchanger'
        if self.hot_fluid_phase_change == True:
            
            if 'HEOS' not in self.hot_fluid and 'INCOMP' in self.cold_fluid:
                
                'For phase change temp_dist_diff must be splited to three different table for each section heat exchanger.'
                'Olny the first and last row are significant for preliminary design.'
                
                coller = self.temp_dist_df[(self.temp_dist_df.Q_hot >=1.0)].iloc[[0, -1]].reset_index(drop=True)
                condenser = self.temp_dist_df[(self.temp_dist_df.Q_hot.between(0.0,1.0))].iloc[[0, -1]].reset_index(drop=True)
                subcoller = self.temp_dist_df[(self.temp_dist_df.Q_hot <= min(condenser.Q_hot))].iloc[[0, -1]].reset_index(drop=True)
                
                'Calculation of LMTD for each section'
                self.LMTD_COL   = self.log_mean_temp(T_hot_in=coller.T_hot[0], T_hot_out=coller.T_hot[1], T_cold_in =coller.T_cold[1], T_cold_out = coller.T_cold[0], flow_arrangement=self.flow_arrangement)
                self.LMTD_COND  = self.log_mean_temp(T_hot_in=condenser.T_hot[0], T_hot_out=condenser.T_hot[1], T_cold_in =condenser.T_cold[1], T_cold_out = condenser.T_cold[0], flow_arrangement=self.flow_arrangement)
                self.LMTD_SUB   = self.log_mean_temp(T_hot_in=subcoller.T_hot[0], T_hot_out=subcoller.T_hot[1], T_cold_in =subcoller.T_cold[1], T_cold_out = subcoller.T_cold[0], flow_arrangement=self.flow_arrangement)
                
                'Calcualtion of heat transfered in each section'
                self.Q_COL = coller.Q[1]-coller.Q[0]
                self.Q_COND = condenser.Q[1]-condenser.Q[0]
                self.Q_SUB = subcoller.Q[1]-subcoller.Q[0]
                
                'Calculation of heat transfer area for each section and heat exhcanger as a whole.'
                self.A_COL  = (self.Q_COL*1000)/(343.01*self.LMTD_COL)
                self.A_COND = (self.Q_COND*1000)/(623.47*self.LMTD_COND)
                self.A_SUB  = (self.Q_SUB*1000)/(623.47*self.LMTD_SUB)
                self.A = self.A_COL+self.A_COND+self.A_SUB
            
        elif self.cold_fluid_phase_change == True:
            
            if 'HEOS' in self.hot_fluid and 'INCOMP' not in self.cold_fluid:
                
                
                'For phase change temp_dist_diff must be splited to three different table for each section heat exchanger.'
                'Olny the first and last row are significant for preliminary design.'
                
                table = self.temp_dist_df 
                
                superheater = self.temp_dist_df[(self.temp_dist_df.Q_cold >=1.0)].iloc[[0, -1]].reset_index(drop=True)
                evaporator  = self.temp_dist_df[(self.temp_dist_df.Q_cold.between(0.0,1.0))].iloc[[0, -1]].reset_index(drop=True)
                preheater   = self.temp_dist_df[(self.temp_dist_df.Q_cold <= min(evaporator.Q_cold))].iloc[[0, -1]].reset_index(drop=True)
                
                'Calculation of LMTD for each section'
                self.LMTD_SUP   = self.log_mean_temp(T_hot_in=superheater.T_hot[0], T_hot_out=superheater.T_hot[1], T_cold_in =superheater.T_cold[1], T_cold_out = superheater.T_cold[0], flow_arrangement=self.flow_arrangement)
                self.LMTD_EVAP  = self.log_mean_temp(T_hot_in=evaporator.T_hot[0], T_hot_out=evaporator.T_hot[1], T_cold_in =evaporator.T_cold[1], T_cold_out = evaporator.T_cold[0], flow_arrangement=self.flow_arrangement)
                self.LMTD_PRE   = self.log_mean_temp(T_hot_in=preheater.T_hot[0], T_hot_out=preheater.T_hot[1], T_cold_in =preheater.T_cold[1], T_cold_out = preheater.T_cold[0], flow_arrangement=self.flow_arrangement)
                
                'Calcualtion of heat transfered in each section'
                self.Q_SUP = superheater.Q[1]-superheater.Q[0]
                self.Q_EVAP = evaporator.Q[1]-evaporator.Q[0]
                self.Q_PRE = preheater.Q[1]-preheater.Q[0]
                
                'Calculation of heat transfer area for each section and heat exhcanger as a whole.'
                self.A_SUP  = (self.Q_SUP*1000)/(343.01*self.LMTD_SUP)
                self.A_EVAP = (self.Q_EVAP*1000)/(623.47*self.LMTD_EVAP)
                self.A_PRE  = (self.Q_PRE*1000)/(623.47*self.LMTD_PRE)
                self.A = self.A_SUP+self.A_EVAP+self.A_PRE
                
        else:
            self.LMTD = self.log_mean_temp(self.T_hot_in, self.T_hot_out, self.T_cold_in, self.T_cold_out, self.flow_arrangement)
            
            if  'HEOS' or 'Air' in self.hot_fluid and 'INCOMP' in self.cold_fluid:
                self.U = 465.40 # W/m^2*K
                self.A = (self.Q_cold*1000)/(self.U*self.LMTD)
        
    def calc_temp_dist(self, points = 4):
        
        """
        Preliminary ploting temperature distribution in function of transfered heat stream.
        Simplifyng assumption: no pressure loses in heat exchanger.
        
        ponits - number of characteristic points on graph 
        """
        
        def check_phase(T,p, fluid, Quality = None):
            
            """
            Wrapper function to check phase of fluid in simple and more complex cases e.g. incompressible fluids, phase change (evaporation, condensation), saturation
            
            T - temperature, C.deg
            p - pressure, kPa(a)
            fluid - fluid name
            
            INCOMP - incompressible fluid, usually liquid
            HEOS - user defined fluid - mixture e.g. exhaust gas, air or other gaseous mixture, usually gas
            """
        
            if 'INCOMP' in fluid:
                
                return 'liquid'
            
            elif 'HEOS' in fluid:
            
                return CP.PhaseSI('T|gas', T+273.15, 'P', p*1000, fluid)
            
            else:
                
                phase = CP.PhaseSI('T', T+273.15, 'P', p*1000, fluid)
                
                if 'Saturation pressure' in str(phase):
                    
                    phase = 'twophase'
        
                return phase
        
        def vapour_quality(h, p, fluid, T = None):
            
            """
            Wrapper function to check vapour quality Q
            
            Q - vapour quality 0-1, 0 saturated liquid, 1 - dry vapour, -1 - subcoled liquid, 2 - superheated gas
            h - enthalpy, kJ/kg*K
            p - pressure, kPa(a)
            fluid - fluid name
            
            INCOMP - incompressible fluid, usually liquid
            HESO   - user defined fluid, usually gas
            """
        
            if 'INCOMP' in fluid:
                
                return -1.0
            
            elif 'HEOS' in fluid:
                
                return 1
            
            else:
                
                Quality = CP.PropsSI('Q', 'Hmass', h*1000, 'P', p*1000, fluid)
                
                if Quality < 0:
                    
                    phase = CP.PhaseSI('T', T+273.15, 'P', p*1000, fluid)
                    
                    if phase == 'gas':
                        
                        Quality = 2.0
                    
                    elif phase == 'liquid':
                    
                        Quality = -1.0
                        
                return Quality
                    
        
        self.points = points
        
        p_hot_out = self.p_hot_in
        p_cold_out = self.p_cold_in
        
        'Checking if hot fluid is a user defined mixture'
        if 'HEOS' in self.hot_fluid:
            T_param = 'T|gas'
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
        
        'Checking phase of hot fluid'
        phase_hot_in = check_phase(self.T_hot_in,self.p_hot_in, self.hot_fluid)
        phase_hot_out = check_phase(self.T_hot_out,p_hot_out, self.hot_fluid)
        
        'Checking phase of hot fluid'
        phase_cold_in = check_phase(self.T_cold_in,self.p_cold_in, self.cold_fluid)
        phase_cold_out = check_phase(self.T_cold_out,p_cold_out, self.cold_fluid)
        
        'Calculating stauration points of cold fluid if it is changid phase in heat exchanger'
        if phase_cold_in != phase_cold_out:
            
            self.cold_fluid_phase_change = True
            
            h_cold_evap0 = CP.PropsSI('Hmass', 'Q', 0, 'P', self.p_cold_in*1000, self.cold_fluid)/1000
            h_cold_evap1 = CP.PropsSI('Hmass', 'Q', 1, 'P', self.p_cold_in*1000, self.cold_fluid)/1000
            
            self.T_cold_evap = CP.PropsSI('T', 'Q', 1, 'P', self.p_cold_in*1000, self.cold_fluid)-273.15
            
            self.Q_cold_superheat = self.m_cold*(self.h_cold_out-h_cold_evap1)
            self.Q_cold_evaporation =  self.m_cold*(h_cold_evap1-h_cold_evap0)
            self.Q_cold_preheat = self.m_cold*(h_cold_evap0-self.h_cold_in)
        
            self.Q_intervals = np.append(self.Q_intervals, [self.Q_cold_superheat, self.Q_cold_superheat+self.Q_cold_evaporation])
            self.Q_intervals = np.sort(self.Q_intervals)
        
        else:
            
            self.cold_fluid_phase_change = False
        
        'Calculating stauration points of hot fluid if it is changid phase in heat exchanger'
        if phase_hot_in != phase_hot_out:
            
            self.hot_fluid_phase_change = True
            
            h_hot_evap0 = CP.PropsSI('Hmass', 'Q', 0, 'P', self.p_hot_in*1000, self.hot_fluid)/1000
            h_hot_evap1 = CP.PropsSI('Hmass', 'Q', 1, 'P', self.p_hot_in*1000, self.hot_fluid)/1000
            
            self.T_hot_evap = CP.PropsSI('T', 'Q', 1, 'P', self.p_hot_in*1000, self.hot_fluid)-273.15
            
            self.Q_hot_cool = self.m_hot*(self.h_hot_in-h_hot_evap1)
            self.Q_hot_condesation =  self.m_hot*(h_hot_evap1-h_hot_evap0)
            self.Q_hot_subcooling = self.m_hot*(h_hot_evap0-self.h_hot_out)
        
            self.Q_intervals = np.append(self.Q_intervals, [self.Q_hot_cool, self.Q_hot_cool+self.Q_hot_condesation])
            self.Q_intervals = np.sort(self.Q_intervals)
        
        else: 
            
            self.hot_fluid_phase_change = False
        
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
            
            self.temp_dist_df.loc[counter, 'p_hot'] = self.p_hot_in
            self.temp_dist_df.loc[counter, 'p_cold'] = self.p_cold_in
            
            self.temp_dist_df.loc[counter, 'h_hot'] = h_hot
            self.temp_dist_df.loc[counter, 'h_cold'] = h_cold
            
            self.temp_dist_df.loc[counter, 'hot_fluid'] = self.hot_fluid
            self.temp_dist_df.loc[counter, 'cold_fluid'] = self.cold_fluid
            
            self.temp_dist_df.loc[counter, 'Q_hot'] = vapour_quality(h_hot,  self.p_hot_in, self.hot_fluid, T_hot)#CP.PropsSI('Q', 'Hmass', h_hot*1000, 'P', self.p_hot_in*1000, self.hot_fluid)
            self.temp_dist_df.loc[counter, 'Q_cold'] = vapour_quality(h_cold,  self.p_cold_in, self.cold_fluid, T_cold) #CP.PropsSI('Q', 'Hmass', h_cold*1000, 'P', self.p_cold_in*1000, self.cold_fluid)
            
            self.temp_dist_df.loc[counter, 'phase_hot_fluid'] =  check_phase(T_hot,self.p_hot_in, self.hot_fluid)
            self.temp_dist_df.loc[counter, 'phase_cold_fluid'] = check_phase(T_cold,self.p_cold_in, self.cold_fluid)
            
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
        
        def HTA_LMTD(self):
            
            """
            Calculation of heat transfer area (A) using LMTD method.
            
            U - overal heat transfer coefficient, W/m^2*K
            A - heat transfer area, m^2
            
            T_hot_evap - temperature of evaportation (phase change) of hotter fluid, deg.C or K
            T_cold_ecvap - temperature of evaporation (phase change) of colder fluid, deg.C or K
            
            Q_cool - heat flux during cooling of hot fluid, kW
            Q_condensation - heat flux during condesnation of hot fluid, kW
            Q_subcooling - heat flux during subcooling of hot fluid, kW
            
            Q_preheat - heat flux during preheat of cold fluid, kW
            Q_evaporation - heat flux during evaporation of cold fluid, kW
            Q_superheat - heat flux during superheating of cold fluid, kW
            
            Q_hot - heat flux transfered in heat exchanger from hot fluid to cold fluid (Q_hot=Q_cold), kW
            """
            
            if self.cold_fluid_phase_change == True or self.hot_fluid_phase_change == True:
                pass
            else:
                pass
    
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
