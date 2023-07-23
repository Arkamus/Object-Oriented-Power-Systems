# -*- coding: utf-8 -*-
"""
Created on Wed Apr  5 11:51:24 2023

@author: Arek
"""
import math
import CoolProp.CoolProp as CP
import matplotlib.pyplot as plt
import numpy as np

class TurbineDesign:
    
    """
    Class to handle turbine calculations. It is one of elements that power systems are build like heat exchangers and pumps.
    
    Obligatory inputs:
    T_in - temperature of vapour at turbine inlet, C
    p_in - pressure of vapour at turbine inlet, kPa(a)
    p_out - pressure of vapour at turbine outlet, kPa(a)
    mf - mass flow of working fluid, kg/s
    medium - working fluid ex. Water, Toluene, R1233zd(E), etc.
    """
    
    def __init__(self, T_in, p_in, p_out, mf,medium):
        self.T_in = T_in
        self.p_in = p_in
        self.p_out = p_out
        self.mf = mf
        self.medium = medium
        
        """Calculation of basic thermodynamic properties for use in other calculations in class.
        
        Sources:
        * https://en.wikipedia.org/wiki/List_of_thermodynamic_properties
        * https://en.wikipedia.org/wiki/State_function
        
        Outputs:
        h_in - specific enthalpy at turbine inlet, kJ/kg
        s_in - specific entropy at turbine inlet, kJ/(kg*K) or kJ/(kg*C) where K is temperature in K-Kelvins or C - Celciuss degres
        
        h_out_s - specific enthalpy at turbine outlet after isentropic expansion (eta_iT=1.0), kJ/kg
        s_out_s - specific entropy at turbine outlet after isentropic expansion (eta_iT=1.0), kJ/(kg*K) or kJ/(kg*C) where K is temperature in K-Kelvins or C - Celciuss degres
        T_out_s - temperature at turbine outlet after isentropic expansion (eta_iT=1.0), C
        """
        
        'Parameters at turbine inlet'
        self.h_in = CP.PropsSI('Hmass', 'T', self.T_in+273.15, 'P', self.p_in*1000, self.medium)/1000
        self.s_in = CP.PropsSI('Smass', 'T', self.T_in+273.15, 'P', self.p_in*1000, self.medium)/1000
        
        'Parameters at turbine outlet after isentropic expansion (eta_iT=1.0), kJ/(kg*K) or kJ/(kg*C) where K is temperature in K-Kelvins or C - Celciuss degres'
        self.s_out_s = self.s_in 
        self.h_out_s = CP.PropsSI('Hmass', 'Smass', self.s_out_s*1000, 'P', self.p_out*1000, self.medium)/1000
        self.T_out_s = CP.PropsSI('T', 'Smass', self.s_out_s*1000, 'P', self.p_out*1000, self.medium)-273.15
        
    def efficiency_estimation(self, n_stages):
        
        """
        Simplified, preliminary estimation of design turbine internal efficiency. Change in inlet and outlet parameters affects value of efficiency.
        
        Inputs:
        n_stages - number of turbine stages (1-3)
        
        Outputs:
        rho_turbine_in - density of vapour at turbine inlet, kg/m^3
        V_turbine_in - vapour volumetric flow at turbine inlet, m^3/s
        
        rho_turbine_out - density of vapour at turbine outlet, kg/m^3
        rho_turbine_out_s - density of vapour at outlet after isentropic expansion (eta_iT=1.0), kg/m^3, kg/m^3
        V_turbine_out_s - vapour volumetric flow at turbine outlet after isentropic expansion (eta_iT=1.0), m^3/s
        
        delta_h_s - enthalpy drop in turbine after isentropic expansion (eta_iT=1.0), kJ/kg
        eta_iT - estimated internal efficiency of turbine 0.0-1.0, - 
        SP - dimensionless size parameter, -
        V_r - dimensionless volumetric flow ration, -
        
        Sources:
        * Macchi E, Astolfi M. Organic rankine cycle (ORC) power systems: Technologies and
          applications. Duxford, United Kingdom: Woodhead Publishing is an imprint of Elsevier;
          2017.
          
        """
        self.n_stages = n_stages 
        
        self.rho_turbine_in = CP.PropsSI('Dmass', 'T', self.T_in+273.15, 'P', self.p_in*1000, self.medium)
        self.rho_turbine_out_s = CP.PropsSI('Dmass', 'Smass', self.s_out_s*1000, 'P', self.p_out*1000, self.medium)
        self.delta_h_s = (self.h_in-self.h_out_s)*1000
        
        self.V_turbine_in = self.mf/self.rho_turbine_in
        self.V_turbine_out_s =self.mf/self.rho_turbine_out_s
        
        self.SP = (self.V_turbine_out_s**(1/2))/(self.delta_h_s**(1/4))

        self.V_r = self.V_turbine_out_s/self.V_turbine_in
                
        if n_stages == 1:                                                      
            Ai = np.array([0.90831500, -0.05248690, -0.04799080, -0.01710380, -0.00244002, 0, 0.04961780, -0.04894860, 0.01171650, -0.00100473, 0.05645970, -0.01859440, 0.01288860, 0.00178187, -0.00021196, 0.00078667])
        elif n_stages == 2:
            Ai = [0.923406, -0.0221021, -0.0233814, -0.00844961, -0.0012978, -0.00069293, 0.0146911, -0.0102795, 0, 0.000317241, 0.0163959, -0.00515265, 0.00358361, 0.000554726, 0, 0.000293607]
        elif n_stages == 3:
            Ai = [0.932274, -0.01243, -0.018, -0.00716, -0.00118, -0.00044, 0, 0, -0.0016, 0.000298, 0.005959, -0.00163, 0.001946, 0.000163, 0, 0.000211]
        else:
            raise Exception("Number of turbine stages (n_stages) in this method should be 1,2 or 3.")
        
        self.eta_iT = (1 * Ai[0]) + \
                      (math.log(self.SP)*(Ai[1])) + \
                      (((math.log(self.SP))**2)*(Ai[2])) + \
                      (((math.log(self.SP))**3)*(Ai[3])) + \
                      (((math.log(self.SP))**4)*(Ai[4])) + \
                      (self.V_r*Ai[5]) + \
                      (math.log(self.V_r)*(Ai[6])) + \
                      ((math.log(self.V_r)**2)*(Ai[7])) + \
                      ((math.log(self.V_r)**3)*(Ai[8])) + \
                      ((math.log(self.V_r)**4)*(Ai[9])) + \
                      (math.log(self.V_r)*math.log(self.SP)*Ai[10]) + \
                      ((math.log(self.V_r)**2)*math.log(self.SP)*(Ai[11])) + \
                      (math.log(self.V_r)*(math.log(self.SP)**2)*(Ai[12])) + \
                      ((math.log(self.V_r)**3)*math.log(self.SP)*(Ai[13])) + \
                      ((math.log(self.V_r)**3)*(math.log(self.SP)**2)*(Ai[14])) + \
                      ((math.log(self.V_r)**2)*(math.log(self.SP)**3)*(Ai[15]))        
                                                                
        self.power(self.eta_iT)
                                                                 
    def power(self, eta_iT):
        
        """
        Calculates internal power output of turbine according to given internal efficiency.
        In case of using iterative calculations, the change of inlet and outlet parameters will not affect the interal turbine efficiency.
        
        Inputs:
        T_in - vapour temperature at turbine inelt, C
        p_in - vapour pressure at turbine inlet, kPa(a)
        p_out - presure at turbine outlet, kPa(a)
        mf - medium (vapour) flow through turbine, kg/s
        medium - type of medium flowing through ex. 'Toluene', 'Water'
        eta_iT - fixed internal efficeincy of turbine 0.0-1.0, -
        
        Main outputs:
        N_iT - internal power generated by turbine, kW
        T_out - temperature of medium at turbine ouitlet, C
        
        Other outputs:
        h_in - enthalpy of vapour at turbine inlet, kJ/kg
        s_in - entropy of vapour at turbine inlet, kJ/(kg*C) or kJ/(kg*K)
        h_out - enthalpy of vapour at turbine outlet, kJ/kg
        h_out_s - enthalpy of vapour at tyrbine outlet after ideal (eta=1.0) expansion in turbine, kJ/kg
        s_out - entropy of vapur at turbine outlet
        s_out_s - entropy of vapour at turbine outlet after ideal (eta = 1.0) expansion in turbine, kJ/(kg*C) or kJ/(kg*K)
        """
        self.eta_iT = eta_iT
        self.h_out = self.h_in-self.eta_iT*(self.h_in-self.h_out_s)
        self.T_out = CP.PropsSI('T', 'Hmass', self.h_out*1000, 'P', self.p_out*1000, self.medium)-273.15
        self.s_out = CP.PropsSI('Smass', 'T', self.T_out*1000, 'P', self.p_out*1000, self.medium)/1000
        self.N_iT = self.mf*(self.h_in-self.h_out)
        
    def Ts_diagram(self):
        
        """
        Plots Temperature - specific entropy (T-s) diagram of vapour expansion in turbine.
        !!!WARNING!!!
        Drawing can be done only if previous calculations with fixed or estimated efficiency were done.
        """
        
        fig, ax = plt.subplots(nrows = 1, ncols = 1, dpi = 400)
        
        'Diagram description'
        ax.set_ylabel('Temperature, $^\circ$C')                                 # y-axis label
        ax.set_xlabel('Specific entropy, ' + r'$\frac{kJ}{kg \cdot^\circ C}$')  # x-axis label
        
        'Turbine inlet'
        ax.scatter(self.s_in, self.T_in, marker = 'x', c = 'red')
        ax.text(self.s_in, self.T_in, '$T_{in}$')
        
        'Turbine outlet - isentropic, after expansion in trubine with eta_iT = 1.0'
        ax.scatter(self.s_out_s, self.T_out_s, marker = 'x', c = 'red')
        ax.text(self.s_out_s, self.T_out_s, '$T_{out,s}$')
        
        'Turbine outlet'
        ax.scatter(self.s_out, self.T_out, marker = 'x', c = 'red')
        ax.text(self.s_out, self.T_out, '$T_{out}$')
        
        'Isentropic expansion line'
        ax.plot([self.s_in, self.s_out_s],[self.T_in, self.T_out_s], label = 'Isentropic expansion line', linestyle = ':')
        
        'Estimated expansion line'
        ax.plot([self.s_in, self.s_out],[self.T_in, self.T_out], label = 'Estimated expansion line', linestyle = ':')
        
        'General diagram settings'
        yticks = ax.get_yticks()
        xticks = ax.get_xticks()
        
        max_y_tick = max(yticks)
        max_x_tick = max(xticks)
        
        ax.set_ylim(bottom = 0.0, top = max_y_tick)
        ax.set_xlim(left = 0.0, right = max_x_tick)
        
        yticks = ax.get_yticks()
        xticks = ax.get_xticks()
        
        max_y_tick = max(yticks)
        max_x_tick = max(xticks)
        
        min_y_tick = min(yticks)
        min_x_tick = min(xticks)
        
        len_yticks = len(yticks)
        len_xticks = len(xticks)
        
        y_tick_step = (max_y_tick-min_y_tick)/(len_yticks-1)
        x_tick_step = (max_x_tick-min_x_tick)/(len_xticks-1)
        
        ax.set_ylim(bottom = 0.0, top = max_y_tick+y_tick_step)
        ax.set_xlim(left = 0.0, right = max_x_tick+x_tick_step)
        
        'Ploting inlet presure line'
        p_out_line = np.full(len_xticks, self.p_out)
        s_out_line = xticks
        T_out_line = CP.PropsSI('T', 'P',  p_out_line*1000, 'S', s_out_line*1000, self.medium)-273.15
        
        'Legend'
        ax.legend(title = f'medium = {self.medium}\n' + f'number of stages = {self.n_stages} \n'+'$\eta_{iT}$='+f'{round(self.eta_iT,4)}', loc='lower right')
        
        'Title'
        ax.set_title(label = 'Turbine T-s diagram')
        
        'Adding grid'
        ax.grid(True, color = "grey", linestyle = "--", axis = 'both')
        
        plt.tight_layout()
    
    def picture(self):
        """
        Simple visualization of turbine.
        """
        
        pass

test_turbine1 = TurbineDesign(155, 3620.0, 1568.5, 237.71, 'R125')
test_turbine1.efficiency_estimation(n_stages= 1)
test_eta_iT1 = test_turbine1.eta_iT
test_turbine1.Ts_diagram()
        