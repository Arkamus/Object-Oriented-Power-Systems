# -*- coding: utf-8 -*-
"""
Created on Wed Apr  5 11:51:24 2023

@author: Arek
"""

import CoolProp.CoolProp as CP

class TurbineDesign:  #This class was originaly named turbine_design, but it was changed to TurbineDesign to be compliant to classes naming standards (convenction). This allows every programer to know, just by the name of it that this is a class and it can have methods (according to Complete Python Developer: Zero to Mastery course at Udemy)
    def __init__(self,T_in, p_in, p_out, mf,medium):
        self.T_in = T_in
        self.p_in = p_in
        self.p_out = p_out
        self.mf = mf
        self.medium = medium
        
    def h_in(self):
        """
        Calculates enthalpy at turbine outlet.
        """
        self.h_in = CP.PropsSI('Hmass', 'T', self.T_in+273.15, 'P', self.p_in*1000, self.medium)/1000
    
    def s_in(self):
        """
        Calculates entropy at turbine inlet and entrop
        """
        self.s_in = CP.PropsSI('Smass', 'T', self.T_in+273.15, 'P', self.p_in*1000, self.medium)/1000
        self.s_out_s = self.s_in
    
    def outlet_params(self, eta_iT):
        """
        Calculates temperature, entalpy and entropy at turbine outlet
        """
        pass
        
    def power(eta_iT):
        
        """
        Calculates internal power output of turbine according to given internal efficiency.
        
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
        X_out - vapour quality at turbine outlet 0.0-1.0, -
        
        Other outputs:
        h_in - enthalpy of vapour at turbine inlet, kJ/kg
        s_in - entropy of vapour at turbine inlet, kJ/(kg*C) or kJ/(kg*K)
        h_out - enthalpy of vapour at turbine outlet, kJ/kg 
        s_out - entropy of vapur at turbine outlet
        s_out_s - entropy of vapour at turbine outlet after ideal (eta = 1.0) expansion in turbine, kJ/(kg*C) or kJ/(kg*K)
        """
        self.N_iT = mf*(h_in) 
        return 1000

T_in = 210 # C
p_in = 845 # kPa(a)
p_out = 20 #kPa(a)
eta = 0.7  #kPa(a)
mf = 0.65
medium = 'Toluene'

turbina = TurbineDesign(T_in, p_in, p_out, mf, medium)