# -*- coding: utf-8 -*-
"""
Created on Wed Apr  5 11:51:24 2023

@author: Arek
"""
import math
import CoolProp.CoolProp as CP

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
    
    def __init__(self,T_in, p_in, p_out, mf,medium):
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
        
    def efficiency_estimation(self):
        
        """
        Simplified, preliminary estimation of design turbine internal efficiency. Change in inlet and outlet parameters affects value of efficiency.
        
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
        """

        self.rho_turbine_in = CP.PropsSI('Dmass', 'T', self.T_in+273.15, 'P', self.p_in*1000, self.medium)
        self.rho_turbine_out_s = CP.PropsSI('Dmass', 'Smass', self.s_out_s*1000, 'P', self.p_out*1000, self.medium)
        self.delta_h_s = self.h_in-self.h_out_s
        
        self.V_turbine_in = self.mf/self.rho_turbine_in
        self.V_turbine_out_s =self.mf/self.rho_turbine_out_s
        
        self.SP = (self.V_turbine_out_s**(1/2))/(self.delta_h_s**(1/4))

        self.V_r = self.V_turbine_out_s/self.V_turbine_in
        
        self.eta_iT = (1 * 0.9090831500) + \
             (math.log(self.SP)*(-0.05248690)) + \
                 (((math.log(self.SP))**2)*(-0.04799080)) + \
                    (((math.log(self.SP))**3)*(-0.01710380)) + \
                         (((math.log(self.SP))**4)*(-0.00244002)) + \
                             (math.log(self.V_r)*(0.04961780)) + \
                             ((math.log(self.V_r)**2)*(-0.04894860)) + \
                                 ((math.log(self.V_r)**3)*(0.01171650)) + \
                                     ((math.log(self.V_r)**4)*(-0.00100473)) + \
                                         (math.log(self.V_r)*math.log(self.SP)*0.05645970) + \
                                            ((math.log(self.V_r)**2)*math.log(self.SP)*(-0.01859440)) + \
                                                 (math.log(self.V_r)*(math.log(self.SP)**2)*(0.01288860)) + \
                                                     ((math.log(self.V_r)**3)*math.log(self.SP)*(0.00178187)) + \
                                                         ((math.log(self.V_r)**3)*(math.log(self.SP)**2)*(-0.00021196)) + \
                                                            ((math.log(self.V_r)**2)*(math.log(self.SP)**3)*(0.00078667))
                                                                
        self.h_out = self.h_in-self.eta_iT*(self.h_in-self.h_out_s)
        self.T_out = CP.PropsSI('T', 'Hmass', self.h_out*1000, 'P', self.p_out*1000, self.medium)-273.15
        self.N_iT = self.mf*(self.h_in-self.h_out)
        
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
        X_out - vapour quality at turbine outlet 0.0-1.0, -
        
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
        self.N_iT = self.mf*(self.h_in-self.h_out)
        
    def plot_Ts():
        pass
        
    def off_design(self, T_in_off, p_in_off, mf_off, p_out_off):
        pass
    
        #m_off  = m_nom*(((p1_off**2-p2**2)/(p1_nom**2-p2_nom**2))*(p1_nom*ve1_nom)/(p1_off*ve1_off))**0.5
        # def turbine_off(eta_off,eta_nom, T1_off, p1_off, p2,  m_off, medium):
        
        # h1_off = CP.PropsSI('Hmass','T',T1_off + 273.15,'P',p1_off  * 1000, medium)/1000                #kJ/kg
        # s1_off = CP.PropsSI('Smass', 'P',p1_off  * 1000, 'T', T1_off + 273.15, medium)/1000             #kJ/(kg*K)
        # h2s_off = CP.PropsSI('Hmass', 'P',p2  * 1000, 'Smass', s1_off * 1000, medium)/1000          #kJ/(kg*K)
        
        # h1_nom = CP.PropsSI('Hmass','T',211.8 + 273.15,'P',845  * 1000, medium)/1000
        # h2_nom = CP.PropsSI('Hmass','T',136.7 + 273.15,'P',24  * 1000, medium)/1000
        # dh_nom = h1_nom - h2_nom
        
    
        # rho2_nom = CP.PropsSI('Dmass','T',136.7 + 273.15,'P',24  * 1000, medium)
        # V2_nom = 0.388/rho2_nom
        
        # a = 0.248
        # b = 1.632
        # c = -1.940
        # d = 0.033
        # e = -1.085
        # f = 2.112
        
        # h2_off = h1_off - eta_off*(h1_off-h2s_off)
                
        # #Q2_off = CP.PropsSI('Q','Hmass',h2_off * 1000 ,'P',p2_nom  * 1000, medium)
        
        # T2_off = CP.PropsSI('T','Hmass',h2_off * 1000,'P',p2  * 1000, medium) - 273.15 #C
        # rho2_off =  CP.PropsSI('Dmass','T',T2_off + 273.15,'P',p2  * 1000, medium)   #kg/m^3
        # V2_off = m_off/rho2_off
                    
        # dh_off = h1_off - h2_off
                    
        # eta_off_calc = (a + b * (dh_off/dh_nom) + c * (dh_off/dh_nom)**2 + d*(V2_off/V2_nom) + e*(V2_off/V2_nom)**2 + f*(dh_off/dh_nom)*(V2_off/V2_nom))*eta_nom
                
        # d_eta = math.fabs(eta_off - eta_off_calc)

T_in = 210 # C
p_in = 845 # kPa(a)
p_out = 20 #kPa(a)
eta = 0.7  #kPa(a)
mf = 0.65
medium = 'Toluene'

turbina = TurbineDesign(T_in, p_in, p_out, mf, medium)
turbina.efficiency_estimation()
eta_iT = turbina.eta_iT
print(eta_iT)