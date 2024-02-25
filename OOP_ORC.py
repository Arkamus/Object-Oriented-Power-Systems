import math
import CoolProp.CoolProp as CP
from scipy.optimize import least_squares
from OOPS import Turbine, Pump, HeatExchanger

class ORC:
    
    def __init__(self, 
                 p_max, p_min, T_sat_min,
                 dT_evap, dT_cond,
                 evap_arangement, 
                 dT_sup, dT_sub, 
                 T_heat_source, T_heat_source_min, p_heat_source, m_heat_source,
                 T_coolant_in,  p_coolant_in,
                 heat_source_carrier, working_fluid, coolant,
                 eta_iWFP, eta_iCP,
                 supercritical): #, m_fg, T_fg_hot, T_fg_min, z_CO2, z_O2,z_N2,z_H2O, dT_evap):
        
        """
        Simple ORC layout - no intermediate loop, no regeneration.
        
        p_max           - Maximal allowed pressure in ORC system, kPa(a)
        p_min           - Minimal allowerd presurre in ORC system, kPa(a)
        p_fluid_max     - Maximal pressure of fluid with know properties, kPa(a)
        

        T_fluid_max     - Maximal temperature of fluid with known properties, C.deg
        
        dT_sup          - Superheating of working fluid in evaporator, C.deg
        dT_sub          - Subcooliong of working fluid in condenser, C.deg
        
        T_sat_min       - Minimal stauration temperature in condenser, C.deg
        T_sat1          - Saturation temperature of working fluid at pressure p_max, C.deg
        T_sat2          - Saturation temperature of working fluid at pressure p_min, C.deg
        T_sup           - Temperature of fluid vapour after superheating, C.deg
        
        supercritical   - True or False. If true, than supercritical parameters of fluid are allowed (aboce critical point)
        """
        
        'Verifying working fluid properties accorcindg to input data - p_max, p_min, p_crit, T_crit'
        

        
        self.p_max = p_max
        self.p_min = p_min
        self.T_sat_min = T_sat_min
        self.dT_evap = dT_evap
        self.dT_cond = dT_cond
        self.evap_arangement = evap_arangement
        self.dT_sup = dT_sup
        self.dT_sub  = dT_sub
        self.T_heat_source = T_heat_source
        self.T_heat_source_min = T_heat_source_min
        self.p_heat_source = p_heat_source
        self.m_heat_source = m_heat_source
        self.T_coolant_in = T_coolant_in 
        self.p_coolant_in = p_coolant_in
        self.heat_source_carrier = heat_source_carrier
        self.working_fluid = working_fluid
        self.coolant = coolant
        self.eta_iWFP = eta_iWFP
        self.eta_iCP = eta_iCP
        self.supercritical = supercritical
        
        
    def solve(self):
        
        def solve_by_mass(guess, dT_min, heat_exchanger, he_type):
        
            """
            Calcualtion of correct working fluid mass flow (mf) and outlet temperature of coler (T_cold_out) or hotter fluid (T_hot_out)
            
            mf              - mass floww of working fluid, kg/s
            dT_min          - desired design value of minimal temperature difference in heat exchanger, C
            heat_exchanger  - object representing heat exchanger in ORC system to be updated during calculations in function
            he_type         - evap or condenser
            """
            
            if he_type == 'evap':
            
                heat_exchanger.T_hot_out = guess[0]
                heat_exchanger.m_cold = guess[1]
                
            elif he_type =='cond':
                
                heat_exchanger.T_cold_out = guess[0]
                heat_exchanger.m_cold = guess[1]
                
            heat_exchanger.calc_temp_dist(points = 4)
            dQ = heat_exchanger.Q_hot - heat_exchanger.Q_cold
            dT = heat_exchanger.dT_min - dT_min
            
            roots = [dQ, dT]
            
            return roots
        
        
        "Solving ORC with given parameters"
    
        if self.supercritical == False:
            
            self.T_crit = CP.PropsSI('Tcrit', self.working_fluid)-273.15
            self.p_crit = CP.PropsSI('Pcrit', self.working_fluid)/1000
            
            if self.p_max > self.p_crit:
                
                self.p_max = self.p_crit
        
        self.T_sat_evap = CP.PropsSI('T','P', self.p_max*1000, 'Q', 1, self.working_fluid) - 273.15
        self.T_sup = self.T_sat_evap+self.dT_sub
        
        self.T_fluid_max = CP.PropsSI('Tmax', self.working_fluid)-273.15
        self.p_fluid_max = CP.PropsSI('pmax', self.working_fluid)/1000
        
        'Checking if p_min is allowerd for specified working fluid and T_sat_min'
        
        self.T_sat_cond = CP.PropsSI('T','P', self.p_min*1000, 'Q', 1, self.working_fluid) - 273.15
        self.T_sub = self.T_sat_cond-self.dT_sub
        
        if self.T_sat_cond<self.T_sat_min:
            self.p_cond =  CP.PropsSI('P','T', self.T_sat_min+273.15, 'Q', 1, self.working_fluid)/1000
        else:
            self.p_cond = self.p_min
        
        'Calculation of temperature at evaporator inlet (pump outlet) by creating working fluid pump'
        medium_pump = Pump(T_in=self.T_sub,p_in = self.p_cond, p_out = self.p_max, mf =1, medium = self.working_fluid) #This have to change. Pump must have option to calculate temperature at outlet without specifying the mass flow (it is not required in this aproach)
        medium_pump.calc_T_out(self.eta_iWFP)
        self.T_pump_out = medium_pump.T_out
        
        'Calculating mass flow using heat exchnager class'
        self.EVAP = HeatExchanger(T_hot_in = self.T_heat_source, T_hot_out = self.T_heat_source_min, T_cold_in = self.T_pump_out, T_cold_out = self.T_sup, 
                                      m_hot = self.m_heat_source, m_cold = self.m_heat_source, p_hot_in = self.p_heat_source, p_cold_in = self.p_max, 
                                      hot_fluid = self.heat_source_carrier, cold_fluid = self.working_fluid, 
                                      flow_arrangement = self.evap_arangement)
        
        self.T_evap_hot_out, self.mf = least_squares(fun=solve_by_mass, x0=[self.T_heat_source,self.m_heat_source], bounds = ([self.T_heat_source_min,0],[self.T_heat_source,math.inf]), args = ([self.dT_evap,self.EVAP, 'evap'])).x
        
        self.dT_evap_calc = self.EVAP.dT_min
        
        'Recalculation of working fluid pump power'
        medium_pump.mf = self.mf
        medium_pump.power(self.eta_iWFP)
        self.N_iWFP = medium_pump.N_iP
        
        'Turbine calculation'
        turbine = Turbine(self.T_sup, self.p_max, self.p_cond, self.mf, self.working_fluid)
        turbine.efficiency_estimation(n_stages= 1)
        self.eta_iT = turbine.eta_iT
        
        turbine.power(self.eta_iT)
        
        self.N_iT = turbine.N_iT
        self.T_turbine_out = turbine.T_out
        
        'Condenser - finding coolant mass flow and temperature at condenser outlet'
        self.COND = HeatExchanger(T_hot_in = self.T_turbine_out, T_hot_out = self.T_sub, T_cold_in = self.T_coolant_in, T_cold_out = self.T_coolant_in, 
                             m_hot = self.mf, m_cold = self.mf, p_hot_in = self.p_cond, p_cold_in = self.p_coolant_in, 
                             hot_fluid = self.working_fluid, cold_fluid = self.coolant, 
                             flow_arrangement = 'counter')
        
        'Precalculation of guess values for condenser'
        self.T_cond_cold_out, self.mc = least_squares(fun=solve_by_mass, x0=[40,90], bounds = ([self.T_coolant_in,0],[90.0,math.inf]), args = ([self.dT_cond,self.COND, 'cond'])).x
        
        'Coolant pump power consumption'
        coolant_pump = Pump(T_in=self.T_coolant_in,p_in = self.p_coolant_in-300, p_out = self.p_coolant_in, mf =self.mc, medium = self.coolant) #This have to change. Pump must have option to calculate temperature at outlet without specifying the mass flow (it is not required in this aproach)
        coolant_pump.power(self.eta_iCP)
        self.N_iCP = coolant_pump.N_iP
        
        'Net power output'
        self.N_net = self.N_iT-self.N_iWFP - self.N_iCP
        
        self.K_T = turbine.cost()
        
        self.EVAP.temp_dist_df
        self.COND.temp_dist_df
        
        

