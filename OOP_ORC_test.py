import unittest
import time
from OOP_ORC import ORC
from scipy.optimize import least_squares, differential_evolution, shgo

# Setting flue gas composition
z_CO2 = 0.0670
z_O2  = 0.0610
z_H2O = 0.1310
z_N2  = 0.7410

T_fg_hot = 510 # Temperature of flue gas, st.C
T_fg_min = 125 # Minimal allowed temperature of flue gas to avoid condensation of water (based on OEM of gas fired engine), st.V
m_fg = 0.65    # Mass flow of flue gas, kg/s

p_max = 1500         # Maximal pressure in ORC system, kPa(a)
p_min = 20           # Minimal presurre in ORC system, kPa(a)
p_flue_gas = 101.325 # Pressure of heat source, at ORC heat exchenger inlet

T_sat_min = 45 # Minimal stauration temperature in condenser, st.C
T_coolant_in = 20 # Temperature of collant at condenser inlet
p_coolant_in = 101.325+300


# Design equipement parameters
dT_eco  = 40 # Minimal temperature difference in economizer, st.C
dT_evap = 80 # Minimal temperature difference in evaporator, st.C
dT_reg  = 50 # Minimal temperature difference in regenerator, st.C
dT_cond = 20 # Minimal temperature difference in condesnator, st.C

dT_sup = 5   # Superheating of working fluid in evaporator, st.C
dT_sub = 5   # Subcooliong of working fluid in condenser, st.C

eta_iT   = 0.75 # Isentropic efficiency of turbine, -
eta_iWFP = 0.60 # Isentropic efficeincy of working fluid pump, -
eta_iCP  = 0.60 # Isentropic efficiency of coolant pump, -
eta_iTOP = 0.45 # Isentropic efficiency of thermal oil pump, -

working_fluids = ['Toluene','MM','Benzene','m-Xylene','p-Xylene'] # Later it will be a list of fluids
flue_gas =  f'HEOS::CO2[{z_CO2}]&O2[{z_O2}]&N2[{z_N2}]&H2O[{z_H2O}]'
evap_flow = 'counter'
coolant = 'INCOMP::MEG-40%'

toluene_ORC1 = ORC(p_max,p_min, T_sat_min, 
                dT_evap,dT_cond,
                evap_flow,
                dT_sup, dT_sub, 
                T_fg_hot, T_fg_min, p_flue_gas, m_fg,
                T_coolant_in, p_coolant_in,
                flue_gas, 'Toluene', coolant,
                eta_iWFP = eta_iWFP, eta_iCP = eta_iCP,
                supercritical = False)

# MM_ORC1 =   ORC(p_max,p_min, T_sat_min, 
#                 dT_evap,dT_cond,
#                 evap_flow,
#                 dT_sup, dT_sub, 
#                 T_fg_hot, T_fg_min, p_flue_gas, m_fg,
#                 T_coolant_in, p_coolant_in,
#                 flue_gas, 'MM', coolant,
#                 eta_iWFP = eta_iWFP, eta_iCP = eta_iCP,
#                 supercritical = False)

def solve_ORC_by_pressure(p_max, ORC):
    
    ORC.p_max = p_max
    ORC.solve()
    N_net = ORC.N_net
    print(N_net)
    return 1/N_net

#test_ORC1.solve()
# start_time = time.time()
# print('least_squares')
# p_opt = least_squares(solve_ORC_by_pressure,x0 = p_max, bounds= (p_min,p_max), args = [toluene_ORC1]).x

# end_time = time.time()

# lest_squares_time =  end_time - start_time

# hours, rem = divmod(lest_squares_time, 3600)
# minutes, seconds = divmod(rem, 60)
# print("{:0>2}:{:0>2}:{:05.2f}".format(int(hours),int(minutes),seconds))



# start_time = time.time()

# p_opt = shgo(func=solve_ORC_by_pressure, bounds= [(p_min,p_max)], args = (toluene_ORC1,)).x

# end_time = time.time()

# shgo_time = end_time - start_time

# hours, rem = divmod(shgo_time, 3600)
# minutes, seconds = divmod(rem, 60)
# print("{:0>2}:{:0>2}:{:05.2f}".format(int(hours),int(minutes),seconds))

start_time = time.time()

p_opt = differential_evolution(solve_ORC_by_pressure, [(p_min,p_max)], args = [toluene_ORC1], seed = 42)

end_time = time.time()

diff_evo_time = end_time - start_time          

hours, rem = divmod(end_time-start_time, 3600)
minutes, seconds = divmod(rem, 60)
print("{:0>2}:{:0>2}:{:05.2f}".format(int(hours),int(minutes),seconds))

T_pump_out = toluene_ORC1.T_pump_out
T_sat_cond = toluene_ORC1.T_sat_cond
T_sup = toluene_ORC1.T_sup
mf = toluene_ORC1.mf
dT_evap_calc = toluene_ORC1.dT_evap_calc
eta_iT = toluene_ORC1.eta_iT
N_iT = toluene_ORC1.N_iT
T_tubine_out = toluene_ORC1.T_turbine_out
mc = toluene_ORC1.mc
T_cond_cold_out = toluene_ORC1.T_cond_cold_out
N_iCP = toluene_ORC1.N_iCP
N_ORC = toluene_ORC1.N_net
K_T = toluene_ORC1.K_T
temp_dist_diff = toluene_ORC1.EVAP.temp_dist_df
#print(toluene_ORC1.EVAP.temp_dist_df)
toluene_ORC1.EVAP.plot_temp_dist_df()
toluene_ORC1.COND.plot_temp_dist_df()

