import numpy as np

"""
Example cost functions for equipement of ORCs with strictly defined validity range.

Sources:
    1. Self-gathered data
    2. 
"""

def PLCCS_cost(NOM, a=1.22583568e+04, b=6.98958079e-01,c=-4.63641907e+04,low_bound = 45, high_bound = 64):
    """
    Total cost of PLC driver and overal control system according to number of measurements.
    
    Inputs:
        NOM - numbers of individual measurements (temperature, pressure, electric signals, eletric signals for control)
        a, b, c - default or user defined cost coeffictients for cost function, -
        low_bound - default or user defined lower boundary of NOM
        high_bound - default or user defined higher boundary of NOM
    
    Output:
        cost of PLC driver and related control systems in PLN
    """
    
    if NOM < low_bound or NOM > high_bound:
        raise ValueError(f'Not valid outside range of {low_bound}-{high_bound} NOM')
    
    return a*NOM**b+c

def ups_cost(N_net, a=183.43987295, b=29906.36344314, low_bound = 32.199, high_bound = 246.4):
    """
    Cost of emergency energy source for safety systems in case of failure as function of ORC net power output.
    
    Inputs:
        N - nett ORC power output, kWel
        a, b - default or user defined cost coefficients of UPS
        low_bound - default or user defined lower boundary of N_net
        high_bound - default or user defined higher boundary of N_net
        
    Output:
        cost of UPS in PLN
    """
    
    if N_net < low_bound or N_net > high_bound:
        raise ValueError(f"Not valid in range of {low_bound}-{high_bound} kWel of ORC electric net power output")
        
    return a*N_net+b

def glycol_amount(Q_cond, a=2.54016375e+02, b=2.81242000e-01, c=-2.22153727e-27, low_bound = 0.0, high_bound = 1889.93):
    """
    Amount of cooling loop as function of heat distpatched in condenser.
    The function is calibrated with 40% water-glycol mixture (INCOMP::MEG-40% using CoolProp thermophysical properties library)
    
    Inputs:
        Q_cond - heat transfered in condenser, kWt
        a, b, c - default of user defined coeficents of cooling loop capacity
        low_bound - default or user defined lower boundary of Q_cond
        high_bound - default or user defined higher boundary of Q_cond
        
    Output:
        capcaity of colling loop in kg
    """
    
    if Q_cond <= low_bound or Q_cond > high_bound:
        raise ValueError(f"Not valid outside range of {low_bound}-{high_bound} kWt")
    
    return a*Q_cond**b+c

def glycol_cost(m, a=-7.97844951e+03, b=7.99166251e+03, c=9.99837733e-01, d=-1.40552932e+02, low_bound = 0.0, high_bound = 2120.0):
    """
    Cost of glycol mixture as function of mass flow.
    
    Inputs:
        m - mass flow, kg/s
        a, b, c, d - default or user defined coeficents of cost
        low_bound - default or user defined lower boundary of m
        high_bound - default or user defined higher boundary of m
    
    Output:
        cost of glycol in PLN
    """
    
    if m < low_bound or m > high_bound:
        raise ValueError("Not valid outside range of 0-2120.0 kg")
        
    return a*m+b*m**c+d

def oil_amount(Q_source, a=38.55960845, b=0.57066217, low_bound = 0.0, high_bound = 2500.0):
    """
    Amount of thermal oil required for intermediate loop between heat source and ORC.
    
    Inputs:
        Q_source - heat harvested from heat source, kWt
        a, b - default or user defined coefficients of oil amount
        low_bound - default or user defined lower boundary of Q_source
        high_bound - default or user defined higher boundary of Q_source
    
    Output:
        amount of oil in intermediate loop, kg
    """
    
    if Q_source<low_bound or Q_source>high_bound:
        raise ValueError(f"Not valid outsie range of {low_bound}-{high_bound} kWt")  
    
    return a*Q_source**b

def oil_cost(OA,a=1.34447159e+01,b=1.14638131e+00,c=2.71741773e-23, low_bound = 0.0, high_bound = 3120.0):
    """
    Estimation of thermal oil cost.
    
    Inputs:
        OA - amount of thermal oil, kg
        a, b, c, d - default or user defined coefficents for cost of thermal oil
        low_bound - default or user defined lower boundary of OA, kg
        high_bound - default or user defined higher boundary of OA, kg
        
    Output:
        cost of oil, PLN
    """
    
    if OA < low_bound or OA > high_bound:
        raise ValueError(f"Not valid outsie range of {low_bound}-{high_bound} kg")  
    
    return a*OA**b+c

def container_cost(N_iT,a = 1.00002265, b = 0.13550264, low_bound = 8.5, high_bound = 335.75):
    """
    Cost of cntainer for ORC
    
    Inputs:
        N_iT - turbine power, kWel
        a, b - default or user defined ceofficents for container cost
        low_bound - default or user defined lower boundary of N_iT, kg
        high_bound - default or user defined higher boundary of N_iT, kg
        
    Output:
        cost of container, PLN
    """
    
    if N_iT < low_bound or N_iT > high_bound:
        raise ValueError("Not valid outsie range of {low_bound}-{high_bound} kW") 
    
    return np.log(b*N_iT)/np.log(a)

def working_fluid_amount(N_iT, a=3.21189450e+02, b=4.14529560e-01, c=-5.79889171e+02, low_bound = 8.5, high_bound = 335.75):
    """
    Amount of thermal oil required for intermediate loop between heat source and ORC.
    
    Inputs:
        Q_source - heat harvested from heat source, kWt
        a, b - default or user defined coefficients of oil amount
        low_bound - default or user defined lower boundary of Q_source
        high_bound - default or user defined higher boundary of Q_source
    
    Output:
        amount of oil in intermediate loop, kg
    """
    
    if N_iT < low_bound or N_iT > high_bound:
        raise ValueError("Not valid outsie range of {low_bound}-{high_bound} kW") 
    
    return a*N_iT**b+c

def working_fluid_cost(m,fluid):
    """
    
    Inputs:
        m - amount of working fluid, kg
    
    Outputs:
        cost of fluid, PLN
    """
    if fluid == 'Toluene':
          a = 28.06405694
          b = 160.7619703
          c = 0.98033196
    elif fluid == 'CycloHexane':
          a = 44.36708842
          b = 232.06086828
          c = 0.97992304
    elif fluid == 'p-Xylene' or fluid == 'm-Xylene' or fluid == 'o-Xylene':
          a = 32.32102097
          b = 195.47108813
          c = 0.97908132
    elif fluid == 'Benzene':
         a = 43.91778304
         b = 255.14015112
         c = 0.98037581
    elif fluid == 'MM':
        a = 1.98635886e+03
        b = -9.31434506e+02
        c = 6.51198175e-01
    elif fluid == 'R1233zd(E)':
        return 186.04212307823948*m
    else:
        raise Exception('Medium not in list of avaliable fluids')
    return a*m**c+b

def funkcja_kosztu_turbiny(X,a = 1.0008018e+06,b = 1.4400000e-01,c = -1.3782345e+06):
    """
    X -Moc turbiny
    """
    return a*X**b+c

def koszt_generatora(N,a = 375.84657759,b = 1.22678157):
    """
    N - moc turbiny
    """
    return a*(N)**b

def koszt_pomp(N, p, K1 = 3.4528671, K2 = 0.86809638, K3 = -0.09892602):
        
        """
        C - korekta ze względu na ciśnienie
        K - Korekta ze względu na moc
        """

        if p<10:

            C1, C2, C3 = 0,0,0

        else:

            C1 = -0.3935
            C2 = 0.3957
            C3 = -0.00226

        B1 = 1.89
        B2 = 1.35

        Fm = 1.50

        CP = 10**(K1+K2*np.log10(N)+K3*(np.log10(N))**2) # człon ceny modułu nie uwzględniający podwyższonego ciśnienia, a jdeynie podstawowy parametr rozmiaru
        Fp = 10**(C1+C2*np.log10(p)+C3*(np.log10(p))**2) # korekta ze względu na podwyższone ciśnienia
        CBM = CP*(B1+B2*Fm*Fp)                            # uwzględnienie kosztów materiałowych

        return CBM

def koszt_par_skr(X, z1 = 1.52579904,z2 = 3.00280577,z3 = -0.63071587):
    """
    Funkcja kosztów parowników i skraplaczy
    f - powierzchnia wymiany ciepła w m^2
    p - nadcisnienie, bar
    """
    f = X[0]
    p = X[1]

    if p < 5.0:   
        C1, C2, C3 = 0,0,0
    else:
        C1 = 0.03881
        C2 = -0.11272
        C3 = 0.08183

    K1 = z1 #4.8306
    K2 = z2 #-0.8509
    K3 = z3 #0.3187

    B1 = 1.63
    B2 = 1.66
    Fm = 1.30

    CP = 10**(K1+K2*np.log10(f)+K3*(np.log10(f))**2)
    Fp = 10**(C1+C2*np.log10(p)+C3*(np.log10(p))**2)
    CBM = CP*(B1+B2*Fm*Fp)

    return CBM

def koszt_eco_reg(X, z1=2.70218796, z2=2.73148019, z3 = -0.5472148, z4 = -1.75237871):#z1 = 2.58688192,z2 = 2.56867884,z3 = -0.50953362,z4 = -5.83094089,z5 = 2.83869729,z6 = 2.33118288):
    
    """
    Funkcja kosztu ekonomizerów i regeneratorów
    F - pole powierzchni, m^2
    P - nadcisnienie, bar
    """
    
    f = X[0]
    p = X[1]
    
    if p < 5.0:   
        C1, C2, C3 = 0,0,0
    else:
        C1 = 0.03881
        C2 = -0.11272
        C3 = 0.08183

    K1 = z1
    K2 = z2
    K3 = z3

    B1 = z4
    B2 = 1.66
    Fm = 1.30

    CP = 10**(K1+K2*np.log10(f)+K3*(np.log10(f))**2)
    Fp = 10**(C1+C2*np.log10(p)+C3*(np.log10(p))**2)
    CBM = CP*(B1+B2*Fm*Fp)

    return CBM

def koszt_reg(X, z1=3.56953291, z2=1.06179318, z3 = -0.15206137, z4 = 1.633511):#z1 = 2.58688192,z2 = 2.56867884,z3 = -0.50953362,z4 = -5.83094089,z5 = 2.83869729,z6 = 2.33118288):
    
    """
    Funkcja kosztu ekonomizerów i regeneratorów
    F - pole powierzchni, m^2
    P - nadcisnienie, bar
    """
    
    f = X[0]
    p = X[1]
    
    if p < 5.0:   
        C1, C2, C3 = 0,0,0
    else:
        C1 = 0.03881
        C2 = -0.11272
        C3 = 0.08183

    K1 = z1
    K2 = z2
    K3 = z3

    B1 = z4
    B2 = 1.66
    Fm = 1.30

    CP = 10**(K1+K2*np.log10(f)+K3*(np.log10(f))**2)
    Fp = 10**(C1+C2*np.log10(p)+C3*(np.log10(p))**2)
    CBM = CP*(B1+B2*Fm*Fp)

    return CBM

def oil_tank_volume(Q_source,a=0.00128034,b=0.90043023):
    """
    Thermal oil tank capacity, as fucntion of heat transfered (harvested) from heat source.
    Default coefficents are valid for Therminol66
    
    Inputs:
        Q_source - heat harvested from heat source, kWt
        a,b - default or user defined coefficents of oil loop capacity
        
    Output:
        Oil tank volume, m^3
    """
    
    if Q_source<0.0 or Q_source>2500.0:
        raise ValueError("Not valid outsie range of 0.0-2500.0 kWt")    
        
    return a*Q_source**b

def fluid_tank_volume(N,a=0.1103605,b=-0.11168448,c=0.98588005,d=0.14599678):
    """
    Working fluid tank volume as a function of difference between turbine power and working fluid pump power consumption
    N - difference between turbine power output and working fluid pump power consumption, kWel
    """
    return a*N+b*N**c+d

def tank_cost(V,a=1.17672243e+04,b=4.48970830e-03,c=1.56753725e+01,d=2.49338070e+04):
    """
    Cost of working medium tank as function of tank volume.
    
    Inputs:
        V - tank volume, m^3
        a, b, c, d - 
    """
    return a*V+b*V**c+d
