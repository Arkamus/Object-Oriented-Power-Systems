# -*- coding: utf-8 -*-
"""
Created on Thu Aug 24 19:06:36 2023

@author: Arek
"""

import unittest
import numpy as np
import matplotlib.pyplot as plt
import CoolProp.CoolProp as CP

from OOP_pump import Pump

class TestOOPPump(unittest.TestCase):
    
    def setUp(self):
        
        """
        Test parameters for pump 
        Test parameters are obtained from:
            Ahlgren F, Mondejar ME, Genrup M, Thern M. Waste Heat Recovery in a Cruise Vessel
            in the Baltic Sea by Using an Organic Rankine Cycle: A Case Study. Journal of
            Engineering for Gas Turbines and Power 2016;138:11702, doi:10.1115/1.4031145
            
        It is possible to obtain more test cases form abovementioned publication.
        """
        
        self.dT_sub  = 5.0       # Subcoling temperature in conderser, C
        self.p_in    = 5.0       # Design/test inlet pressure, kPa(a)
        self.p_out   = 1446.0    # Design/test outlet pressure, kPa(a)
        self.mf      = 3.296     # Design/test mass flow, kg/s
        self.medium  = 'Toluene' # Design/test medium
        self.eta     = 0.8       # Design/test pump efficiency
        
        self.T_in   = CP.PropsSI('T', 'P', self.p_in*1000, 'Q', 0, self.medium)-self.dT_sub-273.15     # Design/test inlet temperature, C
        print(self.T_in)
        
    def test_power(self):
        """
        Drawing test T-s diagrams for test case.
        """
        
        print('Runing test T-s diagrams')

        test_pump1 = Pump(self.T_in, self.p_in, self.p_out, self.mf, self.medium)
        test_pump1.power(self.eta)
        T_out = test_pump1.T_out

    def test_case_off_design1(self):
        
        print('Runing off-design test 1 - ')
        test_pump1 = Pump(self.T_in, self.p_in, self.p_out, self.mf, self.medium)
        test_pump1.power(self.eta)
        
        init_mf_off_ar   = np.arange(2.0,4.687,0.001)
        init_p_out_off_ar = np.arange(9.96*100,18.72*100, 1)
        
        N_iP_off_ar = []
        eta_iP_off_ar = []
        mf_off_ar = []
        p_out_off_ar = []
        
        for mf_off in init_mf_off_ar:
            test_pump1.off_design(self.T_in, self.p_in, self.p_out, mf_off)
            N_iP = test_pump1.N_iP
            eta_iP = test_pump1.eta_iP
            
            N_iP_off = test_pump1.N_iP_off
            eta_iP_off = test_pump1.eta_iP_off
            
            if eta_iP_off>=0:
                mf_off_ar.append(mf_off)
                N_iP_off_ar.append(N_iP_off)
                eta_iP_off_ar.append(eta_iP_off)
          
        fig1, axs = plt.subplots(nrows=2, ncols=1, layout = 'constrained', dpi = 400)
        
        ((ax1, ax2)) = axs
        
        color1 = 'blue'
        
        ax1.plot(mf_off_ar, N_iP_off_ar, c= color1)
        ax1.scatter(self.mf,N_iP, marker = 'x', c =color1)
        
        xticks = np.arange(2.0,5.0,0.25)
        ax1.set_xticks(xticks)
        ax1.set_xlim(2.0,4.75)
        ax1.set_xlabel('Off-design mass flow $mf^{off}$, kg')
        
        ax1.set_ylim(4.0,14.0)
        ax1.set_ylabel(r'$N_{iP}^{off}=f(mf^{off})$', c = color1)
        
        ax1.grid(axis ='both', linestyle = ':')
        ax1.tick_params(axis = 'y', labelcolor = color1)
        ax1.set_title(f'{self.medium}')
        
        color2 = 'red'
        
        ax1a = ax1.twinx()
        
        ax1a.plot(mf_off_ar, eta_iP_off_ar, c = color2)
        ax1a.scatter(self.mf,eta_iP, marker = 'x', c = color2)
        
        ax1a.set_ylim(0.0,1.0)
        ax1a.set_ylabel(r'$\eta_{iP}^{off} = f(mf^{off})$', c = color2)
        ax1a.tick_params(axis = 'y', labelcolor = color2)
        
        N_iP_off_ar = []
        eta_iP_off_ar = []
        p_out_off_ar = []
        
        for p_out_off in init_p_out_off_ar:
            test_pump1.off_design(self.T_in, self.p_in, p_out_off, self.mf)
            N_iP   = test_pump1.N_iP
            eta_iP = test_pump1.eta_iP
            
            N_iP_off   = test_pump1.N_iP_off
            eta_iP_off = test_pump1.eta_iP_off
            
            print(p_out_off, N_iP_off, self.mf, eta_iP_off)
            if eta_iP_off>=0:
                p_out_off_ar.append(p_out_off)
                N_iP_off_ar.append(N_iP_off)
                eta_iP_off_ar.append(eta_iP_off)
        
        ax2.plot(p_out_off_ar, N_iP_off_ar, c = color1)
        ax2.scatter(self.p_out,N_iP, marker = 'x', c = color1)
        
        xticks = np.arange(900,2000,100)
        ax2.set_xticks(xticks)
        ax2.set_xlim(900,1900)
        ax2.set_xlabel('Off-design outlet pressure $p_{out}^{off}$, kPa(a)')
        
        ax2.set_ylim(4.0,14.0)
        ax2.set_ylabel(r'$N_{iP}^{off}=f(mf^{off})$', c = color1)
        
        ax2.grid(axis ='both', linestyle = ':')
        ax2.tick_params(axis = 'y', labelcolor = color1)
            
if __name__=='__main__':
 	
    unittest.main()