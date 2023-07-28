# -*- coding: utf-8 -*-
"""
Created on Wed Jul 19 01:45:43 2023

@author: Arek
"""

import unittest
from OOP_turbine import TurbineDesign

class TestOOPTurbine(unittest.TestCase):
    
    # def test_power(self):
    #     test_T_in = 210 # C
    #     test_p_in = 845 # kPa(a)
    #     test_p_out = 20 #kPa(a)
    #     test_eta_iT = 0.7  #kPa(a)
    #     test_mf = 0.65
    #     test_medium = 'Toluene'
        
    #     test_turbine = TurbineDesign(test_T_in, test_p_in, test_p_out, test_mf, test_medium)
    #     test_turbine.power(test_eta_iT)
        
    def test_estimated_efficiency(self):
        
        """
        This test is based on data from:
            
        Macchi E, Astolfi M. Organic rankine cycle (ORC) power systems: Technologies and
        applications. Duxford, United Kingdom: Woodhead Publishing is an imprint of Elsevier;
        2017.
        """
        
        print('Runing test efficiency estimations')

        test_turbine1 = TurbineDesign(155, 3620.0, 1568.5, 237.71, 'R125')
        test_turbine1.efficiency_estimation(n_stages= 1)
        test_eta_iT1 = test_turbine1.eta_iT
        
        test_turbine2 = TurbineDesign(155, 3620.0, 1568.5, 11.89, 'R125')
        test_turbine2.efficiency_estimation(n_stages = 1)
        test_eta_iT2 = test_turbine2.eta_iT
    
        test_turbine3 = TurbineDesign(234.67, 829.0,25.0, 40.78, 'Hexane')
        test_turbine3.efficiency_estimation(n_stages = 1)
        test_eta_iT3 = test_turbine3.eta_iT
    
        test_turbine4 = TurbineDesign(234.67, 829.0,25.0, 2.04, 'Hexane')
        test_turbine4.efficiency_estimation(n_stages = 1)
        test_eta_iT4 = test_turbine4.eta_iT
        
        self.assertEqual(test_eta_iT1, 0.9028177495017103)
        self.assertEqual(test_eta_iT2, 0.8770645657764151)
        self.assertEqual(test_eta_iT3, 0.835051048808327)
        self.assertEqual(test_eta_iT4, 0.8065933042772973)
        
if __name__=='__main__':
	
    unittest.main()