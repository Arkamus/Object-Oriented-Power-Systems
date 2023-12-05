import unittest

from OOP_heat_exchanger import HeatExchanger

class TestOOP_heat_exchanger(unittest.TestCase):

    """
    Test case 0 - check if flow_arrangement has 
    """
    def test_ValueError(self):

      with self.assertRaises(Exception) as context:
          HeatExchanger(T_hot_in = 510, T_hot_out = 222.24, T_cold_in = 182.24, T_cold_out = 310, 
                                   m_hot = 0.711582, m_cold = 0.719, p_hot_in = 101.325, p_cold_in = 300, 
                                   hot_fluid = 'Air', cold_fluid = 'INCOMP::T66', 
                                   flow_arrangement = 'kąter')
    
    
    def test_HeatExchanger(self):
        
        """
        Test case 1 - heat transfer in counter flow without phase change eg. gas-to-gas, liquid-to-liquid, gas-to-liquid or liquid-to-gas. 
        In other words there is no phase change. In this test case the heat exchange occurs between two predefined fluids: air and thermal oil Therminol66.
        In this case air composition is constant, is predefined in CoolProp and does not require to specify content of oxygen and nitrogen.
        """
        
        test_HE1 = HeatExchanger(T_hot_in = 510, T_hot_out = 222.24, T_cold_in = 182.24, T_cold_out = 310, 
                                  m_hot = 0.7118, m_cold = 0.719, p_hot_in = 101.325, p_cold_in = 300, 
                                  hot_fluid = 'Air', cold_fluid = 'INCOMP::T66', 
                                  flow_arrangement = 'counter')
        
        
        test_HE1.calc_temp_dist(points= 10)
        test_HE1_fig, test_HE1_ax = test_HE1.plot_temp_dist_df()
        test_HE1_ax.set_title('Test HE1')
        
        self.assertEqual(round(test_HE1.dT_min,2), 40.0)
        
        """
        Test case 2 - heat transfer in counter flow without phase change eg. gas-to-gas, liquid-to-liquid, gas-to-liquid or liquid-to-gas.
        In other words there is no phase change. In this case heat transfer occurs between user defined mixture (usually exchaust gases from various proceses) 
        and predefined heat transfer oil Therminol66 
        """
        z_CO2 = 0.0670 		
        z_O2  = 0.0610	
        z_H2O = 0.1310
        z_N2  = 0.7410
        
        test_HE2 = HeatExchanger(T_hot_in = 510, T_hot_out = 222.24, T_cold_in = 182.24, T_cold_out = 310, 
                                  m_hot = 0.65, m_cold = 0.71895, p_hot_in = 100, p_cold_in = 300, 
                                  hot_fluid = f'HEOS::CO2[{z_CO2}]&O2[{z_O2}]&N2[{z_N2}]&H2O[{z_H2O}]', cold_fluid = 'INCOMP::T66', 
                                  flow_arrangement = 'counter')
        
        test_HE2.calc_temp_dist(points= 10)
        test_HE2_fig, test_HE2_ax = test_HE2.plot_temp_dist_df()
        test_HE2_ax.set_title('Test HE2')

        self.assertEqual(round(test_HE2.dT_min,2), 40.0)
        
        """
        Test case 3 - heat transfer with phase change on one of the sides eg. evaporation
        """
        
        test_HE3 = HeatExchanger(T_hot_in = 310, T_hot_out = 182.36, T_cold_in = 86.25, T_cold_out = 209.37, 
                                  m_hot = 0.71874, m_cold = 0.395, p_hot_in = 300, p_cold_in = 810, 
                                  hot_fluid = 'INCOMP::T66', cold_fluid = 'Toluene', 
                                  flow_arrangement = 'counter')
        
        test_HE3.calc_temp_dist(points= 4)
        test_HE3_fig, test_HE3_ax = test_HE3.plot_temp_dist_df()
        test_HE3_ax.set_title('Test HE3')
        
        self.assertEqual(round(test_HE3.dT_min,2), 40.0)
        
        """
        Test case 4- heat transfer with phase change on one of the sides eg. condensation.
        """
        
        test_HE4 = HeatExchanger(T_hot_in = 107.95, T_hot_out = 56.92, T_cold_in = 25.00, T_cold_out = 56.10, 
                                  m_hot = 0.395, m_cold = 1.633, p_hot_in = 20, p_cold_in = 300, 
                                  hot_fluid = 'Toluene', cold_fluid = 'INCOMP::MEG-40%', 
                                  flow_arrangement = 'counter')
        
        test_HE4.calc_temp_dist(points= 4)
        test_HE4_fig, test_HE4_ax = test_HE4.plot_temp_dist_df()
        test_HE4_ax.set_title('Test HE4')
        
        self.assertEqual(round(test_HE4.dT_min,2), 10.0)
        
        
        """
        Test case 4 - heat transfer with phase change on both sides.
        """

if __name__=='__main__':
	unittest.main()
    

