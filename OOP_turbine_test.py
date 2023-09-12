import unittest
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import CoolProp.CoolProp as CP

from OOP_turbine import Turbine
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.lines import Line2D

class TestOOPTurbine(unittest.TestCase):
      
    def test_Ts_diagram(self):
        """
        Drawing test T-s diagrams for test cases.
        """
        
        print('Runing test T-s diagrams')

        test_turbine1 = Turbine(155, 3620.0, 1568.5, 237.71, 'R125')
        test_turbine1.efficiency_estimation(n_stages= 1)
        test_turbine1.Ts_diagram()
        
        test_turbine2 = Turbine(155, 3620.0, 1568.5, 11.89, 'R125')
        test_turbine2.efficiency_estimation(n_stages = 1)
        test_turbine2.Ts_diagram()
        
        test_turbine3 = Turbine(234.67, 829.0,25.0, 40.78, 'Hexane')
        test_turbine3.efficiency_estimation(n_stages = 1)
        test_turbine3.Ts_diagram()
        
        test_turbine4 = Turbine(234.67, 829.0,25.0, 2.04, 'Hexane')
        test_turbine4.efficiency_estimation(n_stages = 1)
        test_turbine4.Ts_diagram()
         
    def test_estimated_efficiency(self):
        
        """
        This test is based on data from:
            
        Macchi E, Astolfi M. Organic rankine cycle (ORC) power systems: Technologies and
        applications. Duxford, United Kingdom: Woodhead Publishing is an imprint of Elsevier;
        2017.
        """
        
        print('Runing test efficiency estimations')

        test_turbine1 = Turbine(155, 3620.0, 1568.5, 237.71, 'R125')
        test_turbine1.efficiency_estimation(n_stages= 1)
        test_eta_iT1 = test_turbine1.eta_iT
        
        test_turbine2 = Turbine(155, 3620.0, 1568.5, 11.89, 'R125')
        test_turbine2.efficiency_estimation(n_stages = 1)
        test_eta_iT2 = test_turbine2.eta_iT
    
        test_turbine3 = Turbine(234.67, 829.0,25.0, 40.78, 'Hexane')
        test_turbine3.efficiency_estimation(n_stages = 1)
        test_eta_iT3 = test_turbine3.eta_iT
    
        test_turbine4 = Turbine(234.67, 829.0,25.0, 2.04, 'Hexane')
        test_turbine4.efficiency_estimation(n_stages = 1)
        test_eta_iT4 = test_turbine4.eta_iT
        
        self.assertEqual(test_eta_iT1, 0.9028177495017103)
        self.assertEqual(test_eta_iT2, 0.8770645657764151)
        self.assertEqual(test_eta_iT3, 0.835051048808327)
        self.assertEqual(test_eta_iT4, 0.8065933042772973)

    def test_case_off_design1(self):
        
        """
        Testing of off-design calculations by ploting calculated off-design calculations.
        
        Test assumes constant pressure at turbine outlet, temperature at turbine inlet 
        """
        
        print('Running off-design test 1')
        
        'Determining design parameters of turbine'
        T_in = 155.0
        p_in = 3620.0
        p_out = 1568.5
        mf = 237.71
        medium = 'R125'
        nos = 1
        
        'Calculating designed parameters of turbine'
        test_turbine1 = Turbine(T_in, p_in, p_out, mf, medium)
        test_turbine1.efficiency_estimation(n_stages = nos)
        test_eta_iT1 = test_turbine1.eta_iT
        test_N_iT1 = test_turbine1.N_iT
        test_T_out1 = test_turbine1.T_out
        
        'Creating sample array of parameters for off design calculations'
        num = 100
        
        'Creating min and max values of pressure for creating off-design dataset'
        
        p_in_off_ar = np.linspace(1568.5, 3617.7, num=num, endpoint = True)
        
        delta_T_sup = 80
        
        T_in_off_min = CP.PropsSI('T', 'P', min(p_in_off_ar)*1000, 'Q', 1, 'R125')-273.15 + delta_T_sup
        
        T_in_off_ar = np.linspace(T_in_off_min, T_in, num, endpoint = True) #CP.PropsSI('T', 'P', p_in_off_ar*1000, 'Q', np.full(num,1), 'R125')-273.15 + delta_T_sup
        
        'Creating list to collect results of off-design calculations'
        eta_iT_off_ar = []
        N_iT_off_ar = []
        T_out_off_ar = []
        mf_off_ar = []
        
        for T_in_off, p_in_off in zip(T_in_off_ar,p_in_off_ar):
            
            test_turbine1.off_design(T_in_off, p_in_off, p_out)
            eta_iT_off_ar.append(test_turbine1.eta_iT_off)
            N_iT_off_ar.append(test_turbine1.N_iT_off)
            T_out_off_ar.append(test_turbine1.T_out_off)
            mf_off_ar.append(test_turbine1.mf_off)
        
        'Ploting set of diagrams to visualize '
        fig1, axs = plt.subplots(nrows = 4, ncols = 1, dpi = 400, figsize = (5,10), sharex = True)
        
        ((ax1,ax2, ax3,ax4)) = axs
        
        'Ploting off-design efficiency'
        ax1.plot(p_in_off_ar, eta_iT_off_ar, label = 'eta')
        ax1.scatter(p_in, test_eta_iT1, marker ='x', c = 'red')
        ax1.set_ylabel('Efficiency $\eta_{iT}$, -')
        ax1.set_ylim(0.0,1.0)
        ax1.set_yticks(np.arange(0.0,1.1,0.1))
        
        'Ploting off-design power'
        ax2.plot(p_in_off_ar, N_iT_off_ar, label = 'N')
        ax2.set_ylabel('Power $N_{iT}$, kW')
        ax2.scatter(p_in, test_N_iT1, marker ='x', c = 'red')
        
        'Ploting off-design temperature at turbine outlet'
        ax3.plot(p_in_off_ar, T_out_off_ar, label = 'T_out')
        ax3.set_ylabel(r'Outlet temperature $T_{out}$, $^\circ$C')
        ax3.scatter(p_in, test_T_out1, marker ='x', c = 'red')
        
        'Ploting off-design mass flow'
        ax4.plot(p_in_off_ar, mf_off_ar, label = 'mf_off')
        ax4.set_ylabel('Mass flow, kg/s')
        ax4.set_xlabel(r'Inlet pressure $p_{in}$, kPa(a)')
        ax4.scatter(p_in, mf, marker ='x', c = 'red')
        ax4.set_xticks(np.arange(1500,3800,100))
        ax4.set_xticklabels(ax4.get_xticks(),rotation = 90)
        
        'List of minimal values of power, temperature at turbine outlet and mass flow'
        min_vals = [min(N_iT_off_ar), min(T_out_off_ar), min(mf_off_ar)]
        
        'Corection of yticks'
        for ax, min_val in zip(axs[1:], min_vals):
            yticks = ax.get_yticks()
            corected_yticks = [y for y in yticks if y>=min_val]
            ax.set_yticks(corected_yticks)
            
        plt.suptitle('Turbine off-design characteritics')
        
        for ax in axs:
            ax.grid(axis ='x', linestyle = '--')
            ax.grid(axis ='y', linestyle = ':')
            
        plt.tight_layout()

    def test_case_off_design2(self):
        
        """
        Testing of off-design calculations by ploting calculated off-design calculations.
        
        Test assumes constant pressure at turbine outlet, temperature at turbine inlet
        
        Useful links:
            * https://stackoverflow.com/questions/23876588/matplotlib-colorbar-in-each-subplot
        """
        
        print('Running off-design test 2')
        
        'Determining design parameters of turbine'
        T_in = 155.0
        p_in = 3620.0
        p_out = 1568.5
        mf = 237.71
        medium = 'R125'
        nos = 1
        
        'Calculating designed parameters of turbine'
        test_turbine1 = Turbine(T_in, p_in, p_out, mf, medium)
        test_turbine1.efficiency_estimation(n_stages = nos)
        test_eta_iT1 = test_turbine1.eta_iT
        test_N_iT1 = test_turbine1.N_iT
        test_T_out1 = test_turbine1.T_out
        
        'Creating sample array of parameters for off design calculations'
        num = 100
        
        'Calculating minimal and maximal pressure for off-design calculations'
        p_min = 1500
        p_max = 1.1*p_in # It is assumed that parameters of turbine can be exceeded by aprox. 10% - this is simplification
        
        'Aray of off design inlet and outlet pressure'
        p_in_off_ar  = np.linspace(p_min, p_max, num=num, endpoint = True) #np.linspace(1568.5, 3617.7, num=num, endpoint = True)
        p_out_off_ar = p_in_off_ar #np.linspace(1568.5, 3617.7, num=num, endpoint = True)
        
        'Superheating of vapour'
        delta_T_sup = 80
        
        'Calculating minimal and maximal temperature'
        T_in_off_min = CP.PropsSI('T', 'P', p_out*1000, 'Q', 1, 'R125')-273.15 + delta_T_sup
        T_in_off_max = 1.1*T_in     # It is assumed that parameters of turbine can be exceeded by aprox. 10% - this is simplification
        
        'Aray of off design inlet temperature'
        T_in_off_ar = np.linspace(T_in_off_min, T_in, num=num, endpoint = True) #CP.PropsSI('T', 'P', p_in_off_ar*1000, 'Q', np.full(num,1), 'R125')-273.15 + delta_T_sup
        
        'Combining arays to create frame of changing prameters'
        p_in_off_ar = np.tile(p_in_off_ar, num)
        T_in_off_ar =  np.tile(T_in_off_ar, num)
        p_out_off_ar = np.repeat(p_out_off_ar, num)
        
        'Creating dict of parameters for further data frame'
        off_design_dict = {'T_in_off' : T_in_off_ar,
                           'p_in_off' : p_in_off_ar,
                           'p_out_off' : p_out_off_ar}

        off_design_frame = pd.DataFrame(off_design_dict)
        
        'Ensuring that pressure at inlet always exceeds pressure at inlet'
        off_design_frame = off_design_frame[off_design_frame['p_out_off']<=off_design_frame['p_in_off']]
        
        'Creating list to collect results of off-design calculations'
        eta_iT_off_ar = []
        N_iT_off_ar = []
        T_out_off_ar = []
        mf_off_ar = []
        
        for T_in_off, p_in_off, p_out_off in zip(off_design_frame.T_in_off, off_design_frame.p_in_off, off_design_frame.p_out_off):
            
            test_turbine1.off_design(T_in_off, p_in_off, p_out_off)
            eta_iT_off_ar.append(test_turbine1.eta_iT_off)
            N_iT_off_ar.append(test_turbine1.N_iT_off)
            T_out_off_ar.append(test_turbine1.T_out_off)
            mf_off_ar.append(test_turbine1.mf_off)
        
        'Ploting set of diagrams to visualize '
        fig1, axs = plt.subplots(nrows = 4, ncols = 1, dpi = 400, figsize = (5,10), sharex = True, layout='constrained')
        
        ((ax1,ax2, ax3,ax4)) = axs
        
        'Ploting off-design efficiency'
        heatmap1 = ax1.scatter(off_design_frame.p_in_off, off_design_frame.p_out_off, c = eta_iT_off_ar, cmap = 'jet', marker = '.', label = 'eta')
        
        'Ploting off-design power'
        heatmap2 = ax2.scatter(off_design_frame.p_in_off, off_design_frame.p_out_off, c = N_iT_off_ar, cmap = 'jet', marker = '.', label = 'eta')
        
        'Ploting off-design temperature at turbine outlet'
        heatmap3 = ax3.scatter(off_design_frame.p_in_off, off_design_frame.p_out_off, c = T_out_off_ar, cmap = 'jet', marker = '.', label = 'eta')
        
        'Ploting off-design mass flow'
        heatmap4 = ax4.scatter(off_design_frame.p_in_off, off_design_frame.p_out_off, c = mf_off_ar, cmap = 'jet', marker = '.', label = 'eta')
        
        'Customizing shared x-axis'
        xticks = [p_min]+list(ax4.get_xticks())+[p_max]
        xticks.sort()
        
        ax4.set_xticks(xticks)
        ax4.set_xlim(p_min, p_max)
        ax4.set_xlabel(r'Inlet pressure $p_{in}$, kPa(a)')
        #ax4.set_xticklabels(labels = [1000.0, 1500.0, 1568.5, 2000.0, 2500.0, 3000.0, 3500.0, '\n3620.0', 4000.0])

        plt.suptitle('Turbine off-design characteritics')
        
        """
        Common customization for all subplots:
            * formating y-label
            * formating x lims
            * adding colorbars with labels to subplots
        """
        
        'List of subplots to iterate'
        heatmaps = [heatmap1, heatmap2, heatmap3, heatmap4]
        
        'Names of colorbars labels'
        labels = ['Efficiency $\eta_{iT}$, -', 'Power $N_{iT}$, kW', 'Outlet temperature $T_{out}$, $^\circ$C', 'Mass flow, kg/s']
        
        'List of minimal values of colorbar ticks'
        cbar_minticks = np.array([min(eta_iT_off_ar), float(min(N_iT_off_ar)), min(T_out_off_ar), min(mf_off_ar)])
        cbar_minticks.flatten()
        
        'List of maximal values of colorbar ticks'
        cbar_maxticks = np.array([max(eta_iT_off_ar), float(max(N_iT_off_ar)), max(T_out_off_ar), max(mf_off_ar)])
        cbar_maxticks.flatten()
        
        'Creating lengend element for plots'
        legend_element = [Line2D([0],[0], c = 'black', linestyle = '--', linewidth = 0.5)]
        
        for ax, hmap, label, mintick, maxtick in zip(axs, heatmaps, labels, cbar_minticks, cbar_maxticks):
            ax.set_yticks(list(ax.get_xticks())+[p_min, p_max])
            ax.set_ylim(p_min, p_max)
            ax.set_ylabel('Outlet pressure $p_{out}$, \n kPa(a)')
            ax.grid(axis ='both', linestyle = '--')
            ax.axhline(p_out, c = 'black', linestyle = '--', linewidth = 0.5)
            ax.axvline(p_in, c = 'black', linestyle = '--', linewidth = 0.5)
            ax.set_axisbelow(True)
            
            fig1.colorbar(hmap, ticks = np.linspace(mintick, maxtick, 10, True), ax = ax, label = label)
        
        fig1.legend(loc = 'outside lower center', handles = legend_element, labels = ['Design point at crossing lines',])
        
            
if __name__=='__main__':
 	
    unittest.main()