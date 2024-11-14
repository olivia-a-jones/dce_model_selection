# -*- coding: utf-8 -*-
"""
Simulate and Model Select Script
This script generates synthetic DCE-MRI data, applies noise, and fits various
pharmacokinetic models to estimate Ktrans.

Make sure to check the parameters for the tissue properties and acquisition parameters. 

@author: Olivia A Jones
"""

import numpy as np
from numpy import loadtxt
from QbiMadym import madym_DCE_lite
from QbiPy.dce_models import tissue_concentration
from QbiMadym.madym_DWI_lite import add_rician_noise
import os
import time

class SimulateAndModelSelect:
    def __init__(self, Fp, PS, ve, vp, temp_res, n_times, TR, FA, T1_assumption, r1_const, dose, hct, aif_name, max_sigma):
        """
        Initialize the SimulatorNoiseKtrans class with user-defined parameters.

        Parameters:
        Fp (np.ndarray): Blood flow per unit volume of tissue (mL/min/mL)
        PS (np.ndarray): Permeability-surface area product (mL/min/mL)
        ve (np.ndarray): Extravascular extracellular volume fraction
        vp (np.ndarray): Vascular volume fraction
        temp_res (float): Temporal resolution of the DCE-MRI acquisition (minutes)
        TR (float): Repetition time of the DCE-MRI acquisition (milliseconds)
        FA (float): Flip angle of the DCE-MRI acquisition (degrees)
        T1_assumption (float): Assumed T1 relaxation time of the tissue (milliseconds)
        r1_const (float): Relaxivity of the contrast agent (mM^-1 s^-1)
        dose (float): Contrast agent dose (mmol/kg)
        hct (float): Hematocrit fraction
        aif_name (str): Path to the arterial input function (AIF) data file
        """
        self.Fp = Fp
        self.PS = PS
        self.ve = ve
        self.vp = vp
        self.temp_res = temp_res
        self.TR = TR
        self.FA = FA
        self.T1_assumption = T1_assumption
        self.r1_const = r1_const
        self.dose = dose
        self.hct = hct
        self.aif_name = aif_name

        self.n_samples = len(self.PS)
        self.n_times = n_times
        self.t = np.arange(self.n_times) * self.temp_res
        self.noise_list = np.linspace(0, max_sigma, self.n_samples)
        
        self.load_aif()

    def load_aif(self):
        """Load the arterial input function (AIF) from the specified file."""
        self.xdata, self.ydata = loadtxt(self.aif_name, dtype='float', unpack=True)
    
    def generate_tissue_signals(self):
        """
        Generate synthetic tissue concentration and signal time-series for each combination
        of PS value and noise level.
    
        Returns:
        np.ndarray: Noisy tissue signal time-series (noise_levels x samples x time_points)
        np.ndarray: Clean tissue signal time-series (samples x time_points)
        """
        # Create dummy input time-series
        Ct_input = np.zeros((self.n_samples, self.n_times))
    
        # Concatenate the model parameters into an initial parameters array
        init_params = np.concatenate(
            (self.Fp.reshape(-1, 1),
             self.PS.reshape(-1, 1),
             self.ve.reshape(-1, 1),
             self.vp.reshape(-1, 1),
             np.zeros(self.n_samples).reshape(-1, 1)), axis=1
        )
    
        # Generate the clean tissue concentration time-series
        C_t = madym_DCE_lite.run(
            model='2CXM',
            input_data=Ct_input,
            dyn_times=self.t,
            no_optimise=True,
            init_params=init_params,
            injection_image=8,
            dose=self.dose,
            hct=self.hct,
            aif_name=self.aif_name
        )[4]
    
        # Convert the tissue concentration to signal
        S_t = tissue_concentration.concentration_to_signal(
            C_t=C_t,
            T1_0=np.ones(self.n_samples) * self.T1_assumption,
            M0=np.ones(self.n_samples) * 100,
            FA=self.FA,
            TR=self.TR,
            relax_coeff=self.r1_const,
            use_M0_ratio=8
        )
    
        # Create a 3D array to store noisy signals for each noise level and PS value
        S_tn = np.zeros((len(self.noise_list), self.n_samples, self.n_times))
        
        # Apply different noise levels to each PS value's signal
        for i, noise_sigma in enumerate(self.noise_list):
            S_tn[i, :, :] = add_rician_noise(S_t, noise_sigma)
    
        return S_tn, S_t

    def fit_pharmacokinetic_models(self, S_tn):
        """
        Fit various pharmacokinetic models to the noisy tissue signals.

        Parameters:
        S_tn (np.ndarray): Noisy tissue signal time-series

        Returns:
        dict: Fitting results for each model
        """
        model_names = ['NLM', 'PTM', 'ETM']
        fit_results_allsigma = []

        for sigma in range(self.n_samples):
            fit_results = {}
            for model_name in model_names:
                if model_name == 'NLM':
                    model = 'PATLAK'
                    fixed_params = [1]
                    fixed_values = [0]
                    upper_bounds = np.ones((self.n_samples, 3)) * (1, 1, 20)
                    lower_bounds = np.ones((self.n_samples, 3)) * (0, 0, -20)
                    init_params = np.ones((self.n_samples, 3)) * (0, 2.000000e-02, 0)
                elif model_name == 'PTM':
                    model = 'PATLAK'
                    fixed_params = None
                    fixed_values = None
                    upper_bounds = np.ones((self.n_samples, 3)) * (1, 1, 20)
                    lower_bounds = np.ones((self.n_samples, 3)) * (-1.000000e-03, 0, -20)
                    init_params = np.ones((self.n_samples, 3)) * (1.000000e-04, 2.000000e-02, 0)
                else:
                    model = 'ETM'
                    fixed_params = None
                    fixed_values = None
                    upper_bounds = np.ones((self.n_samples, 4)) * (1, 1, 1, 20)
                    lower_bounds = np.ones((self.n_samples, 4)) * (-1.000000e-03, 0, 0, -20)
                    init_params = np.ones((self.n_samples, 4)) * (1.000000e-04, 1.000000e-01, 2.000000e-02, 0)

                fit_results[model_name] = madym_DCE_lite.run(
                    model=model,
                    input_data=S_tn[sigma, :, :],
                    dyn_times=self.t,
                    input_Ct=False,
                    T1=np.ones(self.n_samples) * self.T1_assumption,
                    FA=self.FA,
                    r1_const=self.r1_const,
                    TR=self.TR,
                    aif_name=self.aif_name,
                    dose=self.dose,
                    injection_image=8,
                    hct=self.hct,
                    fixed_params=fixed_params,
                    fixed_values=fixed_values,
                    init_params=init_params,
                    max_iter=1000,
                    opt_type='BLEIC'
                )
            fit_results_allsigma.append(fit_results)

        return fit_results_allsigma

    def evaluate_model_performance(self, fit_results_allsigma):
        """
        Evaluate the performance of the fitted pharmacokinetic models.

        Parameters:
        fit_results_allsigma (list): Fitting results for each noise level

        Returns:
        np.ndarray: Best fitting model index for each noise level and parameter set
        np.ndarray: Ktrans values for the best fitting model for each noise level and parameter set
        np.ndarray: Ktrans values for the Patlak model for each noise level and parameter set
        np.ndarray: Ktrans values for the Extended Tofts model for each noise level and parameter set
        """
        best_fits = np.zeros((self.n_samples, self.n_samples))
        best_Ktrans = np.zeros((self.n_samples, self.n_samples))
        PTM_Ktrans = np.zeros((self.n_samples, self.n_samples))
        ETM_Ktrans = np.zeros((self.n_samples, self.n_samples))

        for noises in range(self.n_samples):
            for params in range(self.n_samples):
                Akaike_corr = {}
                Ktrans_fit = {}
                for model_name in ['NLM', 'PTM', 'ETM']:
                    RSS = fit_results_allsigma[noises][model_name][1][params]
                    model_params = fit_results_allsigma[noises][model_name][0][params]
                    N = self.n_times
                    Ktrans_fit[model_name] = model_params.T[0]
                    if model_name == 'NLM':
                        K = 1
                    elif model_name == 'PTM':
                        K = 2
                    elif model_name == 'ETM':
                        K = 3
                    Akaike = (2 * (K + 1)) + (N * np.log(RSS))
                    Akaike_corr[model_name] = Akaike + ((2 * K) * (K + 1)) / (N - K - 1)

                delta_NLM = Akaike_corr['NLM'] - min(Akaike_corr['PTM'], Akaike_corr['ETM'])
                delta_PTM = Akaike_corr['PTM'] - min(Akaike_corr['NLM'], Akaike_corr['ETM'])
                delta_ETM = Akaike_corr['ETM'] - min(Akaike_corr['NLM'], Akaike_corr['PTM'])
                edelta_NLM = np.exp(-delta_NLM / 2)
                edelta_PTM = np.exp(-delta_PTM / 2)
                edelta_ETM = np.exp(-delta_ETM / 2)
                Ak_weights = [edelta_NLM / (edelta_PTM + edelta_ETM), edelta_PTM / (edelta_NLM + edelta_ETM), edelta_ETM / (edelta_NLM + edelta_PTM)]
                best_fits[noises, params] = Ak_weights.index(max(Ak_weights))
                if best_fits[noises, params] == 0:
                    best_Ktrans[noises, params] = Ktrans_fit['NLM']
                elif best_fits[noises, params] == 1:
                    best_Ktrans[noises, params] = Ktrans_fit['PTM']
                elif best_fits[noises, params] == 2:
                    best_Ktrans[noises, params] = Ktrans_fit['ETM']
                PTM_Ktrans[noises, params] = Ktrans_fit['PTM']
                ETM_Ktrans[noises, params] = Ktrans_fit['ETM']

        return best_fits, best_Ktrans, PTM_Ktrans, ETM_Ktrans

    def save_results(self, repeat, best_fits, best_Ktrans, PTM_Ktrans, ETM_Ktrans):
        """
        Save the simulation results to text files.

        Parameters:
        repeat (int): The current repeat number
        best_fits (np.ndarray): Best fitting model index for each noise level and parameter set
        best_Ktrans (np.ndarray): Ktrans values for the best fitting model for each noise level and parameter set
        PTM_Ktrans (np.ndarray): Ktrans values for the Patlak model for each noise level and parameter set
        ETM_Ktrans (np.ndarray): Ktrans values for the Extended Tofts model for each noise level and parameter set
        """
        os.mkdir('output')
        np.savetxt('output/ground_truth_PS.txt', self.PS, fmt='%f')
        np.savetxt('output/applied_sigma.txt', self.noise_list, fmt='%f')
        np.savetxt(f'output/BFModel_Repeat_{repeat}.txt', best_fits + 1, fmt='%f')
        np.savetxt(f'output/BFKtrans_Repeat_{repeat}.txt', best_Ktrans, fmt='%f')
        np.savetxt(f'output/PTMKtrans_Repeat_{repeat}.txt', PTM_Ktrans, fmt='%f')
        np.savetxt(f'output/ETMKtrans_Repeat_{repeat}.txt', ETM_Ktrans, fmt='%f')

    def run_simulation(self, n_repeats):
        """
        Run the simulation for the specified number of repeats.

        Parameters:
        n_repeats (int): The number of simulation repeats to perform
        """
        start = time.time()

        for repeat in range(n_repeats):
            print(f'Processing repeat {repeat} ...')

            # Generate noisy and clean tissue signals
            S_tn, S_t = self.generate_tissue_signals()

            # Fit the pharmacokinetic models
            fit_results_allsigma = self.fit_pharmacokinetic_models(S_tn)

            # Evaluate the model performance
            best_fits, best_Ktrans, PTM_Ktrans, ETM_Ktrans = self.evaluate_model_performance(fit_results_allsigma)

            # Save the results
            self.save_results(repeat, best_fits, best_Ktrans, PTM_Ktrans, ETM_Ktrans)

        end = time.time()
        print(f'Took {end - start} seconds for {n_repeats} repeats.')
        
# Define simulated tissue parameters, e.g. for NAWM
PS = np.concatenate((np.linspace(1e-4,1e-3,9,endpoint=False), np.linspace(1e-3,1e-2, 10))) # PS samples in min-1
n_samples = len(PS)
Fp = np.ones(n_samples) * 0.4              # F_p in min-1
ve = np.ones(n_samples) * 0.2              # v_e
vp = np.ones(n_samples) * 0.01             # v_p
T1_assumption = 1000                # assumed T1 in ms

# Define scan parameters
temp_res = 7.64 / 60                # temporal resolution in min
n_times = 160                       # No of dynamics
TR = 2.4                            # repetition time in ms
FA = 10                             # flip angle in degrees

# Define contrast agent parameters
r1_const = 3.4                      # r1 constant in s-1mM-1
dose = 0.1                          # dose in mM/kg

# Define VIF parameters
hct = 0.42                          # Hematocrit
aif_name = 'group_average_VIF.txt'  # VIF filename

# Define maximum sigma of Rician noise
max_sigma = 6

# Create the Sim23NoiseKtrans object
sim = SimulateAndModelSelect(Fp, PS, ve, vp, temp_res, n_times, TR, FA, T1_assumption, r1_const, dose, hct, aif_name, max_sigma)

# Run the simulation for 250 repeats
sim.run_simulation(n_repeats=1000)
