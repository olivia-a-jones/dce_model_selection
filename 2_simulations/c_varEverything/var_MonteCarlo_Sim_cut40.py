#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 19 12:38:51 2025

@author: p15094oj

Monte Carlo DCE simulation with random parameter sampling
"""

import numpy as np
from QbiMadym import madym_DCE_lite
from QbiPy.dce_models import tissue_concentration
from QbiMadym.madym_DWI_lite import add_rician_noise
import time
import concurrent.futures
import os

def log_uniform_sample(low, high, size):
    """Sample from log-uniform distribution"""
    return np.exp(np.random.uniform(np.log(low), np.log(high), size))

def uniform_sample(low, high, size):
    """Sample from uniform distribution"""
    return np.random.uniform(low, high, size)

def sample_ground_truth_params(n_samples):
    """
    Sample ground truth parameters using the ranges specified:
    - T1_0: 800-2000 ms (uniform)
    - PS: 10^-4 to 10^-2 (log-uniform)
    - v_p: 1-10% (uniform, converted to fraction)
    - v_e: 5-25% (uniform, converted to fraction) 
    - F_p: 0.1-0.8 (uniform)
    """
    
    # Sample parameters
    T1_0 = uniform_sample(800, 2000, n_samples)
    PS = log_uniform_sample(1e-4, 1e-2, n_samples)
    v_p = uniform_sample(0.01, 0.10, n_samples)  # 1-10% as fractions
    v_e = uniform_sample(0.05, 0.25, n_samples)  # 5-25% as fractions
    F_p = uniform_sample(0.1, 0.8, n_samples)
    tau_a = np.zeros(n_samples)  # Keep tau_a as zero
    
    # Stack parameters for madym input (F_p, PS, v_e, v_p, tau_a)
    init_params = np.column_stack((F_p, PS, v_e, v_p, tau_a))
    
    return T1_0, init_params

def process_model_fit(args):
    S_tn_sample, model_name, t, T1_0, FA, r1_const, TR, n_samples, n_times = args
    
    # Pre-define model configurations
    model_configs = {
        'NLM': {
            'model': 'PATLAK',
            'fixed_params_uncut': [1],
            'fixed_values_uncut': [0],
            'fixed_params_cut': [1, 3],
            'fixed_values_cut': [0, 0],
            'init_params': np.ones((n_samples, 3)) * [0, 2.000000e-02, 0],
            'bounds': (
                np.ones((n_samples, 3)) * [0, 0, -20],
                np.ones((n_samples, 3)) * [1, 0.55, 20] # vp limit is 1 - Hct
            )
        },
        'PTM': {
            'model': 'PATLAK',
            'fixed_params_uncut': None,
            'fixed_values_uncut': None,
            'fixed_params_cut': [3],
            'fixed_values_cut': [0],
            'init_params': np.ones((n_samples, 3)) * [1.000000e-03, 2.000000e-02, 0],
            'bounds': (
                np.ones((n_samples, 3)) * [-1, 0, -20],
                np.ones((n_samples, 3)) * [1, 0.55, 20] # vp limit is 1 - Hct
            )
        },
        'ETM': {
            'model': 'ETM',
            'fixed_params_uncut': None,
            'fixed_values_uncut': None,
            'fixed_params_cut': [4],
            'fixed_values_cut': [0],
            'init_params': np.ones((n_samples, 4)) * [1.000000e-03, 2.000000e-01, 2.000000e-02, 0], # v_e init 20%, v_p init 2%
            'bounds': (
                np.ones((n_samples, 4)) * [-1, 0, 0, -20],
                np.ones((n_samples, 4)) * [1, 0.55, 0.55, 20] # vp limit is 1 - Hct
            )
        }
    }
    
    config = model_configs[model_name]
    
    # First run
    results_uncut = madym_DCE_lite.run(
        model=config['model'],
        input_data=S_tn_sample,
        dyn_times=t,
        input_Ct=False,
        T1=T1_0,
        FA=FA,
        r1_const=r1_const,
        TR=TR,
        aif_name='/Users/user/Documents/Python_scripts/MODELLING_SIMULATIONS/controls_aif_avg.txt',
        dose=0.1,
        injection_image=8,
        hct=0.42,
        fixed_params=config['fixed_params_uncut'],
        fixed_values=config['fixed_values_uncut'],
        init_params=config['init_params'],
        max_iter=1000,
        opt_type='BLEIC'
    )
    
    # Second run with cut data
    results_cut = madym_DCE_lite.run(
        model=config['model'],
        input_data=results_uncut[5],
        dyn_times=t,
        input_Ct=True,
        T1=T1_0,
        FA=FA,
        r1_const=r1_const,
        TR=TR,
        aif_name='/Users/user/Documents/Python_scripts/MODELLING_SIMULATIONS/controls_aif_avg.txt',
        dose=0.1,
        injection_image=8,
        hct=0.42,
        fixed_params=config['fixed_params_cut'],
        fixed_values=config['fixed_values_cut'],
        init_params=config['init_params'],
        max_iter=1000,
        opt_type='BLEIC',
        first_image=14
    )
    
    return results_cut

def main():
    # Monte Carlo simulation parameters
    n_samples = 10000  # Number of random parameter combinations per repeat
    n_repeats = 100   # Number of Monte Carlo repeats
    
    # Fixed simulation parameters
    tissue = 'mc'  # Monte Carlo - mixed tissue types due to random T1
    
    # Time points
    n_times = 160
    temp_res_mins = 7.64/60
    t = np.arange(n_times) * temp_res_mins
    
    # Signal parameters
    r1_const = 3.4
    FA = 10
    TR = 2.4
    S_t0 = 100
    
    model_names = ['NLM', 'PTM', 'ETM']
    
    # Create output directory
    base_path = f'/Users/user/Documents/Python_scripts/MODELLING_SIMULATIONS/Paper_Draft_2/MonteCarlo_cut40_{tissue}/'
    if not os.path.isdir(base_path):
        os.makedirs(base_path)
    
    for repeat in range(n_repeats):
        print(f'Starting Monte Carlo repeat {repeat + 1}/{n_repeats}...')
        
        # Sample ground truth parameters
        T1_0_samples, init_params = sample_ground_truth_params(n_samples)
        
        # Extract individual parameters for saving
        F_p_samples = init_params[:, 0]
        PS_samples = init_params[:, 1] 
        v_e_samples = init_params[:, 2]
        v_p_samples = init_params[:, 3]
        tau_a_samples = init_params[:, 4]
        
        # Generate concentration time series for all samples
        Ct_input = np.zeros((n_samples, n_times))
        
        C_t = madym_DCE_lite.run(
            model='2CXM',
            input_data=Ct_input,
            dyn_times=t,
            no_optimise=True,
            init_params=init_params,
            injection_image=8,
            dose=0.1,
            hct=0.42,
            aif_name='/Users/user/Documents/Python_scripts/MODELLING_SIMULATIONS/controls_aif_avg.txt',
        )[4]
        
        # Convert to signal using sampled T1 values
        S_t0_all = np.ones(n_samples) * S_t0
        S_t = tissue_concentration.concentration_to_signal(
            C_t=C_t,
            T1_0=T1_0_samples,
            M0=S_t0_all,
            FA=FA,
            TR=TR,
            relax_coeff=r1_const,
            use_M0_ratio=8
        )
        
        # Generate noisy signals
        S_tn = add_rician_noise(S_t, 3.0)  # Fixed noise level
        
        # Process models in parallel
        fit_results_allsigma = []
        
        # Prepare arguments for parallel processing
        args_list = [(S_tn, model_name, t, T1_0_samples, FA, r1_const, TR, n_samples, n_times) 
                    for model_name in model_names]
        
        # Process models in parallel
        with concurrent.futures.ProcessPoolExecutor() as executor:
            results = list(executor.map(process_model_fit, args_list))
        
        # Store results in dictionary format
        fit_results = {model_name: result for model_name, result in zip(model_names, results)}
        
        # Calculate metrics and determine best fit models
        best_fits = np.zeros(n_samples)
        best_Ktrans = np.zeros(n_samples)
        best_vp = np.zeros(n_samples)
        best_ve = np.zeros(n_samples)
        
        # Store all fitted parameters for each model
        NLM_Ktrans = np.zeros(n_samples)
        PTM_Ktrans = np.zeros(n_samples) 
        ETM_Ktrans = np.zeros(n_samples)
        NLM_vp = np.zeros(n_samples)
        PTM_vp = np.zeros(n_samples)
        ETM_vp = np.zeros(n_samples)
        ETM_ve = np.zeros(n_samples)
        
        for sample_idx in range(n_samples):
            Akaike_corr = {}
            Ktrans_fit = {}
            vp_fit = {}
            ve_fit = {}
            
            for model_name in model_names:
                results = fit_results[model_name]
                RSS = results[1][sample_idx]
                model_params = results[0][sample_idx]
                N = n_times - 15  # Number of data points after cutting
                
                Ktrans_fit[model_name] = model_params[0]
                
                # Extract vp and ve parameters based on model
                if model_name == 'PTM':
                    vp_fit[model_name] = model_params[1]  # PTM fits vp (Patlak: Ktrans, vp)
                elif model_name == 'ETM':
                    ve_fit[model_name] = model_params[1]
                    vp_fit[model_name] = model_params[2]
                else:  # NLM
                    vp_fit[model_name] = model_params[1]  # NLM fits vp (Patlak: Ktrans, vp)
                
                # Store fitted parameters for each model
                if model_name == 'NLM':
                    NLM_Ktrans[sample_idx] = Ktrans_fit[model_name]
                    NLM_vp[sample_idx] = vp_fit[model_name]
                elif model_name == 'PTM':
                    PTM_Ktrans[sample_idx] = Ktrans_fit[model_name] 
                    PTM_vp[sample_idx] = vp_fit[model_name]
                elif model_name == 'ETM':
                    ETM_Ktrans[sample_idx] = Ktrans_fit[model_name]
                    ETM_vp[sample_idx] = vp_fit[model_name]
                    ETM_ve[sample_idx] = ve_fit[model_name]
                
                # Calculate AIC
                K = {'NLM': 1, 'PTM': 2, 'ETM': 3}[model_name]
                Akaike = (2 * (K + 1)) + (N * np.log(RSS))
                Akaike_corr[model_name] = Akaike + ((2 * K) * (K + 1)) / (N - K - 1)
            
            # Calculate Akaike weights and determine best model
            min_aic = min(Akaike_corr.values())
            deltas = {k: v - min_aic for k, v in Akaike_corr.items()}
            edeltas = {k: np.exp(-v/2) for k, v in deltas.items()}
            sum_edeltas = sum(edeltas.values())
            weights = {k: v/sum_edeltas for k, v in edeltas.items()}
            
            best_model = max(weights.items(), key=lambda x: x[1])[0]
            best_fits[sample_idx] = {'NLM': 0, 'PTM': 1, 'ETM': 2}[best_model]
            best_Ktrans[sample_idx] = Ktrans_fit[best_model]
            best_vp[sample_idx] = vp_fit[best_model]
            if best_model == 'ETM':
                best_ve[sample_idx] = ve_fit[best_model]
            else:
                best_ve[sample_idx] = np.nan  # NLM and PTM don't have ve
        
        # Save ground truth parameters
        np.savetxt(f'{base_path}ground_truth_T1_Repeat{repeat}.txt', T1_0_samples, fmt='%.6f')
        np.savetxt(f'{base_path}ground_truth_Fp_Repeat{repeat}.txt', F_p_samples, fmt='%.6f')
        np.savetxt(f'{base_path}ground_truth_PS_Repeat{repeat}.txt', PS_samples, fmt='%.6e')
        np.savetxt(f'{base_path}ground_truth_ve_Repeat{repeat}.txt', v_e_samples, fmt='%.6f')
        np.savetxt(f'{base_path}ground_truth_vp_Repeat{repeat}.txt', v_p_samples, fmt='%.6f')
        np.savetxt(f'{base_path}ground_truth_tau_a_Repeat{repeat}.txt', tau_a_samples, fmt='%.6f')
        
        # Save best fit results
        best_fits_showing = best_fits + 1  # Convert to 1-based indexing for saving
        np.savetxt(f'{base_path}BFModel_Repeat{repeat}.txt', best_fits_showing, fmt='%d')
        np.savetxt(f'{base_path}BFKtrans_Repeat{repeat}.txt', best_Ktrans, fmt='%.6e')
        np.savetxt(f'{base_path}BFvp_Repeat{repeat}.txt', best_vp, fmt='%.6f')
        np.savetxt(f'{base_path}BFve_Repeat{repeat}.txt', best_ve, fmt='%.6f')
        
        # Save all fitted parameters for each model
        np.savetxt(f'{base_path}NLMKtrans_Repeat{repeat}.txt', NLM_Ktrans, fmt='%.6e')
        np.savetxt(f'{base_path}PTMKtrans_Repeat{repeat}.txt', PTM_Ktrans, fmt='%.6e') 
        np.savetxt(f'{base_path}ETMKtrans_Repeat{repeat}.txt', ETM_Ktrans, fmt='%.6e')
        np.savetxt(f'{base_path}NLMvp_Repeat{repeat}.txt', NLM_vp, fmt='%.6f')
        np.savetxt(f'{base_path}PTMvp_Repeat{repeat}.txt', PTM_vp, fmt='%.6f')
        np.savetxt(f'{base_path}ETMvp_Repeat{repeat}.txt', ETM_vp, fmt='%.6f')
        np.savetxt(f'{base_path}ETMve_Repeat{repeat}.txt', ETM_ve, fmt='%.6f')
        
        print(f'Completed repeat {repeat + 1}/{n_repeats}')

if __name__ == '__main__':
    # Set random seed for reproducibility (optional)
    np.random.seed(42)
    main()