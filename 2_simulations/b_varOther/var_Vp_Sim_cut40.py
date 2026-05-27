#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 17 16:35:18 2025

@author: p15094oj
"""

import numpy as np
from QbiMadym import madym_DCE_lite
from QbiPy.dce_models import tissue_concentration
from QbiMadym.madym_DWI_lite import add_rician_noise
import time
import concurrent.futures
import os

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
    # Initialize parameters
    n_samples = 5
    PS = np.array((1e-4, 5e-4, 1e-3, 5e-3, 1e-2))  # vary PS from 0.0001 to 0.01
    v_p_list = np.array((0.01, 0.02, 0.04, 0.06, 0.1))  # vary v_p from 1% (low vascularity) to 10% (very high vascularity)
    
    ### GM
    ## Fixed parameters
    # tissue = 'gm'
    # T1_assumption = 1500
    # F_p = 0.4  # Fixed Fp
    # v_e = np.ones(n_samples) * 0.20
    # tau_a = np.ones(n_samples) * 0.00

    ### WM
    ## Fixed parameters
    tissue = 'wm'
    T1_assumption = 1000
    F_p = 0.4  # Fixed Fp
    v_e = np.ones(n_samples) * 0.20
    tau_a = np.ones(n_samples) * 0.00
    
    # Time points
    n_times = 160
    temp_res_mins = 7.64/60
    t = np.arange(n_times) * temp_res_mins
    
    # Signal parameters
    T1_0 = np.ones(n_samples) * T1_assumption
    r1_const = 3.4
    FA = 10
    TR = 2.4
    S_t0 = np.ones(n_samples) * 100
    
    model_names = ['NLM', 'PTM', 'ETM']
    
    for repeat in range(0, 100000):
        # Generate concentration time series for entire grid
        init_params_all = []
        for v_p in v_p_list:
            F_p_array = np.ones(n_samples) * F_p  # Fixed Fp
            v_p_array = np.ones(n_samples) * v_p  # Varying vp
            init_params_row = np.column_stack((F_p_array, PS, v_e, v_p_array, tau_a))
            init_params_all.append(init_params_row)
        
        init_params = np.vstack(init_params_all)  # Stack all rows
        n_total = len(init_params)  # Total number of combinations (25)
        
        # Generate concentration time series for all combinations at once
        Ct_input = np.zeros((n_total, n_times))
        T1_0_all = np.ones(n_total) * T1_assumption
        
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
        
        # Convert to signal
        S_t0_all = np.ones(n_total) * 100
        S_t = tissue_concentration.concentration_to_signal(
            C_t=C_t,
            T1_0=T1_0_all,
            M0=S_t0_all,
            FA=FA,
            TR=TR,
            relax_coeff=r1_const,
            use_M0_ratio=8
        )
        
        # Generate noisy signals - different approach: apply same noise to all
        S_tn = add_rician_noise(S_t, 3.0)  # Fixed noise level
        
        # Process the entire grid in parallel - one call per model
        fit_results_allsigma = []
        
        # Prepare arguments for parallel processing - each model gets the entire grid
        args_list = [(S_tn, model_name, t, T1_0_all, FA, r1_const, TR, n_total, n_times) 
                    for model_name in model_names]
        
        # Process models in parallel - each processes the entire 25-combination grid
        with concurrent.futures.ProcessPoolExecutor() as executor:
            results = list(executor.map(process_model_fit, args_list))
        
        # Store results in dictionary format
        fit_results = {model_name: result for model_name, result in zip(model_names, results)}
        
        # Calculate metrics and save results
        best_fits = np.zeros((n_samples, n_samples))
        best_Ktrans = np.zeros((n_samples, n_samples))
        best_vp = np.zeros((n_samples, n_samples))
        best_ve = np.zeros((n_samples, n_samples))
        PTM_Ktrans = np.zeros((n_samples, n_samples))
        ETM_Ktrans = np.zeros((n_samples, n_samples))
        NLM_vp = np.zeros((n_samples, n_samples))
        PTM_vp = np.zeros((n_samples, n_samples))
        ETM_vp = np.zeros((n_samples, n_samples))
        ETM_ve = np.zeros((n_samples, n_samples))
        
        for vp_idx in range(n_samples):
            for ps_idx in range(n_samples):
                params_idx = vp_idx * n_samples + ps_idx  # Linear index into results
                
                Akaike_corr = {}
                Ktrans_fit = {}
                vp_fit = {}
                ve_fit = {}
                
                for model_name in model_names:
                    results = fit_results[model_name]
                    RSS = results[1][params_idx]
                    model_params = results[0][params_idx]
                    N = n_times - 15
                    Ktrans_fit[model_name] = model_params[0]
                    
                    # Extract vp and ve parameters based on model
                    if model_name == 'PTM':
                        vp_fit[model_name] = model_params[1]  # PTM fits vp (Patlak: Ktrans, vp)
                    elif model_name == 'ETM':
                        ve_fit[model_name] = model_params[1]
                        vp_fit[model_name] = model_params[2]
                    else:  # NLM
                        vp_fit[model_name] = model_params[1]  # NLM fits vp (Patlak: Ktrans, vp)
                    
                    K = {'NLM': 1, 'PTM': 2, 'ETM': 3}[model_name]
                    Akaike = (2 * (K + 1)) + (N * np.log(RSS))
                    Akaike_corr[model_name] = Akaike + ((2 * K) * (K + 1)) / (N - K - 1)
                
                # Calculate Akaike weights
                deltas = {
                    'NLM': Akaike_corr['NLM'] - min(Akaike_corr['PTM'], Akaike_corr['ETM']),
                    'PTM': Akaike_corr['PTM'] - min(Akaike_corr['NLM'], Akaike_corr['ETM']),
                    'ETM': Akaike_corr['ETM'] - min(Akaike_corr['NLM'], Akaike_corr['PTM'])
                }
                
                edeltas = {k: np.exp(-v/2) for k, v in deltas.items()}
                weights = {
                    'NLM': edeltas['NLM'] / (edeltas['PTM'] + edeltas['ETM']),
                    'PTM': edeltas['PTM'] / (edeltas['NLM'] + edeltas['ETM']),
                    'ETM': edeltas['ETM'] / (edeltas['NLM'] + edeltas['PTM'])
                }
                
                best_model = max(weights.items(), key=lambda x: x[1])[0]
                best_fits[vp_idx, ps_idx] = {'NLM': 0, 'PTM': 1, 'ETM': 2}[best_model]
                best_Ktrans[vp_idx, ps_idx] = Ktrans_fit[best_model]
                best_vp[vp_idx, ps_idx] = vp_fit[best_model]
                if best_model == 'ETM':
                    best_ve[vp_idx, ps_idx] = ve_fit[best_model]
                else:
                    best_ve[vp_idx, ps_idx] = np.nan  # NLM and PTM don't have ve
                PTM_Ktrans[vp_idx, ps_idx] = Ktrans_fit['PTM']
                ETM_Ktrans[vp_idx, ps_idx] = Ktrans_fit['ETM']
                NLM_vp[vp_idx, ps_idx] = vp_fit['NLM']
                PTM_vp[vp_idx, ps_idx] = vp_fit['PTM']
                ETM_vp[vp_idx, ps_idx] = vp_fit['ETM']
                ETM_ve[vp_idx, ps_idx] = ve_fit['ETM']
        
        # Save results
        base_path = f'/Users/user/Documents/Python_scripts/MODELLING_SIMULATIONS/Paper_Draft_2/varVpcut40_{tissue}/'
        if os.path.isdir(base_path) == False:
            os.mkdir(base_path)
            
        # Create input parameter grids for saving
        ground_truth_PS = np.tile(PS, (n_samples, 1))
        ground_truth_vp = np.tile(v_p_list, (n_samples, 1)).T
        best_fits_showing = best_fits + 1
        
        # Save all results
        np.savetxt(f'{base_path}ground_truth_PS.txt', ground_truth_PS, fmt='%f')
        np.savetxt(f'{base_path}ground_truth_vp.txt', ground_truth_vp, fmt='%f')
        np.savetxt(f'{base_path}BFModel_Repeat{repeat}.txt', best_fits_showing, fmt='%f')
        np.savetxt(f'{base_path}BFKtrans_Repeat{repeat}.txt', best_Ktrans, fmt='%f')
        np.savetxt(f'{base_path}BFvp_Repeat{repeat}.txt', best_vp, fmt='%f')
        np.savetxt(f'{base_path}BFve_Repeat{repeat}.txt', best_ve, fmt='%f')
        np.savetxt(f'{base_path}PTMKtrans_Repeat{repeat}.txt', PTM_Ktrans, fmt='%f')
        np.savetxt(f'{base_path}ETMKtrans_Repeat{repeat}.txt', ETM_Ktrans, fmt='%f')
        np.savetxt(f'{base_path}NLMvp_Repeat{repeat}.txt', NLM_vp, fmt='%f')
        np.savetxt(f'{base_path}PTMvp_Repeat{repeat}.txt', PTM_vp, fmt='%f')
        np.savetxt(f'{base_path}ETMvp_Repeat{repeat}.txt', ETM_vp, fmt='%f')
        np.savetxt(f'{base_path}ETMve_Repeat{repeat}.txt', ETM_ve, fmt='%f')
        
        print(f'Processed repeat {repeat}...')

if __name__ == '__main__':
    main()