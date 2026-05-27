# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 15:30:59 2025

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
            'init_params': np.ones((n_samples, 3)) * [0, 3.000000e-02, 0],
            'bounds': (
                np.ones((n_samples, 3)) * [0, 0, -20],
                np.ones((n_samples, 3)) * [1, 1, 20]
            )
        },
        'PTM': {
            'model': 'PATLAK',
            'fixed_params_uncut': None,
            'fixed_values_uncut': None,
            'fixed_params_cut': [3],
            'fixed_values_cut': [0],
            'init_params': np.ones((n_samples, 3)) * [1.000000e-03, 3.000000e-02, 0],
            'bounds': (
                np.ones((n_samples, 3)) * [-1, 0, -20],
                np.ones((n_samples, 3)) * [1, 1, 20]
            )
        },
        'ETM': {
            'model': 'ETM',
            'fixed_params_uncut': None,
            'fixed_values_uncut': None,
            'fixed_params_cut': [4],
            'fixed_values_cut': [0],
            'init_params': np.ones((n_samples, 4)) * [1.000000e-03, 3.000000e-01, 2.000000e-02, 0],
            'bounds': (
                np.ones((n_samples, 4)) * [-1, 0, 0, -20],
                np.ones((n_samples, 4)) * [1, 1, 1, 20]
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
    PS_samples = 10
    PS = np.concatenate((np.linspace(1e-4, 1e-3, PS_samples-1, endpoint=False), 
                        np.linspace(1e-3, 1e-2, PS_samples)))
    n_samples = len(PS)
    # # GM
    # Fixed parameters
    tissue = 'gm'
    T1_assumption = 1500
    F_p = np.ones(n_samples) * 0.8
    v_e = np.ones(n_samples) * 0.20
    v_p = np.ones(n_samples) * 0.04
    tau_a = np.ones(n_samples) * 0.00

    # WM
    # Fixed parameters
    # tissue = 'wm'
    # T1_assumption = 1000
    # F_p = np.ones(n_samples) * 0.4
    # v_e = np.ones(n_samples) * 0.20
    # v_p = np.ones(n_samples) * 0.02
    # tau_a = np.ones(n_samples) * 0.00
    
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
    
    for repeat in range(0,10):
        # Generate concentration time series
        Ct_input = np.zeros((n_samples, n_times))
        init_params = np.column_stack((F_p, PS, v_e, v_p, tau_a))
        
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
        S_t = tissue_concentration.concentration_to_signal(
            C_t=C_t,
            T1_0=T1_0,
            M0=S_t0,
            FA=FA,
            TR=TR,
            relax_coeff=r1_const,
            use_M0_ratio=8
        )
        
        # Generate noisy signals
        noise_list = np.linspace(0, 6, n_samples)
        S_tn = np.array([add_rician_noise(S_t, sig) for sig in noise_list])
        
        # Process each noise level in parallel
        fit_results_allsigma = []
        
        for sigma in range(n_samples):
            # Prepare arguments for parallel processing
            args_list = [(S_tn[sigma], model_name, t, T1_0, FA, r1_const, TR, n_samples, n_times) 
                        for model_name in model_names]
            
            # Process models in parallel
            with concurrent.futures.ProcessPoolExecutor() as executor:
                results = list(executor.map(process_model_fit, args_list))
            
            # Store results in dictionary format
            fit_results = {model_name: result for model_name, result in zip(model_names, results)}
            fit_results_allsigma.append(fit_results)
        
        # Calculate metrics and save results
        best_fits = np.zeros((n_samples, n_samples))
        best_Ktrans = np.zeros((n_samples, n_samples))
        PTM_Ktrans = np.zeros((n_samples, n_samples))
        ETM_Ktrans = np.zeros((n_samples, n_samples))
        
        for noises in range(n_samples):
            for params in range(n_samples):
                Akaike_corr = {}
                Ktrans_fit = {}
                
                for model_name in model_names:
                    results = fit_results_allsigma[noises][model_name]
                    RSS = results[1][params]
                    model_params = results[0][params]
                    N = n_times - 15
                    Ktrans_fit[model_name] = model_params[0]
                    
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
                best_fits[noises, params] = {'NLM': 0, 'PTM': 1, 'ETM': 2}[best_model]
                best_Ktrans[noises, params] = Ktrans_fit[best_model]
                PTM_Ktrans[noises, params] = Ktrans_fit['PTM']
                ETM_Ktrans[noises, params] = Ktrans_fit['ETM']
        
        # Save results
        base_path = f'/Users/user/Documents/Python_scripts/MODELLING_SIMULATIONS/Paper_Draft_2/opt_varSim_cut40_{tissue}/'
        if os.path.isdir(base_path)==False:
            os.mkdir(base_path)
        ground_truth_PS = np.tile(PS, (n_samples, 1))
        noise_grid = np.tile(noise_list, (n_samples, 1)).T
        best_fits_showing = best_fits + 1
        
        np.savetxt(f'{base_path}ground_truth_PS.txt', ground_truth_PS, fmt='%f')
        np.savetxt(f'{base_path}applied_noise.txt', noise_grid, fmt='%f')
        np.savetxt(f'{base_path}BFModel_Repeat{repeat}.txt', best_fits_showing, fmt='%f')
        np.savetxt(f'{base_path}BFKtrans_Repeat{repeat}.txt', best_Ktrans, fmt='%f')
        np.savetxt(f'{base_path}PTMKtrans_Repeat{repeat}.txt', PTM_Ktrans, fmt='%f')
        np.savetxt(f'{base_path}ETMKtrans_Repeat{repeat}.txt', ETM_Ktrans, fmt='%f')
        
        print(f'Processed repeat {repeat}...')
if __name__ == '__main__':
    main()