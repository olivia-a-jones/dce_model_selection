% Fit models to DCE-MRI data.
data_dir = 'E:/DCE_Model_Selection/';
ids_file = fullfile(data_dir,'Patient_IDs.txt');
hct_file = fullfile(data_dir,'Hematocrit.txt');

id_list = strsplit(ids_file);
hct_list = strsplit(hct_file);
for person=1:length(id_list)
    id = char(id_list(person))
    hct = double(hct_list(person));
    participant_dir = fullfile(data_dir,id);

    % Intravascular
    run_madym_DCE(...
        'cmd_exe', [local_madym_root 'madym_DCE'],...
        'model', 'PATLAK',... Model to fit, specified by its name in CAPITALS, see notes for options
        'output_dir', fullfile(participant_dir,'dce_fitting','Intravascular'), ...Output path, will use temp dir if empty;
        'no_optimise', 0, ...Flag to switch off optimising, will just fit initial parameters values for model
        'dynamic_basename', fullfile(dce_input_base, 'dynamic_series', 'rdyn_'), ...Template name for dynamic sequences eg. dynamic/dyn_
        'sequence_format', '%01u',...Format for converting dynamic series index to string, eg %01u
        'sequence_start', 2,...Start index of dynamic series
        'sequence_step', 1,...Step size between indexes in dynamic series
        'n_dyns', 159,...Number of dynamic sequence maps to load. If <=0, loads all maps in dynamic dir matching -dyn pattern
        'input_Ct', 0,...Flag specifying input dynamic sequence are concentration (not signal) maps
        'output_Ct_sig', 1,...Flag requesting concentration (derived from signal) are saved to output
        'output_Ct_mod', 1,...Flag requesting modelled concentration maps are saved to output
        'T1_name', fullfile(participant_dir, 'T1_map', 'VFA', 'T1.nii'),...Path to T1 map
        'M0_name', fullfile(participant_dir, 'T1_map', 'VFA', 'M0.nii'),...Path to M0 map
        'r1_const', 3.4,...Relaxivity constant of concentration in tissue (in ms)
        'M0_ratio', 0,...Flag to use ratio method to scale signal instead of supplying M0
        'dose', 0.1,...Concentration dose (mmole/kg)
        'injection_image', 8,...Injection image
        'hct', hct,...Haematocrit correction
        'aif_map', fullfile(participant_dir, 'dynamic_series', 'lrdyn_160.nii'),...Path to mask from which AIF will be computed on the fly
        'init_params', [1.000000e-03,2.000000e-02,0],...Initial values for model parameters to be optimised, either as single vector, or 2D array NSamples x N_params
        'fixed_params', [1],...Parameters fixed to their initial values (ie not optimised)
        'fixed_values', [0],..._values for fixed parameters (overrides default initial parameter values)
        'upper_bounds', [1,1,20],...Upper bounds for each parameter during optimisation
        'lower_bounds', [-1.000000e-03,0,-20],...Lower bounds for each parameter during optimisation
        'overwrite', true,...Set overwrite existing analysis in output dir
        'max_iter', 1000,... Maximum number of iterations in model fit
        'opt_type', 'BLEIC',... Type of optimisation to run
        'dummy_run', false);...Don't run any thing, just print the cmd we'll run to inspect

    % Patlak
    run_madym_DCE(...
        'cmd_exe', [local_madym_root 'madym_DCE'],...
        'model', 'PATLAK',... Model to fit, specified by its name in CAPITALS, see notes for options
        'output_dir', fullfile(participant_dir,'dce_fitting','Patlak'), ...Output path, will use temp dir if empty;
        'no_optimise', 0, ...Flag to switch off optimising, will just fit initial parameters values for model
        'dynamic_basename', fullfile(dce_input_base, 'dynamic_series', 'rdyn_'), ...Template name for dynamic sequences eg. dynamic/dyn_
        'sequence_format', '%01u',...Format for converting dynamic series index to string, eg %01u
        'sequence_start', 2,...Start index of dynamic series
        'sequence_step', 1,...Step size between indexes in dynamic series
        'n_dyns', 159,...Number of dynamic sequence maps to load. If <=0, loads all maps in dynamic dir matching -dyn pattern
        'input_Ct', 0,...Flag specifying input dynamic sequence are concentration (not signal) maps
        'output_Ct_sig', 1,...Flag requesting concentration (derived from signal) are saved to output
        'output_Ct_mod', 1,...Flag requesting modelled concentration maps are saved to output
        'T1_name', fullfile(participant_dir, 'T1_map', 'VFA', 'T1.nii'),...Path to T1 map
        'M0_name', fullfile(participant_dir, 'T1_map', 'VFA', 'M0.nii'),...Path to M0 map
        'r1_const', 3.4,...Relaxivity constant of concentration in tissue (in ms)
        'M0_ratio', 0,...Flag to use ratio method to scale signal instead of supplying M0
        'dose', 0.1,...Concentration dose (mmole/kg)
        'injection_image', 8,...Injection image
        'hct', hct,...Haematocrit correction
        'aif_map', fullfile(participant_dir, 'dynamic_series', 'lrdyn_160.nii'),...Path to mask from which AIF will be computed on the fly
        'init_params', [1.000000e-03,2.000000e-02,0],...Initial values for model parameters to be optimised, either as single vector, or 2D array NSamples x N_params
        'fixed_params', [],...Parameters fixed to their initial values (ie not optimised)
        'fixed_values', [],..._values for fixed parameters (overrides default initial parameter values)
        'upper_bounds', [1,1,20],...Upper bounds for each parameter during optimisation
        'lower_bounds', [-1.000000e-03,0,-20],...Lower bounds for each parameter during optimisation
        'overwrite', true,...Set overwrite existing analysis in output dir
        'max_iter', 1000,... Maximum number of iterations in model fit
        'opt_type', 'BLEIC',... Type of optimisation to run
        'dummy_run', false);...Don't run any thing, just print the cmd we'll run to inspect

    % Extended Tofts
    run_madym_DCE(...
        'cmd_exe', [local_madym_root 'madym_DCE'],...
        'model', 'PATLAK',... Model to fit, specified by its name in CAPITALS, see notes for options
        'output_dir', fullfile(participant_dir,'dce_fitting','Patlak'), ...Output path, will use temp dir if empty;
        'no_optimise', 0, ...Flag to switch off optimising, will just fit initial parameters values for model
        'dynamic_basename', fullfile(dce_input_base, 'dynamic_series', 'rdyn_'), ...Template name for dynamic sequences eg. dynamic/dyn_
        'sequence_format', '%01u',...Format for converting dynamic series index to string, eg %01u
        'sequence_start', 2,...Start index of dynamic series
        'sequence_step', 1,...Step size between indexes in dynamic series
        'n_dyns', 159,...Number of dynamic sequence maps to load. If <=0, loads all maps in dynamic dir matching -dyn pattern
        'input_Ct', 0,...Flag specifying input dynamic sequence are concentration (not signal) maps
        'output_Ct_sig', 1,...Flag requesting concentration (derived from signal) are saved to output
        'output_Ct_mod', 1,...Flag requesting modelled concentration maps are saved to output
        'T1_name', fullfile(participant_dir, 'T1_map', 'VFA', 'T1.nii'),...Path to T1 map
        'M0_name', fullfile(participant_dir, 'T1_map', 'VFA', 'M0.nii'),...Path to M0 map
        'r1_const', 3.4,...Relaxivity constant of concentration in tissue (in ms)
        'M0_ratio', 0,...Flag to use ratio method to scale signal instead of supplying M0
        'dose', 0.1,...Concentration dose (mmole/kg)
        'injection_image', 8,...Injection image
        'hct', hct,...Haematocrit correction
        'aif_map', fullfile(participant_dir, 'dynamic_series', 'lrdyn_160.nii'),...Path to mask from which AIF will be computed on the fly
        'init_params', [1.000000e-03,2.000000e-02,1.000000e-01,0],...Initial values for model parameters to be optimised, either as single vector, or 2D array NSamples x N_params
        'fixed_params', [],...Parameters fixed to their initial values (ie not optimised)
        'fixed_values', [],..._values for fixed parameters (overrides default initial parameter values)
        'upper_bounds', [1,1,1,20],...Upper bounds for each parameter during optimisation
        'lower_bounds', [-1.000000e-03,0,0,-20],...Lower bounds for each parameter during optimisation
        'overwrite', true,...Set overwrite existing analysis in output dir
        'max_iter', 1000,... Maximum number of iterations in model fit
        'opt_type', 'BLEIC',... Type of optimisation to run
        'dummy_run', false);...Don't run any thing, just print the cmd we'll run to inspect
end