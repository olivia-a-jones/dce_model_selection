"""
Model Select

Script to DCE-MRI model select based on the Akaike Information Criterion

Requires the model fit outputs from Madym, and a brain mask: patient_folder/structural/brain_mask_cut.nii

@author: Olivia A Jones
"""

import os
import numpy as np
from nibabel import load as load_nii, save as save_nii, Nifti1Image

def calculate_akaikes(id, model_name, k, n):
    """
    Calculate Akaike Information Criterion for Voxel-wise DCE-MRI.

    Args:
        id (str): Unique identifier for the patient.
        model_name (str): Name of the model, one of 'Intravascular', 'Patlak', or 'ETM'.
        k (int): Number of free parameters in the model.
        n (int): Number of dynamics (time points).

    Returns:
        float: Corrected Akaike Information Criterion.
    """
    rss_file = f'E:/DCE_Model_Selection/{id}/dce_fitting/{model_name}/residuals.nii'
    rss = np.double(load_nii(rss_file).get_fdata())

    akaike = 2 * (k + 1) + n * np.log(rss)
    akaike_corrected = akaike + ((2 * k) * (k + 1)) / (n - k - 1)
    return akaike_corrected

def make_best_fit_map(id):
    """
    Generate a map of the best fitting DCE-MRI model for each voxel.

    Args:
        id (str): Unique identifier for the patient.

    Returns:
        np.ndarray: Map of the best fitting model for each voxel.
        np.ndarray: Map of the Ktrans values for the best fitting model.
    """
    # Calculate Akaike Weights for each model
    n_timepoints = 160  # Replace with actual number of timepoints in your data
    akaikes_iv = calculate_akaikes(id, 'Intravascular', 1, n_timepoints)
    akaikes_ptm = calculate_akaikes(id, 'Patlak', 2, n_timepoints)
    akaikes_etm = calculate_akaikes(id, 'ETM', 3, n_timepoints)

    # Calculate delta AIC and weights
    delta_iv = akaikes_iv - min(akaikes_ptm, akaikes_etm)
    delta_ptm = akaikes_ptm - min(akaikes_iv, akaikes_etm)
    delta_etm = akaikes_etm - min(akaikes_iv, akaikes_ptm)

    e_delta_iv = np.exp(-delta_iv / 2)
    e_delta_ptm = np.exp(-delta_ptm / 2)
    e_delta_etm = np.exp(-delta_etm / 2)

    akaike_weight_iv = e_delta_iv / (e_delta_ptm + e_delta_etm)
    akaike_weight_ptm = e_delta_ptm / (e_delta_iv + e_delta_etm)
    akaike_weight_etm = e_delta_etm / (e_delta_iv + e_delta_ptm)

    # Save Akaike weight maps
    output_dir = f'E:/DCE_Model_Selection/{id}/dce_model_selection'
    os.makedirs(output_dir, exist_ok=True)
    
    save_nii(Nifti1Image(akaike_weight_iv, np.eye(4)), f'{output_dir}/Akaike_Weight_Intravascular.nii')
    save_nii(Nifti1Image(akaike_weight_ptm, np.eye(4)), f'{output_dir}/Akaike_Weight_Patlak.nii')
    save_nii(Nifti1Image(akaike_weight_etm, np.eye(4)), f'{output_dir}/Akaike_Weight_ETM.nii')

    # Load brain mask and error maps
    brain = np.double(load_nii(f'E:/DCE_Model_Selection/{id}/structural/brain_roi_cut.nii').get_fdata())
    error_map = np.zeros_like(brain)
    
    for model in ['Intravascular', 'Patlak', 'ETM']:
        error_data = load_nii(f'E:/DCE_Model_Selection/{id}/dce_fitting/{model}/error_tracker.nii').get_fdata()
        error_map[error_data > 0] = 1
    
    brain[error_map == 1] = np.nan
    save_nii(Nifti1Image(brain, np.eye(4)), f'{output_dir}/brain_mask.nii')

    # Load Ktrans and vp maps
    ktrans_iv = load_nii(f'E:/DCE_Model_Selection/{id}/dce_fitting/Intravascular/Ktrans.nii').get_fdata()
    ktrans_ptm = load_nii(f'E:/DCE_Model_Selection/{id}/dce_fitting/Patlak/Ktrans.nii').get_fdata()
    ktrans_etm = load_nii(f'E:/DCE_Model_Selection/{id}/dce_fitting/ETM/Ktrans.nii').get_fdata()
    vp_iv = load_nii(f'E:/DCE_Model_Selection/{id}/dce_fitting/Intravascular/v_p.nii').get_fdata()

    # Create best fitting model map and Ktrans map
    best_fits = np.full_like(brain, np.nan)
    ktrans_best_fit = np.full_like(brain, np.nan)

    # Determine best fitting model for each voxel
    for i in range(brain.shape[0]):
        for j in range(brain.shape[1]):
            for k in range(brain.shape[2]):
                if brain[i, j, k] == 1:
                    probs = [akaike_weight_iv[i, j, k], akaike_weight_ptm[i, j, k], akaike_weight_etm[i, j, k]]
                    best = np.max(probs)
                    if len(best)==1:
                        best_fit_ak = np.argmax(probs) + 1
                    else:
                        best_fit_ak = np.nan
                    
                    if np.isnan(best):
                        best_fits[i, j, k] = np.nan
                    elif vp_iv[i, j, k] < 0.0001:
                        best_fits[i, j, k] = 1  # Negligible vp = no model
                    else:
                        best_fits[i, j, k] = best_fit_ak
                
                # Assign appropriate Ktrans value based on best fitting model
                if best_fits[i, j, k] == 2:
                    ktrans_best_fit[i, j, k] = ktrans_iv[i, j, k]
                elif best_fits[i, j, k] == 3:
                    ktrans_best_fit[i, j, k] = ktrans_ptm[i, j, k]
                elif best_fits[i, j, k] == 4:
                    ktrans_best_fit[i, j, k] = ktrans_etm[i, j, k]
                elif best_fits[i, j, k] == 1:
                    ktrans_best_fit[i, j, k] = 0

    # Save best fit and Ktrans maps
    save_nii(Nifti1Image(best_fits, np.eye(4)), f'{output_dir}/best_fit.nii')
    save_nii(Nifti1Image(ktrans_best_fit, np.eye(4)), f'{output_dir}/Ktrans_best_fit.nii')
    
    return best_fits, ktrans_best_fit

def process_patient_data():
    """
    Process DCE-MRI data for all patient IDs and generate best fitting model maps.
    """
    with open('E:\DCE_Model_Selection\Participant_IDs.txt') as f:
        id_list = f.read().splitlines()

    for id in id_list:
        print(f"Processing patient {id}")
        make_best_fit_map(id)

if __name__ == "__main__":
    process_patient_data()