import nibabel as nib
import numpy as np
import os
from scipy.ndimage import zoom

def calculate_correlation(moving_path, reference_path):
    """
    Calculates Pearson Correlation, ensuring the reference is downsampled 
    to match the moving image (registration output) resolution.
    """
    if not os.path.exists(moving_path) or not os.path.exists(reference_path):
        return None

    # Load data
    img_mov = nib.load(moving_path).get_fdata()
    img_ref = nib.load(reference_path).get_fdata()

    # 1. Collapse 4D to 3D (Middle volume)
    if img_mov.ndim == 4: img_mov = img_mov[..., img_mov.shape[-1] // 2]
    if img_ref.ndim == 4: img_ref = img_ref[..., img_ref.shape[-1] // 2]

    # 2. Resample REFERENCE to match MOVING resolution (e.g., 256 -> 128)
    # This prevents interpolation artifacts in your DL results.
    if img_ref.shape != img_mov.shape:
        print(f"Resampling reference {img_ref.shape} to match moving {img_mov.shape}...")
        factors = [t / s for t, s in zip(img_mov.shape, img_ref.shape)]
        img_ref = zoom(img_ref, factors, order=1)
        
        # Force exact match if rounding occurred
        if img_ref.shape != img_mov.shape:
            temp = np.zeros(img_mov.shape)
            slices = tuple(slice(0, min(img_ref.shape[i], img_mov.shape[i])) for i in range(len(img_mov.shape)))
            temp[slices] = img_ref[slices]
            img_ref = temp

    # 3. Flatten arrays to 1D vectors
    vec_mov = img_mov.flatten()
    vec_ref = img_ref.flatten()

    # 4. Calculate Pearson Correlation Coefficient
    correlation_matrix = np.corrcoef(vec_mov, vec_ref)
    return correlation_matrix[0, 1]

# --- Implementation ---
# Fixed the patient ID to 024 to match your logs
base_path = r"C:\Users\eic20eh\Downloads\testing_outputs\test_results\iBE-2128-028_niis"
fixed_ref = os.path.join(base_path, "dwi_registered.nii")
unreg     = os.path.join(base_path, "dwi_unregistered.nii")
dl_warped = os.path.join(base_path, "iBE-2128-028_baseline_FULL_WARPED.nii")

print("\n--- Calculating Correlation Coefficients (Downsampled to 128x128) ---")

# Compare DL Result (128x128) to Baseline (downsampled to 128x128)
dl_corr = calculate_correlation(dl_warped, unreg)
if dl_corr: print(f"DL_REG vs UN_REG: {dl_corr:.4f}")

# Compare Iterative Result (256x256) to Baseline (256x256)
# Note: Since these are both likely high-res, no resampling will occur here
op_corr = calculate_correlation(fixed_ref, unreg)
if op_corr: print(f"OP_REG vs UN_REG: {op_corr:.4f}")

# Direct comparison of the two registration methods (Iterative downsampled to 128)
comp_corr = calculate_correlation(dl_warped, fixed_ref)
if comp_corr: print(f"DL_REG vs OP_REG: {comp_corr:.4f}")