import os
import nibabel as nib
import napari
import numpy as np
from dipy.reconst.dti import TensorModel
from dipy.core.gradients import gradient_table
from dipy.io.gradients import read_bvals_bvecs
from scipy.interpolate import interp1d

base = r"C:\Users\eic20eh\Downloads\testing_outputs\test_results\iBE-2128_012_niis"
FA_name = "iBE-2128_012_baseline_FA.nii"
img_name = "iBE-2128_012_baseline_FULL_WARPED.nii"

def run_pipeline():
    # --- 1. DATA LOADING & TENSOR FITTING ---
    print("Loading raw DWI data and fitting tensors...")
    data_img = nib.load(os.path.join(base, img_name))
    data = data_img.get_fdata()
    
    # Load FA for masking
    fa_nii = nib.load(os.path.join(base, FA_name))
    fa_img = fa_nii.get_fdata()

    mask = (fa_img > 0.001)

    # Load bvals/bvecs
    bvals = os.path.join(base, "files.bval")
    bvecs = os.path.join(base, "files.bvec")
   
    # 1. Load original gradients (e.g., length 200)
    bvals, bvecs = read_bvals_bvecs(bvals, bvecs)
    original_length = len(bvals)
    target_length = 144

    # 2. Create the "Time" indices for interpolation
    original_indices = np.linspace(0, 1, original_length)
    target_indices = np.linspace(0, 1, target_length)

    # 3. Interpolate b-values
    # We use 'linear' to mimic MONAI's bilinear spatial stretching
    f_bvals = interp1d(original_indices, bvals, kind='linear')
    bvals = f_bvals(target_indices)

    # 4. Interpolate b-vectors
    f_bvecs = interp1d(original_indices, bvecs, axis=0, kind='linear')
    bvecs_resized = f_bvecs(target_indices)

    # --- CRITICAL STEP: Re-normalization ---
    # Interpolating between two points on a sphere cuts through the sphere.
    # We must push the vectors back out to unit length (norm=1).
    norms = np.linalg.norm(bvecs_resized, axis=1, keepdims=True)
    # Avoid division by zero for b=0 frames
    norms[norms == 0] = 1 
    bvecs = bvecs_resized / norms

    gtab = gradient_table(bvals, bvecs=bvecs)

    # Fit the Tensor and extract V1
    tenmodel = TensorModel(gtab)
    tenfit = tenmodel.fit(data, mask)

    # 1. Extract V1 (direction) and L1 (magnitude/size)
    v1_data = tenfit.evecs[:, :, :, :, 0]  # (X, Y, Z, 3)
    l1_data = tenfit.evals[:, :, :, 0]      # (X, Y, Z) - The principal eigenvalue
    
    # 2. Define your filters
    # You can combine FA and L1 thresholds for a very clean mask
    l1_threshold = 0.0005 # Adjust based on your data's noise floor
    fa_threshold = 0.0015
    
    vis_mask = (fa_img > fa_threshold) & (l1_data > l1_threshold)
    
    # 3. Prepare coordinates and vectors
    coords = np.argwhere(vis_mask)
    
    # Extract unit vectors for the masked area
    unit_vectors = v1_data[vis_mask]
    
    # Scale them by L1 (and an optional visualization constant)
    scale_factor = 1000 # Since L1 is usually very small
    magnitudes = l1_data[vis_mask] * scale_factor
    
    min_length = 1.5  # Adjust this value to remove "short" noise vectors
    length_mask = magnitudes > min_length
    
    # Apply the length mask to everything
    filtered_coords = coords[length_mask]
    filtered_magnitudes = magnitudes[length_mask]
    filtered_unit_vectors = unit_vectors[length_mask]

    # Apply scaling: (N, 3) * (N, 1)
    scaled_vectors = filtered_unit_vectors * filtered_magnitudes[:, np.newaxis]

    # Format for napari: (N, 2, 3)
    vec_display = np.zeros((len(filtered_coords), 2, 3))
    vec_display[:, 0, :] = filtered_coords   # Start point
    vec_display[:, 1, :] = scaled_vectors  # Direction

    # --- 3. VISUALIZATION ---
    print("Launching napari...")
    viewer = napari.Viewer()
    viewer.add_image(fa_img, name='FA Map', colormap='gray')
    
    rgb_colors = np.abs(scaled_vectors)

    # Add vectors layer
    viewer.add_vectors(
        vec_display, 
        name='V1 Glyphs', 
        edge_width=2.0, 
        length=1.0, 
        edge_color=rgb_colors # Or use RGB based on orientation
    )

    napari.run()

if __name__ == "__main__":
    run_pipeline()