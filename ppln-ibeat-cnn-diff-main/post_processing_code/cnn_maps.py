import nibabel as nib
import numpy as np
from dipy.io import read_bvals_bvecs
from dipy.core.gradients import gradient_table
import dipy.reconst.dti as dti
from dipy.io.image import save_nifti
from scipy.interpolate import interp1d
from dipy.io.gradients import read_bvals_bvecs
import glob
import os
import zipfile
from scipy import ndimage

# 1. Define paths
# Make sure these point to your specific subject's bval/bvec
target_dir = r"C:\Users\eic20eh\Downloads\testing_outputs\test_results"
results_dir = r"C:\Users\eic20eh\Downloads\testing_outputs\test_results"
bval_path = r"C:\Users\eic20eh\Downloads\testing_outputs\test_results\iBE-2128_012_niis\files.bval"
bvec_path = r"C:\Users\eic20eh\Downloads\testing_outputs\test_results\iBE-2128_012_niis\files.bvec"

# The ID must match the prefix in your filename (e.g., iBE-2128-032_baseline)
subject_id = "iBE-2128_012_baseline" 

def extract_nested_zips(path):
    for root, dirs, files in os.walk(path):
        for file in files:
            if file.endswith(".zip"):
                zip_path = os.path.join(root, file)
                # Create a folder name based on the zip filename (without .zip)
                folder_name = file.replace(".zip", "")
                extract_path = os.path.join(root, folder_name)
                
                print(f"Extracting: {file} -> {folder_name}")
                
                with zipfile.ZipFile(zip_path, 'r') as zip_ref:
                    zip_ref.extractall(extract_path)
                
                # Optional: Delete the zip file after extraction to save space
                # os.remove(zip_path)

def reassemble_subject(subject_id, results_dir):
    print(f"Reassembling slices for {subject_id}...")
    
    # 1. Use ** to search through all sub-folders recursively
    # 2. Changed suffix to match your actual file: _test.nii.gz
    search_pattern = os.path.join(results_dir, "**", f"{subject_id}_slice*_warped*_test.nii.gz")
    
    # recursive=True is key here
    slice_files = sorted(glob.glob(search_pattern, recursive=True), 
                         key=lambda x: int(x.split('_slice')[-1].split('_')[0]))
    
    if not slice_files:
        print(f"No slices found for {subject_id} in {results_dir}")
        # Debugging tip: print what it WAS looking for
        print(f"Searched for: {search_pattern}") 
        return None

    # Load and stack along Z-axis (axis 2)
    volume_list = [nib.load(f).get_fdata() for f in slice_files]
    full_volume = np.stack(volume_list, axis=2)

    print(f"Reassembled Volume Shape: {full_volume.shape}") 

    ref_img = nib.load(slice_files[0])
    output_path = os.path.join(results_dir, f"{subject_id}_FULL_WARPED.nii.gz")
    
    # Save the reassembled 4D volume
    new_img = nib.Nifti1Image(full_volume, ref_img.affine, ref_img.header)
    nib.save(new_img, output_path)
    return output_path

# --- Execution ---

# 1. Reassemble Slices
full_volume_path = reassemble_subject(subject_id, results_dir)

if full_volume_path:
    # 2. Load the newly created full volume
    img = nib.load(full_volume_path)
    data = img.get_fdata()
    affine = img.affine


    # 1. Load original gradients (e.g., length 200)
    bvals, bvecs = read_bvals_bvecs(bval_path, bvec_path)
    original_length = len(bvals)
    target_length = 144

    # 2. Create the "Time" indices for interpolation
    original_indices = np.linspace(0, 1, original_length)
    target_indices = np.linspace(0, 1, target_length)

    # 3. Interpolate b-values
    # We use 'linear' to mimic MONAI's bilinear spatial stretching
    f_bvals = interp1d(original_indices, bvals, kind='linear')
    bvals_resized = f_bvals(target_indices)

    # 4. Interpolate b-vectors
    f_bvecs = interp1d(original_indices, bvecs, axis=0, kind='linear')
    bvecs_resized = f_bvecs(target_indices)

    # --- CRITICAL STEP: Re-normalization ---
    # Interpolating between two points on a sphere cuts through the sphere.
    # We must push the vectors back out to unit length (norm=1).
    norms = np.linalg.norm(bvecs_resized, axis=1, keepdims=True)
    # Avoid division by zero for b=0 frames
    norms[norms == 0] = 1 
    bvecs_resized = bvecs_resized / norms

    # 5. Create the Table
    gtab = gradient_table(bvals_resized, bvecs=bvecs_resized)

    # 5. Fit DTI Model
    print("Fitting DTI model...")
    tenmodel = dti.TensorModel(gtab)
    tenfit = tenmodel.fit(data)

    # 6. Save Maps
    print("Saving FA and MD maps...")
    save_nifti(os.path.join(results_dir, f"{subject_id}_FA.nii.gz"), tenfit.fa.astype(np.float32), affine)
    save_nifti(os.path.join(results_dir, f"{subject_id}_MD.nii.gz"), tenfit.md.astype(np.float32), affine)
    print("Done!")