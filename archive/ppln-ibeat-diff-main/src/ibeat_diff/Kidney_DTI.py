"""
Kidney DTI-MRI case study processing pipeline.

What does it do?

This module processes a kidney DTI-MRI dataset and corrects motion artifacts using a group-wise registration approach.:

How to use?

1. Select python interpreter (ctr + Shift + P)
2. create a virtual environment using conda env create -f DTI_environment.yml    
3. Activate the environment using conda activate mc_env
4. Set the DATAPATH variable to point to the folder containing the DICOM data.
5. Run the script using python c:/Users/eic20eh/Downloads/ppln-ibeat-diff-main/src/ibeat_diff/Kidney_DTI.py. It will create a build/kidney_DTI folder in the current working directory.
6. All intermediate results and final outputs will be saved in the build/kidney_DTI folder

All steps will be executed sequentially when running the script
Steps 1-3 are automated, while steps 4-5 require user interaction. 
The script will pause and open a napari viewer for manual editing 
of segmentations and drawing of AIF. After editing, close the napari 
window to resume execution. All subsequent steps are automated again.
"""

import numpy as np
import napari
import matplotlib.pyplot as plt
import mdreg
import nibabel as nib
import os
from dipy.core.gradients import gradient_table
from dipy.reconst.dti import TensorModel
from dipy.reconst.dti import fractional_anisotropy
import subprocess
import argparse
import glob
import zipfile
import shutil

COREG_PATH = r"C:\Users\eic20eh\Downloads\ppln-ibeat-diff-main\src\ibeat_diff\Kidney_DTI.py"

def run_dcm2niix(dicom_dir, output_dir, dcm2niix_path):
    print(f"--- Starting DICOM conversion from: {dicom_dir} ---")
    os.makedirs(output_dir, exist_ok=True)
    # Added '-m', 'y' for merging Siemens/IMA stacks
    cmd = [dcm2niix_path, "-o", output_dir, "-f", "%f", "-z", "y", "-m", "y", dicom_dir]

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"ERROR: dcm2niix failed!\n{result.stderr}")
    else:
        print("DICOM conversion successful.")

def dwi_model(dwi_map, bvals_orig=None, bvecs_orig=None):

    if not hasattr(dwi_model, "has_run"):
        print(f"\n--- MODEL DEBUG INFO ---")
        print(f"dwi_map shape (from mdreg): {dwi_map.shape}")
        print(f"bvals_orig count: {len(bvals_orig)}")
        print(f"First 5 b-values: {bvals_orig[:5]}")
        print(f"bvecs_orig shape: {bvecs_orig.shape}")
        dwi_model.has_run = True

    if bvecs_orig.shape[0] == 3 and bvecs_orig.shape[1] != 3:
        bvecs_orig = bvecs_orig.T

    gtab = gradient_table(bvals_orig, bvecs=bvecs_orig)

    ten_model = TensorModel(gtab)
    ten_fit = ten_model.fit(dwi_map)

    return ten_fit

def mdreg_dwi_model(dwi_map, bvals_orig=None, bvecs_orig=None):

    gtab = gradient_table(bvals_orig, bvecs=bvecs_orig)
    ten_fit = dwi_model(dwi_map, bvals_orig, bvecs_orig)

    predicted = ten_fit.predict(gtab)
    
    if not np.all(np.isfinite(predicted)):
        print(f"!!! WARNING: Model predicted {np.isnan(predicted).sum()} NaNs !!!")

    # --- CRITICAL STABILITY FIX ---
    # Replace any NaNs or Infs with 0 or the original data
    predicted = np.nan_to_num(predicted, nan=0.0, posinf=0.0, neginf=0.0)
    # Clip values to ensure they stay within a realistic physical range
    predicted = np.clip(predicted, 0, np.max(dwi_map))
    
    params = ten_fit.model_params
    return predicted, params

def handle_zipped_dicoms(input_path, extract_to):
    # 1. If it's a zip, unzip it first
    if input_path.endswith('.zip'):
        unzip_dir = os.path.join(extract_to, "temp_extracted")
        if os.path.exists(unzip_dir):
            shutil.rmtree(unzip_dir)
        
        with zipfile.ZipFile(input_path, 'r') as zip_ref:
            zip_ref.extractall(unzip_dir)
        search_path = unzip_dir
    else:
        # It's already a folder, so we just search within it
        search_path = input_path

    for root, dirs, files in os.walk(search_path):
        if os.path.basename(root).upper() == "DICOM":
            # If 'DICOM' has a subfolder called 'files', go deeper
            if "files" in dirs:
                root = os.path.join(root, "files")
            return root

def save_nifti(data, path, affine):
    
    os.makedirs(os.path.dirname(path), exist_ok=True)
    data = np.flip(data, axis=1)
    img_to_save = nib.Nifti1Image(data.astype(np.float32), affine)
    nib.save(img_to_save, path)


def motion_correction(nii_obj, raw_data, coreg_path, bvals_orig, bvecs_orig):
    
    original_zooms = list(nii_obj.header.get_zooms())
    pixel_spacing = [original_zooms[1], original_zooms[0]]

    z=10
    raw_data = np.ascontiguousarray(raw_data[:, :, z, :])
    bvals_orig = bvals_orig[:]
    bvecs_orig = bvecs_orig[:, :]

    if bvecs_orig.shape[0] == 3:
        bvecs_orig = bvecs_orig.T

    print("Starting GROUP-WISE DWI motion correction")

    print(f"Input Shape: {raw_data.shape}") # Should be (172, 172, 145)

    coreg_affine, fit_affine, transfo_aff, pars_aff = mdreg.fit(
    raw_data,
    fit_image={
            'func': mdreg_dwi_model,
            'bvals_orig': bvals_orig,
            'bvecs_orig': bvecs_orig,
        },
        fit_coreg={
            'package': 'elastix',
            'method': 'affine',
            'spacing': pixel_spacing,
        },
        maxit=1,
        verbose=2,
        force_2d=True,
    )

    coreg, fit_final, transfo_bs, pars_bs = mdreg.fit(
    coreg_affine,
    fit_image={
            'func': mdreg_dwi_model,
            'bvals_orig': bvals_orig,
            'bvecs_orig': bvecs_orig,
        },
        fit_coreg={
            'package': 'elastix',
            'method': 'bspline',
            'spacing': pixel_spacing,
            'FinalGridSpacingInPhysicalUnits': 50,
        },
        maxit=1,
        verbose=2,
        force_2d=True,
    )


    print("Group-wise motion correction finished")
    
    # Save the coregistered 4D image as NIfTI
    nii_path = coreg_path.replace('.npy', '.nii.gz')
    nii_fit_path = coreg_path.replace('.npy', '_fit.nii.gz')
    save_nifti(coreg, nii_path, nii_obj.affine)
    np.save(coreg_path, coreg) # Saves DTI_coreg.npy
    np.save(coreg_path.replace('.npy', '_fit.npy'), fit_final)   

    return coreg, coreg_affine, fit_final, fit_affine, pars_bs, pars_aff


def view_unregistered_maps(bvals_orig, bvecs_orig, raw_data):

    if bvecs_orig.shape[0] == 3:
        bvecs_orig = bvecs_orig.T

    # Sanity check: print shapes
    # print(dwi.shape)
    # print(bvals.shape)
    # print(bvecs.shape)

    ten_fit = dwi_model(raw_data, bvals_orig, bvecs_orig)

    from dipy.reconst.dti import fractional_anisotropy
    FA = fractional_anisotropy(ten_fit.evals)
    FA = np.clip(FA, 0, 1)

    MD = ten_fit.md

    FA_napari = np.transpose(FA, (2, 1, 0))
    FA_napari = np.flip(FA_napari, axis=-2)
    MD_napari = np.transpose(MD, (2, 1, 0))
    MD_napari = np.flip(MD_napari, axis=-2)

    viewer = napari.Viewer()
    viewer.add_image(FA_napari, name="FA_unregistered", colormap="gray")
    viewer.add_image(MD_napari, name="MD_unregistered", colormap="gray", contrast_limits=[0, 0.002])

    viewer.layers[0].interactive = True
    viewer.layers[1].interactive = True

    napari.run()

def view_registered_maps(bvals_orig, bvecs_orig, coreg_data):

    if bvecs_orig.shape[0] == 3:
        bvecs_orig = bvecs_orig.T
    
    coreg_data = np.transpose(coreg_data, (2, 0, 1, 3))
    
    # Sanity check: print shapes
    # print(dwi.shape)
    # print(bvals.shape)
    # print(bvecs.shape)

    ten_fit = dwi_model(coreg_data, bvals_orig, bvecs_orig)

    from dipy.reconst.dti import fractional_anisotropy
    FA_coreg = np.clip(fractional_anisotropy(ten_fit.evals), 0, 1)

    MD_coreg = ten_fit.md

    FA_coreg_napari = np.transpose(FA_coreg, (1, 2, 0))
    FA_coreg_napari = np.rot90(FA_coreg_napari, k=1, axes=(-2, -1))
    MD_coreg_napari = np.transpose(MD_coreg, (1, 2, 0))
    MD_coreg_napari = np.rot90(MD_coreg_napari, k=1, axes=(-2, -1))

    viewer = napari.Viewer()
    viewer.add_image(FA_coreg_napari, name="FA_registered", colormap="gray")
    viewer.add_image(MD_coreg_napari, name="MD_registered", colormap="gray", contrast_limits=[0, 0.002])
    viewer.layers[0].interactive = True
    viewer.layers[1].interactive = True 

    napari.run()

def view_unregistered_Images(raw_data):

    viewer = napari.Viewer()
    raw_data = np.transpose(raw_data, (3, 2, 1, 0,)) # (Z, Volumes, Y, X)
    raw_data = np.flip(raw_data, axis=-2)
    viewer.add_image(
        raw_data,
        name="DTI_unregistered",
        axis_labels=("x", "y", "z", "volume"),
    )

    napari.run()

    # for i in range(len(bvals)):
    #     print(f"Volume {i}: b = {bvals[i]}, bvec = {bvecs[:, i]}")

    # print("DTI data shape:", data.shape)
    # print("Number of b-values:", len(bvals))

def view_registered_Images(coreg_data):

    viewer = napari.Viewer()
    coreg_data = np.transpose(coreg_data, (3, 0, 2, 1,))
    coreg_data = np.flip(coreg_data, axis=(-2))
    viewer.add_image(
        coreg_data,
        name="DTI_registered",
        axis_labels=("x", "y", "z", "volume"),
    )

    napari.run()

    # for i in range(len(bvals)):
    #     print(f"Volume {i}: b = {bvals[i]}, bvec = {bvecs[:, i]}")

    # print("DTI data shape:", data.shape)
    # print("Number of b-values:", len(bvals))

def view_model(fit_data, coreg_data):

    model_data = fit_data[:, :, :] 
    model_data = np.transpose(model_data, (2, 0, 1))
    model_data = np.rot90(model_data, k=1, axes=(-2,-1))

    coreg_data = np.transpose(coreg_data, (2, 0, 1))
    coreg_data = np.rot90(coreg_data, k=1, axes=(-2,-1))

    viewer = napari.Viewer()

    viewer.add_image(model_data)
    viewer.add_image(coreg_data)

    napari.run()

def map_download(bvals_orig, bvecs_orig, raw_data, coreg_data, fa_path, md_path, affine):

    tenfit_raw = dwi_model(raw_data, bvals_orig, bvecs_orig)
    tenfit_coreg = dwi_model(coreg_data, bvals_orig, bvecs_orig)

    FA_raw = fractional_anisotropy(tenfit_raw.evals)
    FA_coreg = fractional_anisotropy(tenfit_coreg.evals)
    FA_raw = np.clip(FA_raw, 0, 1)
    FA_coreg = np.clip(FA_coreg, 0, 1) 

    MD_raw = tenfit_raw.md
    MD_coreg = tenfit_coreg.md
    
    sub_dir = os.path.dirname(fa_path)

    save_nifti(FA_raw, os.path.join(sub_dir, 'FA_unregistered.nii.gz'), affine)
    save_nifti(MD_raw, os.path.join(sub_dir, 'MD_unregistered.nii.gz'), affine)
    save_nifti(FA_coreg, os.path.join(sub_dir, 'FA_registered.nii.gz'), affine)
    save_nifti(MD_coreg, os.path.join(sub_dir, 'MD_registered.nii.gz'), affine)

def map_comparison(bvals_orig, bvecs_orig, raw_data, coreg_data, fa_path, md_path):

    if bvecs_orig.shape[0] == 3:
        bvecs_orig = bvecs_orig.T

    tenfit_raw = dwi_model(raw_data, bvals_orig, bvecs_orig)
    tenfit_coreg = dwi_model(coreg_data, bvals_orig, bvecs_orig)

    FA_raw = fractional_anisotropy(tenfit_raw.evals)
    FA_coreg = fractional_anisotropy(tenfit_coreg.evals)
    FA_raw = np.clip(FA_raw, 0, 1)
    FA_coreg = np.clip(FA_coreg, 0, 1) 

    MD_raw = tenfit_raw.md
    MD_coreg = tenfit_coreg.md

    # --- ADAPTIVE DIMENSION CHECK ---
    def prepare_for_napari(data):
        # If 3D (X, Y, Z), transpose to (Z, Y, X)
        if data.ndim == 3:
            data = np.transpose(data, (2, 0, 1))
        # Rotate 90 degrees to get the kidneys upright (Radiological view)
        return np.rot90(data, k=1, axes=(-2, -1))

    FA_raw_napari = prepare_for_napari(FA_raw)
    MD_raw_napari = prepare_for_napari(MD_raw)
    FA_coreg_napari = prepare_for_napari(FA_coreg)
    MD_coreg_napari = prepare_for_napari(MD_coreg)

    viewer = napari.Viewer()
    viewer.add_image(FA_raw_napari, name="FA_unregistered", colormap="gray")
    viewer.add_image(FA_coreg_napari, name="FA_registered", colormap="gray")
    viewer.add_image(MD_raw_napari, name="MD_unregistered", colormap="gray",  contrast_limits=[0, 0.002])
    viewer.add_image(MD_coreg_napari, name="MD_registered", colormap="gray",  contrast_limits=[0, 0.002])
    viewer.dims.axis_labels = ("Slice", "Y", "X")

    viewer.layers[0].interactive = True
    viewer.layers[1].interactive = True  
    viewer.layers[2].interactive = True  
    viewer.layers[3].interactive = True  

    viewer.grid.enabled = True
    viewer.grid.shape = (2, 2)

    napari.run()

def image_download(raw_data, coreg_data, image_path, affine):
    sub_dir = os.path.dirname(image_path)
    
    # Save 4D volumes with 'dwi' prefix to distinguish from the FA maps
    save_nifti(coreg_data, os.path.join(sub_dir, 'dwi_registered.nii.gz'), affine)
    save_nifti(raw_data, os.path.join(sub_dir, 'dwi_unregistered.nii.gz'), affine)

def image_comparison(raw_data, coreg_data, image_path):

    viewer = napari.Viewer() #Order = (Slider 1, Slider 2, Y, X)

    if coreg_data.ndim == 4:
        # Standard 4D: Move Vol to front, then Z, Y, X
        axes = (3, 2, 0, 1)
    else:
        # Single slice 3D: Move Vol to front, then Y, X
        axes = (2, 0, 1)
        z=10
        raw_data = raw_data[:, :, z, :]
    
    raw_data = np.transpose(raw_data, axes) 
    raw_data = np.rot90(raw_data, k=1, axes=(-2, -1))
    coreg_data = np.transpose(coreg_data, axes)
    coreg_data = np.rot90(coreg_data, k=1, axes=(-2, -1))

    viewer.add_image(raw_data, name="DTI_unregistered",)
    viewer.add_image(coreg_data, name="DTI_registered",)
    # viewer.dims.axis_labels = ("b-value", "Slice", "Y", "X")

    viewer.grid.enabled = True
    viewer.grid.shape = (1, 2)
    # viewer.reset_view()

    napari.run()

def numerical_check(raw_data, coreg_data):
    
    coreg_data = np.transpose(coreg_data, (2, 1, 0, 3)) 
    diff = np.mean(np.abs(coreg_data - raw_data), axis=(0,1,2))

    plt.plot(diff)
    plt.xlabel("Volume")
    plt.ylabel("Mean |ΔI|")
    plt.title("Motion correction effect per volume")
    plt.show()

def view_all_stages(raw_data, model_fit, affine_data, final_data):
    """
    Opens a 4-pane dashboard in Napari showing the progression of the pipeline.
    """
    viewer = napari.Viewer(title="DTI Motion Correction Pipeline Stages")

    # Helper to format data for Napari (Z, Y, X, Vol)
    def prep(d):
        if d.ndim == 4:
            d = np.transpose(d, (3, 2, 0, 1))
        return np.rot90(d, k=1, axes=(-2, -1))

    viewer.add_image(prep(raw_data), name="1. Input Series (Raw)", colormap='gray')
    viewer.add_image(prep(model_fit), name="2. Model Prediction (Target)", colormap='viridis')
    viewer.add_image(prep(affine_data), name="3. Affine Registered", colormap='gray')
    viewer.add_image(prep(final_data), name="4. Fully Registered (B-Spline)", colormap='gray')

    viewer.grid.enabled = True
    viewer.grid.shape = (2, 2) # 2x2 grid for easy comparison
    napari.run()

if __name__=='__main__':

# 1. DEFINE YOUR SPECIFIC PATHS HERE

    # PATHS for RAW DATA
    RAW_DICOM = r"C:\Users\eic20eh\Downloads\data\iBE-2128-001_followup.zip"
   
    raw_data = os.path.basename(RAW_DICOM).replace('.zip', '')
    BUILD = os.path.join(r"C:\Users\eic20eh\Downloads\ppln-ibeat-diff-main\build", raw_data)

    # ARGUMENT SETUP  
    parser = argparse.ArgumentParser(description="Kidney DTI Pipeline")
    parser.add_argument('--dicom', type=str, default=RAW_DICOM)
    parser.add_argument('--build', type=str, default=BUILD)
    parser.add_argument('--dcm2niix', type=str, default="dcm2niix")
    args = parser.parse_args()

    # CONVERSION
    effective_dicom_path = handle_zipped_dicoms(args.dicom, args.build)
    run_dcm2niix(effective_dicom_path, BUILD, "dcm2niix")

    # PATHS TO INPUT FILES
    dti_path = glob.glob(os.path.join(args.build, "*.nii.gz"))[0]
    bval_path = glob.glob(os.path.join(args.build, "*.bval"))[0]
    bvec_path = glob.glob(os.path.join(args.build, "*.bvec"))[0]
    
    # VARIABLES
    niiobj = nib.load(dti_path)
    data = niiobj.get_fdata()
    bvals = np.loadtxt(bval_path)
    bvecs = np.loadtxt(bvec_path)

    # PATHS TO OUTPUT FILES
    resultspath = os.path.join(BUILD, "DTI_coreg.npy")
    fitpath = os.path.join(BUILD, "DTI_coreg_fit.npy")
    imagepath = os.path.join(BUILD, "DTI_coreg.nii.gz")
    fapath = os.path.join(BUILD, "FA_map.nii.gz")
    mdpath = os.path.join(BUILD, "MD_map.nii.gz")

    # MOTION CORRECTION
    # coreg, coreg_affine, model_fit, fit_affine, pars_bs, pars_aff = motion_correction(niiobj, data, resultspath, bvals, bvecs)
    # save_nifti(coreg, fapath, niiobj.affine)
    coreg = np.load(resultspath)
    fit = np.load(fitpath)

    # CHECKS
    # view_unregistered_maps(bvals, bvecs, data)
    # view_registered_maps(bvals, bvecs, coreg)
    # view_unregistered_Images(data)
    # view_registered_Images(coreg)
    # view_model(fit, coreg)
    # map_download(bvals, bvecs, data, coreg, fapath, mdpath, niiobj.affine)
    map_comparison(bvals, bvecs, data, coreg, fapath, mdpath)
    # image_download(data, coreg, imagepath, niiobj.affine)
    # image_comparison(data, coreg, imagepath)
    # numerical_check(data, coreg)
    # view_all_stages(data, model_fit, coreg_affine, coreg)

# if __name__ == '__main__':
#     parser = argparse.ArgumentParser()
#     # We now accept a single zip file path instead of an input directory
#     parser.add_argument('--input_zip', type=str, required=True)
#     parser.add_argument('--output_dir', type=str, required=True)
#     args = parser.parse_args()

#     # Define paths based on the specific zip file provided by Slurm
#     sub_path = args.input_zip
#     OUTPUT_ROOT = args.output_dir
    
#     # Get the subject ID from the filename (e.g., 'series_1401')
#     sub_id = os.path.basename(sub_path).replace('.zip', '')
#     sub_output_dir = os.path.join(OUTPUT_ROOT, sub_id)
    
#     try:
#         print(f"\n--- Starting Subject: {sub_id} ---")
        
#         # Create the folder for this specific subject
#         os.makedirs(sub_output_dir, exist_ok=True)

#         # Step A: DICOM Extraction
#         effective_dicom_path = handle_zipped_dicoms(sub_path, sub_output_dir)
#         run_dcm2niix(effective_dicom_path, sub_output_dir, "/mnt/parscratch/users/eic20eh/envs/mc_env/bin/dcm2niix")

#         # Step B: Load Data
#         dti_files = glob.glob(os.path.join(sub_output_dir, "*.nii.gz"))
#         bval_files = glob.glob(os.path.join(sub_output_dir, "*.bval"))
#         bvec_files = glob.glob(os.path.join(sub_output_dir, "*.bvec"))

#         if not dti_files or not bval_files or not bvec_files:
#              # Check if it was a 3D vs 4D issue
#              raise FileNotFoundError(f"Missing conversion files. Check if DICOMs are valid 4D DTI.")

#         dti_path = dti_files[0]
#         niiobj = nib.load(dti_path)
#         data = niiobj.get_fdata()

#         # DIMENSION CHECK: Catch the 3D vs 4D error before starting math
#         if data.ndim != 4:
#             raise ValueError(f"Expected 4D data, but got {data.ndim}D shape: {data.shape}")

#         bvals = np.loadtxt(bval_files[0])
#         bvecs = np.loadtxt(bvec_files[0])

#         # Step C: Motion Correction
#         resultspath = os.path.join(sub_output_dir, f"{sub_id}_coreg.npy")
#         coreg, pars = motion_correction(niiobj, data, resultspath, bvals, bvecs)

#         # Step D: Save Final Outputs
#         fapath = os.path.join(sub_output_dir, f"{sub_id}_FA.nii.gz")
#         mdpath = os.path.join(sub_output_dir, f"{sub_id}_MD.nii.gz")
#         map_download(bvals, bvecs, data, coreg, fapath, mdpath, niiobj.affine)
        
#         imgpath = os.path.join(sub_output_dir, f"{sub_id}_coreg.nii.gz")
#         image_download(data, coreg, imgpath, niiobj.affine)
        
#         print(f"SUCCESS: {sub_id} processed completely.")

#     except Exception as e:
#         print(f"ERROR processing {sub_id}: {e}")
#         # Clean up failed folders to keep /mnt/parscratch tidy
#         if os.path.exists(sub_output_dir):
#             shutil.rmtree(sub_output_dir)
        
#         # Log the error in a central file
#         with open(os.path.join(OUTPUT_ROOT, "error_log.txt"), "a") as f:
#             f.write(f"Subject {sub_id} failed: {str(e)}\n")

#         # FIX: The ignore_errors=True prevents the script from crashing during cleanup
#         if os.path.exists(sub_output_dir):
#             shutil.rmtree(sub_output_dir, ignore_errors=True)
            

#VIEW MAPS AND IMAGES
#.\.venv\Scripts\activate
#napari