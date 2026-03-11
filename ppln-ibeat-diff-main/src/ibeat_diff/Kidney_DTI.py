"""
Kidney DTI-MRI case study processing pipeline.

What does it do?

This module processes a kidney DTI-MRI dataset and corrects motion artifacts using a group-wise registration approach.:

How to use?

1. Select python interpreter (ctr + Shift + P)
2. create a virtual environment using conda env create -f DTI_environment.yml    
3. Activate the environment using conda activate kidney_dti_env
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

COREG_PATH = r"C:\Users\eic20eh\Downloads\case-studies-main\code\src\Kidney_DTI.py"

def dwi_model(dwi_map, bvals_orig=None, bvecs_orig=None):

    gtab = gradient_table(bvals_orig, bvecs=bvecs_orig)

    ten_model = TensorModel(gtab)
    ten_fit = ten_model.fit(dwi_map)
    # ten_fit = ten_fit.predict(gtab)

    return ten_fit

def mdreg_dwi_model(dwi_map, bvals_orig=None, bvecs_orig=None):

    gtab = gradient_table(bvals_orig, bvecs=bvecs_orig)

    ten_fit = dwi_model(dwi_map, bvals_orig, bvecs_orig)

    predicted = ten_fit.predict(gtab)
    params = ten_fit.model_params
    return predicted, params

def save_nifti(data, path):
    
    os.makedirs(os.path.dirname(path), exist_ok=True)

    data = np.flip(data, axis=1)

    affine = np.diag([2.0, 2.0, 4.0, 1.0])

    img_to_save = nib.Nifti1Image(data.astype(np.float32), affine)
    
    # 2. Save the object, not the raw array
    nib.save(img_to_save, path)

def run_dcm2niix(dicom_dir, output_dir, dcm2niix_path):
        
    print(f"--- Starting DICOM conversion from: {dicom_dir} ---")
    
    cmd = [dcm2niix_path, "-o", output_dir, "-f", "dti", "-z", "y", dicom_dir]
    
    result = subprocess.run(cmd, capture_output=True, text=True, shell=True)


def motion_correction(raw_data, coreg_path, bvals_orig, bvecs_orig):
    
    assert raw_data.ndim == 4, "Expected 4D DWI data (x, y, z, volumes)"

    #one slice
    z=10
    raw_data = raw_data[:, :, z, :]
    bvals_orig = bvals_orig[:]
    bvecs_orig = bvecs_orig[:, :] 

    print("Starting GROUP-WISE DWI motion correction")

    coreg, fit, transfo, pars = mdreg.fit(
    raw_data,
    fit_image={
            'func': mdreg_dwi_model,
            'bvals_orig': bvals_orig,
            'bvecs_orig': bvecs_orig,
        },
        fit_coreg={
            'package': 'elastix', #defult is b-spline
            'spacing': [2.3255813121796, 2.3255813121796],
            'FinalGridSpacingInPhysicalUnits': 50.0,
            'progress_bar': True,
        },
        maxit=1,
        verbose=2,
        force_2d=True,
    )
    print("Group-wise motion correction finished")

    os.makedirs(os.path.dirname(coreg_path), exist_ok=True)
    np.save(coreg_path, coreg)
    
    # Save the coregistered 4D image as NIfTI
    nii_path = coreg_path.replace('.npy', '.nii.gz')
    nii_fit = coreg_path.replace('.npy', '_fit.nii.gz')
    save_nifti(coreg, nii_path)
    save_nifti(fit, nii_fit)
               
    return coreg

# os.path.join(COREG_PATH, "dwi_coreg.npy"), coreg

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

def map_download(bvals_orig, bvecs_orig, raw_data, coreg_data, fa_path, md_path):

    tenfit_raw = dwi_model(raw_data, bvals_orig, bvecs_orig)
    tenfit_coreg = dwi_model(coreg_data, bvals_orig, bvecs_orig)

    FA_raw = fractional_anisotropy(tenfit_raw.evals)
    FA_coreg = fractional_anisotropy(tenfit_coreg.evals)
    FA_raw = np.clip(FA_raw, 0, 1)
    FA_coreg = np.clip(FA_coreg, 0, 1) 

    MD_raw = tenfit_raw.md
    MD_coreg = tenfit_coreg.md
    
    save_nifti(FA_coreg, fa_path)
    save_nifti(MD_coreg, md_path)
    save_nifti(FA_raw, fa_path.replace('.nii.gz', '_RAW.nii.gz'))
    save_nifti(MD_raw, md_path.replace('.nii.gz', '_RAW.nii.gz'))

def map_comparison(bvals_orig, bvecs_orig, raw_data, coreg_data, fa_path, md_path):

    if bvecs_orig.shape[0] == 3:
        bvecs_orig = bvecs_orig.T

    # Now check the dimensions of the actual data array
    if raw_data.ndim == 4:
        # Your logic for 4D volumes
        pass

    tenfit_raw = dwi_model(raw_data, bvals_orig, bvecs_orig)
    tenfit_coreg = dwi_model(coreg_data, bvals_orig, bvecs_orig)

    FA_raw = fractional_anisotropy(tenfit_raw.evals)
    FA_coreg = fractional_anisotropy(tenfit_coreg.evals)
    FA_raw = np.clip(FA_raw, 0, 1)
    FA_coreg = np.clip(FA_coreg, 0, 1) 

    MD_raw = tenfit_raw.md
    MD_coreg = tenfit_coreg.md

    FA_raw_napari = np.transpose(FA_raw, (2, 1, 0))
    MD_raw_napari = np.transpose(MD_raw, (2, 1, 0))
    FA_coreg_napari = np.transpose(FA_coreg, (2, 1, 0))
    MD_coreg_napari = np.transpose(MD_coreg, (2, 1, 0))

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

def image_download(raw_data, coreg_data, image_path):

    save_nifti(coreg_data, image_path)
    save_nifti(raw_data, image_path.replace('.nii.gz', '_RAW.nii.gz'))

def image_comparison(raw_data, coreg_data, image_path):

    viewer = napari.Viewer() #Order = (Slider 1, Slider 2, Y, X)

    if coreg_data.ndim == 4:
        # Standard 4D: Move Vol to front, then Z, Y, X
        t_axes = (3, 2, 1, 0)
    else:
        # Single slice 3D: Move Vol to front, then Y, X
        t_axes = (2, 1, 0)
        z=10
        raw_data = raw_data[:, :, z, :]
    
    raw_data = np.transpose(raw_data, t_axes)
    raw_data = np.flip(raw_data, axis=(-2))
    coreg_data = np.transpose(coreg_data, t_axes)
    coreg_data = np.flip(coreg_data, axis=(-2))
    
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


if __name__=='__main__':

# 1. DEFINE YOUR SPECIFIC PATHS HERE

    # PATHS for RAW DATA
    RAW_DICOM = r"C:\Users\eic20eh\Downloads\case-studies-main\data\iBE-2128-001_followup\scans\33-DTI_kidneys_cor_oblique_fb_DFC\resources\DICOM"
    BUILD = r"C:\Users\eic20eh\Downloads\case-studies-main\data\iBE-2128-001_conversion"
   
    # ARGUMENT SETUP  
    parser = argparse.ArgumentParser(description="Kidney DTI Pipeline")
    parser.add_argument('--dicom', type=str, default=RAW_DICOM)
    parser.add_argument('--build', type=str, default=BUILD)
    parser.add_argument('--dcm2niix', type=str, default="dcm2niix")
    args = parser.parse_args()

    # CONVERSION
    run_dcm2niix(args.dicom, args.build, args.dcm2niix)

    # PATHS TO INPUT FILES
    dti_path = os.path.join(args.build, "dti.nii.gz")
    bval_path = os.path.join(args.build, "dti.bval")
    bvec_path = os.path.join(args.build, "dti.bvec")
    
    # VARIABLES
    nii_obj = nib.load(dti_path)
    data = nii_obj.get_fdata()
    bvals = np.loadtxt(bval_path)
    bvecs = np.loadtxt(bvec_path)

    # PATHS TO OUTPUT FILES
    resultspath = os.path.join(args.build, "DTI_coreg.npy")
    imagepath = os.path.join(args.build, "DTI_coreg.nii.gz")
    fapath = os.path.join(args.build, "FA_map.nii.gz")
    mdpath = os.path.join(args.build, "MD_map.nii.gz")

    # MOTION CORRECTION
    motion_correction(data, resultspath, bvals, bvecs)

    # RESULTS
    coreg = np.load(resultspath)

    # CHECKS
    # view_unregistered_maps(bvals, bvecs, data)
    # view_registered_maps(bvals, bvecs, coreg)
    # view_unregistered_Images(data)
    # view_registered_Images(coreg)
    map_download(bvals, bvecs, data, coreg, fapath, mdpath)
    # map_comparison(bvals, bvecs, data, coreg, fapath, mdpath)
    image_download(data, coreg, imagepath)
    image_comparison(data, coreg, imagepath)
    # numerical_check(data, coreg)

#VIEW MAPS AND IMAGES
#.\.venv\Scripts\activate
#napari