import nibabel as nib
import glob
import os

# Path to the NIfTI files created by dcm2niix in your previous run
results_path = "/mnt/parscratch/users/eic20eh/data/ibeat_diff/mc_results/series_*/*.nii.gz"

print(f"{'Folder':<20} | {'Shape':<20} | {'Status'}")
print("-" * 60)

for nii in glob.glob(results_path):
    try:
        header = nib.load(nii).header
        shape = header.get_data_shape()
        folder = nii.split('/')[-2]
        
        # DTI must be 4D (X, Y, Z, Directions)
        if len(shape) < 4:
            status = "❌ 3D STACK (SKIP)"
        else:
            status = "✅ VALID 4D"
            
        print(f"{folder:<20} | {str(shape):<20} | {status}")
    except Exception as e:
        pass