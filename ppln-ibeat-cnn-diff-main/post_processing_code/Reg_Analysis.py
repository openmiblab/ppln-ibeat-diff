import os
import nibabel as nib
import matplotlib.pyplot as plt
import numpy as np

def generate_kymograph(file_paths, labels, save_name, slice_idx=86):
    fig, axes = plt.subplots(len(file_paths), 1, figsize=(15, 5 * len(file_paths)))
    
    # Check if we have multiple files or just one to avoid indexing errors
    if len(file_paths) == 1:
        axes = [axes]

    for idx, (path, label) in enumerate(zip(file_paths, labels)):
        if not os.path.exists(path):
            print(f"ERROR: File not found at {path}")
            continue
            
        img = nib.load(path)
        data = img.get_fdata()
        
        # Slicing logic for (X, Y, T) & (X, Y, Z, T)
        # Slicing logic for (85, 172, 172) & (145, 29, 172, 172)
        if data.ndim == 3:
            data = data.transpose (1, 0, 2)
            data = np.flip(data, axis=0)
            kymo = data[:, 46, :] 
        elif data.ndim == 4:
            if data.shape[0] < 172: # warped
                data = data.transpose (1, 0, 2, 3)
                print(f"{data.shape}")
                data = np.flip(data, axis=0)
                kymo = data[:, 50, 23, :]
            else: # registered & unregistered
                data = data.transpose (1, 0, 2, 3)
                print(f"{data.shape}")
                data = np.flip(data, axis=0)
                kymo = data[: ,56, 23, :] 
        else:
            print(f"Unexpected dimensions: {data.ndim}")
            continue
            
        print(f"Processing {label}: {data.shape} -> {kymo.shape}")
        
        axes[idx].imshow(kymo, cmap='gray', aspect='auto', origin='lower')
        axes[idx].set_title(f"{label}", fontsize=20, fontweight='bold', pad = 15)
        axes[idx].set_ylabel("Spatial Position (X)", fontsize = 18, labelpad=10)
        axes[idx].set_xlabel("Time/Slice Index", fontsize = 18, labelpad=10)
            
    plt.tight_layout()
    
    # USE ABSOLUTE PATH TO BE CERTAIN
    full_save_path = os.path.abspath(save_name)
    plt.savefig(full_save_path)
    plt.close(fig) # Releases memory and finalizes the file
    
    print("-" * 30)
    print(f"SUCCESS! Image saved at:\n{full_save_path}")
    print("-" * 30)

# --- YOUR PATHS ---
base = r"C:\Users\eic20eh\Downloads\testing_outputs\test_results\iBE-2128-037_niis" 
files = [
    os.path.join(base, "dwi_unregistered.nii"),
    os.path.join(base, "dwi_registered.nii"),
    os.path.join(base, "iBE-2128-037_baseline_FULL_WARPED.nii")
]

save_path = os.path.join(base, "Motion_Comparison.png")

generate_kymograph(files, ['UN-REG', 'OP-REG', 'DL-REG'], save_path)