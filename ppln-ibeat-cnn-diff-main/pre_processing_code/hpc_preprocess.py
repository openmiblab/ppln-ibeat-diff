import os
import glob
import numpy as np
import shutil
import json
import csv
import random
import nibabel as nib  # Added to check slice counts

def run_pipeline(sourcedir, outputdir):
    # Setup Directories
    train_dir = os.path.join(outputdir, "Train")
    val_dir = os.path.join(outputdir, "Val")    # New Validation Dir
    test_dir = os.path.join(outputdir, "Test")
    
    for d in [train_dir, val_dir, test_dir]:
        os.makedirs(d, exist_ok=True)
    
    report_path = os.path.join(outputdir, "preprocessing_report.csv")
    all_results = []
    passed_pairs = []

    subjects = [f for f in os.listdir(sourcedir) if os.path.isdir(os.path.join(sourcedir, f))]
    print(f"Analyzing {len(subjects)} subjects...")

    for folder in subjects:
        folder_path = os.path.join(sourcedir, folder)
        try:
            reg_files = glob.glob(os.path.join(folder_path, "**/*_coreg.npy"), recursive=True)
            if not reg_files: continue
            
            data = np.load(reg_files[0])
            # Check a central slice across all timepoints
            z_index = data.shape[2] // 2
            slice_data = data[:, :, z_index, :] # Shape is (Height, Width, Volumes)

            num_volumes = slice_data.shape[-1]
            volume_scores = []

            # 2. Use the first volume as the 'Fixed' reference
            fixed_vol = slice_data[:, :, 0].flatten()

            # 3. Compare every other volume to that reference
            for v in range(1, num_volumes):
                moving_vol = slice_data[:, :, v].flatten()
                
                # Calculate correlation
                correlation = np.corrcoef(fixed_vol, moving_vol)[0, 1]
                volume_scores.append(correlation)

            # 4. The final score is the average stability of all volumes at that slice
            score = np.mean(volume_scores)

            # Log it
            status = "PASSED" if score > 0.50 else "FAIL"
            
            unreg_nii = glob.glob(os.path.join(folder_path, "**/dwi_unregistered.nii.gz"), recursive=True)
            reg_nii = glob.glob(os.path.join(folder_path, "**/dwi_registered.nii.gz"), recursive=True)

            all_results.append([folder, f"{score:.4f}", status])

            if status == "PASSED" and unreg_nii and reg_nii:
                # Load header to find out how many slices (Z) we have
                img = nib.load(reg_nii[0])
                num_slices = img.shape[2] # Assuming (X, Y, Z, T)
                
                passed_pairs.append({
                    "id": folder,
                    "moving": unreg_nii[0],
                    "fixed": reg_nii[0],
                    "src_path": folder_path,
                    "slices": num_slices
                })

        except Exception as e:
            print(f"Error skipping {folder}: {e}")

    # 3. 70:15:15 Split by SUBJECT
    random.seed(42)
    random.shuffle(passed_pairs)
    
    num_subs = len(passed_pairs)
    train_end = int(num_subs * 0.70)
    val_end = int(num_subs * 0.85)
    
    final_json_data = {"training": [], "validation": [], "testing": []}

    for i, pair in enumerate(passed_pairs):
        # Determine split
        if i < train_end:
            target_root, key = train_dir, "training"
        elif i < val_end:
            target_root, key = val_dir, "validation"
        else:
            target_root, key = test_dir, "testing"
            
        dest_path = os.path.join(target_root, pair["id"])
        shutil.copytree(pair["src_path"], dest_path, dirs_exist_ok=True)
        
        new_moving = glob.glob(os.path.join(dest_path, "**/dwi_unregistered.nii.gz"), recursive=True)[0]
        new_fixed = glob.glob(os.path.join(dest_path, "**/dwi_registered.nii.gz"), recursive=True)[0]
        
        # Create an entry for EVERY SLICE
        for s in range(pair["slices"]):
            final_json_data[key].append({
                "id": f"{pair['id']}_slice{s}",
                "moving": new_moving,
                "fixed": new_fixed,
                "slice_index": s  # Your training script will use this to slice the 4D file
            })

    # 4. Save Outputs
    with open(report_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["Patient_ID", "Score", "Status"])
        writer.writerows(all_results)

    with open(os.path.join(outputdir, "dataset_monai.json"), 'w') as f:
        json.dump(final_json_data, f, indent=4)

    print(f"Complete. Train slices: {len(final_json_data['training'])} | Val: {len(final_json_data['validation'])} | Test: {len(final_json_data['testing'])}")

if __name__ == "__main__":
    SOURCE = "/mnt/parscratch/users/eic20eh/data/ibeat_diff/mc_results"
    OUTPUT = "/mnt/parscratch/users/eic20eh/data/ibeat_diff/clean_results"
    run_pipeline(SOURCE, OUTPUT)