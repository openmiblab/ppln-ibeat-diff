import torch
import matplotlib.pyplot as plt
import os
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"
import json
import argparse
from monai.networks.nets import LocalNet
from monai.networks.blocks import Warp
from monai.losses import BendingEnergyLoss, MultiScaleLoss
from monai.data import DataLoader, Dataset, CacheDataset, decollate_batch
from monai.transforms import (
    Compose, LoadImaged, EnsureChannelFirstd, 
    Resized, ScaleIntensityd, ToTensord, SaveImaged, MapTransform,
)

parser = argparse.ArgumentParser()
parser.add_argument("--json", type=str, required=True)
parser.add_argument("--out", type=str, required=True)
parser.add_argument("--epochs", type=int, default=20)
args = parser.parse_args()

# Use these in your script
with open(args.json, "r") as f:
    config = json.load(f)
output_root = args.out

class ExtractSliced(MapTransform):
    def __call__(self, data):
        d = dict(data)
        for key in self.keys:
            z_idx = int(d["slice_index"])
            sliced = d[key][..., z_idx]
            d[key] = sliced.permute(1, 2, 0).unsqueeze(0)  # (X, Y, Z, T) -> (T, X, Y)
        return d

# 1. Load your newly created JSON
with open(args.json, "r") as f:
    config = json.load(f)

output_root = args.out
os.makedirs(output_root, exist_ok=True)

val_saver = SaveImaged(
    keys=["moving", "fixed", "warped"], 
    output_dir=os.path.join(output_root, "val_results"),
    output_postfix="val", 
    resample=False
)

test_saver = SaveImaged(
    keys=["moving", "fixed", "warped"], 
    output_dir=os.path.join(output_root, "test_results"),
    output_postfix="test", 
    resample=False
)

# 2. Define Transforms (Resizing to 128x128x64 to save memory)
train_transforms = Compose([
LoadImaged(keys=["moving", "fixed"]),
    EnsureChannelFirstd(keys=["moving", "fixed"]),
    # Use 'bilinear' for images to keep them smooth
    ExtractSliced(keys=["moving", "fixed"]), 
    Resized(keys=["moving", "fixed"], spatial_size=(160, 160, 144), mode="bilinear"),
    # Percentile scaling handles outliers (bright fat/artifacts) better than min-max
    ScaleIntensityd(keys=["moving", "fixed"], minv=0.0, maxv=1.0), 
    ToTensord(keys=["moving", "fixed"]),
])

# 3. Data Loaders
train_ds = CacheDataset(
    data=config["training"], 
    transform=train_transforms, 
    cache_rate=1.0, 
    num_workers=8
)
train_loader = DataLoader(train_ds, batch_size=16, shuffle=True)

val_ds = CacheDataset(data=config["testing"], transform=train_transforms, cache_rate=1.0)
val_loader = DataLoader(val_ds, batch_size=1)

test_ds = CacheDataset(data=config["testing"], transform=train_transforms, cache_rate=1.0)
test_loader = DataLoader(test_ds, batch_size=1)

# 4. Model, Loss, and Optimizer
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
model = LocalNet(
    spatial_dims=3,
    in_channels=2,   # Moving + Fixed concatenated
    out_channels=3, 
    num_channel_initial=16,  # Predicting DVF (dx, dy, dz)
    extract_levels=[0,1,2,3]
).to(device)

warp_layer = Warp().to(device)
optimizer = torch.optim.Adam(model.parameters(), 1e-4)
scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
    optimizer, 
    mode='min',      # 'min' because we want to minimize MSE loss
    factor=0.5,      # Reduce LR by half (0.0001 -> 0.00005)
    patience=5,      # Wait 5 epochs with no improvement before dropping LR
    verbose=True     # Print a message when the LR is updated
)

# MultiScaleLoss compares the warped image to the fixed image at different resolutions
image_loss = MultiScaleLoss(torch.nn.MSELoss())
# BendingEnergyLoss ensures the deformation field is smooth/realistic
regularization_loss = BendingEnergyLoss()

num_epochs = args.epochs
train_loss_history = []
val_mse_history = []

# 5. Training Loop
for epoch in range(num_epochs):
    model.train()
    epoch_loss = 0
    for batch_data in train_loader:
        moving = batch_data["moving"].to(device)
        fixed = batch_data["fixed"].to(device)

        optimizer.zero_grad()
        # Predict Displacement Field
        ddf = model(torch.cat((moving, fixed), dim=1))
        # Warp the moving image using the predicted field
        warped_image = warp_layer(moving, ddf)
        
        # Calculate Loss
        loss = image_loss(warped_image, fixed) + 0.005 * regularization_loss(ddf)
        loss.backward()
        optimizer.step()
        epoch_loss += loss.item()

    avg_train_loss = epoch_loss / len(train_loader)
    train_loss_history.append(avg_train_loss)
    print(f"Epoch {epoch} - Avg Train Loss: {epoch_loss/len(train_loader):.4f}")

    # Run validation every 5 epochs to save time
    if (epoch + 1) % 1 == 0:
        model.eval()
        val_loss = 0
        os.makedirs(output_root, exist_ok=True) 
        with torch.no_grad():
            for i, val_data in enumerate(val_loader):
                v_moving = val_data["moving"].to(device)
                v_fixed = val_data["fixed"].to(device)
                
                # Predict
                v_ddf = model(torch.cat((v_moving, v_fixed), dim=1))
                v_warped = warp_layer(v_moving, v_ddf)
                
                # Calculate metric (MSE)
                v_loss = torch.nn.functional.mse_loss(v_warped, v_fixed)
                val_loss += v_loss.item()

                if (epoch +1) == num_epochs and i < 12:  # Save outputs for the last epoch only

                    val_saver.output_dir = os.path.join(output_root, "val_results")
                    os.makedirs(val_saver.output_dir, exist_ok=True)

                    val_data["warped"] = v_warped.cpu()
                    val_data["moving"] = v_moving.cpu()
                    val_data["fixed"] = v_fixed.cpu() 

                    val_outputs = decollate_batch(val_data)
                    item = val_outputs[0] # Get the first (only) item in the batch    
                    unique_id = item["id"]
                    
                    for k in ["moving", "fixed", "warped"]:
                        # 1. Ensure warped has metadata (copy from fixed)
                        if k == "warped":
                            item["warped"].copy_meta_from(item["fixed"])
                        
                        item[k].meta["filename_or_obj"] = f"{unique_id}_{k}_val.nii.gz"

                    val_saver(item)

        avg_val_mse = val_loss / len(val_loader)
        val_mse_history.append(avg_val_mse)

        scheduler.step(avg_val_mse)  # Adjust learning rate based on validation MSE
                
        print(f"Epoch {epoch} - Avg Valid Loss: {avg_val_mse:.4f}")

torch.save(model.state_dict(), os.path.join(output_root, "localnet_reg_v2.pth"))

# Generate Loss Curves
plt.figure(figsize=(10, 6))
plt.plot(range(len(train_loss_history)), train_loss_history, label='Training Loss (MSE + Reg)')
plt.plot(range(len(val_mse_history)), val_mse_history, label='Validation MSE')
plt.xlabel('Epochs')
plt.ylabel('Loss/MSE')
plt.title('Training and Validation Performance')
plt.legend()
plt.grid(True)

# Save the plot
plot_path = os.path.join(output_root, "training_curves.png")
plt.savefig(plot_path)
print(f"Graph saved to: {plot_path}")

# --- Final Inference on Test Set ---
print("Starting Inference on Test Set...")
model.eval()

test_saver.output_dir = os.path.join(output_root, "test_results")
os.makedirs(test_saver.output_dir, exist_ok=True)

# Update saver to put test results in a specific folder
with torch.no_grad():
    for i, test_data in enumerate(test_loader):
        t_moving = test_data["moving"].to(device)
        t_fixed = test_data["fixed"].to(device)
        
        t_ddf = model(torch.cat((t_moving, t_fixed), dim=1))
        t_warped = warp_layer(t_moving, t_ddf)
        
        # Prepare for saving
        test_data["warped"] = t_warped.cpu()
        test_data["moving"] = t_moving.cpu()
        test_data["fixed"] = t_fixed.cpu()

        # Decollate and fix naming
        test_outputs = decollate_batch(test_data)
        item = test_outputs[0]

        unique_id = item["id"]

        for k in ["moving", "fixed", "warped"]:
            if k == "warped":
                item["warped"].copy_meta_from(item["fixed"])
            
            item[k].meta["filename_or_obj"] = f"{unique_id}_{k}_test.nii.gz"

        test_saver(item)
        print(f"Saved results for test case {unique_id}")