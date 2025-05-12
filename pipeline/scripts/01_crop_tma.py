import pandas as pd
import tifffile
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import NullLocator
import yaml


def plot_overview(dapi_img, metadata):
    # now we take the centroids, create crops with a size of crop_size around them and plot the resulting crops on top of the big image
    downsample = 4
    y_dim, x_dim= dapi_img.shape

    plt.figure(figsize=(8, 24))
    ax = plt.gca()
    ax.imshow(dapi_img[::downsample, ::downsample], extent=[0, x_dim, y_dim, 0])

    w, h = crop_size, crop_size

    for i, row in metadata.iterrows():
        x, y = row.loc["Centroid X µm"] * magnification, row.loc["Centroid Y µm"] * magnification
        sample_id = row["Name"]

        ax.vlines(x - w/2, y - h/2, y + h/2)
        ax.vlines(x + w/2, y - h/2, y + h/2)
        ax.hlines(y - h/2, x - w/2, x + w/2)
        ax.hlines(y + h/2, x - w/2, x + w/2)
        ax.text(x, y, sample_id,
                fontsize=8, color="deepskyblue", ha='center', va='center')

    plt.tight_layout()
    plt.savefig(snakemake.output.overview_plot, bbox_inches='tight', pad_inches=0)

    
def safe_crop(image, center_x, center_y, crop_size):
    """
    Crops a square region of size `crop_size` around (center_x, center_y) from an image.
    Supports both 2D and multi-dimensional arrays with spatial axes at the end.
    """
    half_crop = crop_size // 2
    *prefix_shape, y_dim, x_dim = image.shape

    x_min = int(center_x - half_crop)
    y_min = int(center_y - half_crop)
    x_max = x_min + crop_size
    y_max = y_min + crop_size

    # Adjust crop to stay within bounds
    if x_min < 0:
        x_max += -x_min
        x_min = 0
    if y_min < 0:
        y_max += -y_min
        y_min = 0
    if x_max > x_dim:
        x_min -= (x_max - x_dim)
        x_max = x_dim
    if y_max > y_dim:
        y_min -= (y_max - y_dim)
        y_max = y_dim

    # Clamp again
    x_min = max(x_min, 0)
    y_min = max(y_min, 0)
    x_max = min(x_max, x_dim)
    y_max = min(y_max, y_dim)

    if image.ndim == 2:
        return image[y_min:y_max, x_min:x_max]
    else:
        # Crop along last two axes (Y, X), preserve leading dims (e.g., channel)
        return image[..., y_min:y_max, x_min:x_max]


    
def plot_images(dapi_img, metadata):
    num_images = metadata.shape[0]
    num_rows = (num_images + 9) // 10  # Ensure at least 1 row

    fig, axes = plt.subplots(num_rows, 10, figsize=(20, num_rows * 2))
    axes = axes.flatten()

    for i, row in metadata.iterrows():
        # Load image (only looking at the DAPI staining)
        x, y = row.loc["Centroid X µm"] * magnification, row.loc["Centroid Y µm"] * magnification
        image = safe_crop(dapi_img, x, y, crop_size)
        sample_id = row["Name"]
        
        # Normalize the image to the 5th and 95th percentiles
        min_val = np.percentile(image, 5)
        max_val = np.percentile(image, 95)
        normalized_image = np.clip((image - min_val) / (max_val - min_val), 0, 1)

        # Plot the image
        axes[i].imshow(normalized_image, cmap='gray')
        axes[i].set_title(sample_id)
        axes[i].xaxis.set_major_locator(NullLocator())
        axes[i].yaxis.set_major_locator(NullLocator())

    # Remove empty subplots if necessary
    for i in range(num_images, num_rows * 10):
        fig.delaxes(axes[i])

    plt.tight_layout()
    plt.savefig(snakemake.output.individual_cores_plot, bbox_inches='tight', pad_inches=0)

    
def save_images(img, metadata):
    for i, row in metadata.iterrows():
        # Load image
        x, y = row.loc["Centroid X µm"] * magnification, row.loc["Centroid Y µm"] * magnification
        sample_id = row["Name"]
        
        # crop and store the sample
        img_cropped = safe_crop(img, x, y, crop_size)
        save_path = f"{snakemake.input.cropped_image_path}/{sample_id}/01_image.tiff"
        tifffile.imwrite(save_path, img_cropped)
       
    
# === parameter config ===
with open(snakemake.input.parameter_config_path, "r") as yaml_file:
    yaml_data = yaml.safe_load(yaml_file)
crop_size = yaml_data["crop_size"]
magnification = yaml_data["magnification"]

# === file reading ===
img = tifffile.imread(snakemake.input.tma_image_path)
dapi_img = img[0, :, :]

metadata = pd.read_csv(snakemake.input.tma_annotation_path, sep="\t")

# === plotting ===
plot_overview(dapi_img, metadata)
plot_images(dapi_img, metadata)

# === storing results ===
save_images(img, metadata)
