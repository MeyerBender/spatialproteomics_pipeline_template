import matplotlib.pyplot as plt
import pandas as pd
import spatialproteomics as sp
import tifffile
import yaml

image_path = snakemake.input.image_path
marker_path = snakemake.input.marker_path
segmentation_path = snakemake.input.mask_path
parameter_config_path = snakemake.input.parameter_config_path

# === markers ===
markers = list(pd.read_csv(marker_path, header=None)[0])

# === data loading ===
image = tifffile.imread(image_path)
img = sp.load_image_data(image, channel_coords=markers)

# === parameter config ===
with open(parameter_config_path, "r") as yaml_file:
    yaml_data = yaml.safe_load(yaml_file)
mask_growth = yaml_data["mask_growth_pixels"]
min_cell_size = yaml_data["min_cell_size"]
max_cell_size = yaml_data["max_cell_size"]

# === reading segmentation ===
segmentation = tifffile.imread(segmentation_path)
img = img.pp.add_segmentation(segmentation)

# === filtering segmentation ===
img = img.pp.add_observations("area").pp.filter_by_obs(
    "area", func=lambda x: (x >= min_cell_size) & (x <= max_cell_size)
)

# === growing segmentation ===
img = img.pp.grow_cells(iterations=mask_growth)

# === removing outlier masks ===
img = img.pp.remove_outlying_cells()

# === plotting ===
plt.figure(figsize=(10, 10))
_ = img.pp["DAPI"].pl.colorize("gold").pl.show(render_segmentation=True)
plt.savefig(snakemake.output.plot_path, bbox_inches="tight", pad_inches=0)

# === exporting ===
segmentation_processed = img["_segmentation"].values
tifffile.imwrite(snakemake.output.mask_path, segmentation_processed)
