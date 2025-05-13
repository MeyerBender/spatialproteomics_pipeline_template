import matplotlib.pyplot as plt
import pandas as pd
import spatialproteomics as sp
import tifffile
import yaml

image_path = snakemake.input.image_path
marker_path = snakemake.input.marker_path
segmentation_path = snakemake.input.mask_path

# === markers ===
markers = list(pd.read_csv(marker_path, header=None)[0])

# === data loading ===
image = tifffile.imread(image_path)
img = sp.load_image_data(image, channel_coords=markers)

# === reading segmentation ===
segmentation = tifffile.imread(segmentation_path)
img = img.pp.add_segmentation(segmentation)

# === quantifying protein expression ===
img = img.pp.add_quantification(func="intensity_mean").pp.transform_expression_matrix(
    method="arcsinh"
)

# === running astir ===
with open(snakemake.input.yaml_marker_path, "r") as stream:
    marker_dict = yaml.safe_load(stream)

img = img.tl.astir(marker_dict)

img.drop_encoding().to_zarr(snakemake.output.zarr_path)
