import matplotlib.pyplot as plt
import pandas as pd
import spatialproteomics as sp
import tifffile

image_path = snakemake.input.image_path
marker_path = snakemake.input.marker_path

# === markers ===
markers = list(pd.read_csv(marker_path, header=None)[0])

# === data loading ===
image = tifffile.imread(image_path)
img = sp.load_image_data(image, channel_coords=markers)

# === segmentation ===
img = img.tl.cellpose(channel="DAPI")

# === plotting ===
plt.figure(figsize=(10, 10))
_ = img.pp["DAPI"].pl.colorize("gold").pl.show(render_segmentation=True)
plt.savefig(snakemake.output.plot_path, bbox_inches="tight", pad_inches=0)

# === exporting ===
segmentation_raw = img["_segmentation"].values
tifffile.imwrite(snakemake.output.mask_path, segmentation_raw)
