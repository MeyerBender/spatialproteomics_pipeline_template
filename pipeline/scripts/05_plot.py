import matplotlib
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import spatialproteomics as sp
import xarray as xr
import yaml

# === DATA LOADING ===
# Load the YAML data
with open(snakemake.input.yaml_marker_path, "r") as f:
    config = yaml.safe_load(f)

cell_types = config["cell_type"]

# Generate a color palette
color_map = matplotlib.colormaps.get_cmap("tab10")

# Initialize dictionaries
cell_type_to_color = {}
marker_to_color = {}
marker_to_cell_type = {}

# Assign colors and relationships
for i, (cell_type, markers) in enumerate(cell_types.items()):
    color = color_map(i)
    color_hex = mcolors.to_hex(color)

    cell_type_to_color[cell_type] = color_hex
    for marker in markers:
        marker_to_color[marker] = color_hex
        marker_to_cell_type[marker] = cell_type

ds = xr.open_zarr(snakemake.input.zarr_path)
ds = ds.la.set_label_colors(cell_type_to_color.keys(), cell_type_to_color.values())

# === PLOTTING OF MARKERS AND ASSOCIATED CELL TYPES ===
nrows = len(marker_to_cell_type) + 1
ncols = 2
scaling = 5

fig, ax = plt.subplots(nrows, ncols, figsize=(ncols * scaling, nrows * scaling))
ax = ax.flatten()

# first row
ds.pp[marker_to_color.keys()].pl.colorize(marker_to_color.values()).pl.show(ax=ax[0])
ds.pl.show(render_image=False, render_labels=True, ax=ax[1])
ax[0].set_title("Markers")
ax[1].set_title("Cell Types")

# subsequent rows
i = 2
for marker, color in marker_to_color.items():
    celltype = marker_to_cell_type[marker]
    ds.pp[marker].pl.colorize(color).pl.show(ax=ax[i])
    ax[i].set_title(marker)
    i += 1
    try:
        ds.la[celltype].pl.show(render_image=False, render_labels=True, ax=ax[i])
    except ValueError:
        # Plot an empty black image
        ax[i].imshow(np.zeros((100, 100)), cmap="gray", vmin=0, vmax=1)
    ax[i].set_title(celltype)
    i += 1

for axis in ax:
    axis.axis("off")

plt.tight_layout()
plt.savefig(snakemake.output.plot_path_1, bbox_inches="tight", pad_inches=0)

# === PLOTTING OF HEATMAP ===
adata = ds.tl.convert_to_anndata()

sc.pl.heatmap(
    adata,
    var_names=list(marker_to_cell_type.keys()),  # or a subset of genes/markers
    groupby="_labels",  # this colors the rows by label
    use_raw=False,  # or True if your data is in adata.raw
)
plt.savefig(snakemake.output.plot_path_2, bbox_inches="tight", pad_inches=0)

# === PLOTTING OF UMAPS ===
# Optional: recompute UMAP if not yet done
if "X_umap" not in adata.obsm:
    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)

# First, plot UMAP colored by labels
sc.pl.umap(
    adata,
    color="_labels",
    palette=adata.uns.get("_labels_colors", None),
    title="UMAP by Cell Type",
    frameon=False,
)
plt.savefig(snakemake.output.plot_path_3, bbox_inches="tight", pad_inches=0)

# Now, plot UMAPs colored by each marker
markers = adata.var_names

# Determine subplot layout (e.g. 3 per row)
n_cols = 10
n_rows = (len(markers) + n_cols - 1) // n_cols

# Create figure
fig, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 5 * n_rows))

# Flatten axes for easy iteration
axes = axes.flatten()

for i, marker in enumerate(markers):
    sc.pl.umap(adata, color=marker, ax=axes[i], title=marker, frameon=False, show=False)

# Remove any unused subplots
for j in range(len(markers), len(axes)):
    fig.delaxes(axes[j])

plt.tight_layout()
plt.savefig(snakemake.output.plot_path_4, bbox_inches="tight", pad_inches=0)
