# Overview
Image analysis workflow for quantification of CDK-activity in embryos and gastruloids.
Please refer to [Saykali, et al. Cell Report 2025](https://doi.org/10.1016/j.celrep.2025.115558)

# Image segementation and quantification
Image segmentation and quantifcation utilizes Python (v3.10). Custom [Cellpose](https://github.com/MouseLand/cellpose) models were used for segmentation. Segmentation using GPU-enabled HPC clusters is recommended, but is not essential.
For image visualization and verification, [napari](https://github.com/napari/napari) multi-dimensional image viewer was used, but other software such as [FIJI](https://fiji.sc/) can be substituted.

For live embryo analysis, cell tracks were assigned and edited using [Imaris Bitplane (v10.2)](https://imaris.oxinst.com/), however other software such as FIJI can be used as well.

## live_embryo_analysis.ipynb:

- Open live embryo timelapse .CZI image as a numpy array using aicsimageio.
  + Channel 1: CDK reporter
  + Channel 2: H2B-GFP
- Save each H2B nuclear image as individual timepoint
- Open fixed embryo .CZI image as a numpy array
- Use Cellpose GUI and trained segmentation model for 3D segmentation
  + Label mask is saved as h2b_{timepoint}_cp_mask.tif
- Stitch label mask timepoints into timeseries
- Import label mask with original data into Imaris
  + Use Imaris to create tracks and assign CDK2/SOX2 identity based on fixed image
- Import .IMS file with labels as numpy array using imaris_ims_file_reader
- Use label expand to create a cytoplasmic area for each nuclei
  + Cytoplasmic label is saved as cyto_{timepoint}_{z_plane}.tif
- Quantify image statistics for each nuclear and cytoplasmic region

# Data analysis and plotting
Data analysis scripts utilizes R (v4.2).
