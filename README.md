# Overview
Image analysis workflow for quantification of CDK-activity in embryos and gastruloids.
Please refer to [Saykali, et al. Cell Report 2025](https://doi.org/10.1016/j.celrep.2025.115558)

# Image segementation and quantification
Image segmentation and quantifcation utilizes Python (v3.10). Custom [Cellpose](https://github.com/MouseLand/cellpose) models were used for segmentation. Segmentation using GPU-enabled HPC clusters is recommended, but is not essential.
For image visualization and verification, [napari](https://github.com/napari/napari) multi-dimensional image viewer was used, but other software such as [FIJI](https://fiji.sc/) can be substituted.

For live embryo analysis, cell tracks were assigned and edited using [Imaris Bitplane (v10.2)](https://imaris.oxinst.com/), however other software such as FIJI can be used as well.

Please see the attached 'example_images.zip' file for test images and expected segmentation results.

# Data analysis and plotting
Data analysis scripts utilizes R (v4.2).
