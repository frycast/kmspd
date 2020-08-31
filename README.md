## kmspd: k-means on SPD matrices in R

This repository holds code and data for the paper [k-means on Positive Definite Matrices, and an Application to Clustering in Radar Image Sequences](https://arxiv.org/pdf/2008.03454.pdf) by [Daniel Fryer](www.danielvfryer.com), [Hien Nguyen](www.hiendn.github.io) and [Pascal Castellazzi](https://www.linkedin.com/in/pascal-castellazzi-b9533ba0/).

### Features

The script [R/SAR_application.R](R/SAR_application.R) calls functions that save their intermediary results in Rds format. These functions check for previously saved results before executing. So, most computationally intensive tasks need not be repeated more than once -- making for easier / quicker interactive exploratory analysis! 

### Setup

After cloning this repo, and before running any scripts, download the files from the [Zenodo archive](https://zenodo.org/record/4008883) and place them in [SAR_app_data/MG](SAR_app_data/MG) without changing their file names. Download links:

```R
"https://zenodo.org/record/4008883/files/CC_sub_norm.tif?download=1"
"https://zenodo.org/record/4008883/files/VH_sub_norm.tif?download=1"
"https://zenodo.org/record/4008883/files/VV_sub_norm.tif?download=1"
```

### Structure

* [R/SAR_application.R](R/SAR_application.R) is the script that produces all results in the paper.

* [SAR_app_data/MG](SAR_app_data/MG) contains pre-processed SAR data.
* [results](results) contains figures and other results produced in the code.

### Tools

[`{rsar}`](https://github.com/frycast/rsar) :package: provides tools for working with radar images in `SAR_matrix` format.

[`{tidyverse}`](https://www.tidyverse.org/) :package: for data manipulation and publication graphics.

[`{raster}`](https://cran.r-project.org/web/packages/raster/index.html) :package: for working with raster images.

[`{plotly}`](https://github.com/plotly/) :package: for interactive 3D scatter plots!

### Data usage license

The authors must be acknowledged in any publications or applications that make use of the [data](SAR_app_data/MG) in this repository.  





