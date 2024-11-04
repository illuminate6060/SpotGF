# SpotGF
SpotGF: Denoising Spatially Resolved Transcriptomics Data Using an Optimal Transport-Based Gene Filtering Algorithm
https://doi.org/10.1016/j.cels.2024.09.005

preprint paper:Optimal Transport Method-Based Gene Filter (GF) Denoising Algorithm for Enhancing Spatially Resolved Transcriptomics Data

![GRAPHICAL ABSTRACT_revised](https://github.com/user-attachments/assets/d1e78bc1-85d0-4e06-85be-28da83cbb0e2)


## Installation
We recommend using conda to manage the installation of all dependencies. To do this, simply run:

```
conda create --name SpotGF
conda activate SpotGF

# conda config --add channels conda-forge ##If you do not have this channel added before#
conda install python=3.7 pandas pot=0.8.2 numpy scipy matplotlib descartes scanpy=1.9.2 seaborn
```

Then, download this repo and install it.
```
git clone [repo_path]
cd [path/to/SpotGF]
pip install .
```

The total installation time is around 5 mintunes. If error occuors, please upgrade pip and try again.


## Usage
Once the input data have been processed into the supported format, the full SpotGF workflow can be run by calling the `SpotGF.py` script. The input files can include various formats such as `gem`, `txt`, `csv`, `gem.gz`and others, containing the raw SRT data information. These files must contain specific columns including `geneID`, `x`, `y`, `MIDCount`. 

The denoising resolution can be adjusted using the `binsize` parameter. A smaller `binsize` will result in a finer denoising effect but will also increase the processing time.Please note that the `binsize` parameter is only used when calculating SpotGF scores and it does not affect the final denoised data. 

If the input file contains columns named ‘cen_x’ and ‘cen_y’, the denoising process can be performed based on the cell bin data, achieving denoising at the single-cell level by setting ‘binsize = 1’ or ‘-b 1’.

The `proportion` parameter determines the proportion of genes to be retained in the final denoised data. For example, setting proportion=0.9 will retain only 90% of the effective genes, resulting in a new denoised SRT dataset. 


```
#python [path/to/SpotGF.py] -i [path/to/input.gem] -b [binsize] -p [proportion] -s [spotsize]

cd /path/test

python ../SpotGF.py -i ./demo.gem -o ./result -b 70 -p 0.5 -s 5
```

If you use Spot in Jupyter environment, you can choose blow Usage.

```
import SpotGF	

gem_path = './demo.gem'  #please change this file path base on your work directory

binsize =70
lower= 0
upper = 100000
proportion  = 0.1
max_iterations = 10000
auto_threshold = True
visualize=True
spot_size=5
alpha=0

#initial class
spotgf =SpotGF.SpotGF(gem_path,binsize,proportion,auto_threshold,lower,upper,max_iterations,output,visualize,spot_size,alpha)

#acquire SpotGF socres distribution figure to help choose proportion 
GF_df = spotgf.calculate_GFscore(gem_path,binsize,alpha)

#denoised SRT data
new_gem  = spotgf.generate_GFgem(gem_path,GF_df,proportion,auto_threshold,visualize,spot_size)
```


### Data preparation
SpotGF workflow requires two mandatory files, a `input.gem` (with "\t" or "," as the delimiter) for gene expression data,

|geneID|x|y|MIDCount|
|-----|-----|-----|-----|
|#gene1|23|36|  1|
|#gene1|24|35|  1|
|#gene1|23|30|  1|
|#gene2|20|31|  1|
|#gene2|21|22|  1|


### Output files
"SpotGF_scores.txt": This file contains the SpotGF scores for each gene. The SpotGF score indicates the degree of clustering or diffusion of a gene. A smaller SpotGF score suggests that the gene is more diffuse, while a larger SpotGF score indicates that the gene is more clustered.

"SpotGF_proportion_0.5.gem": A gene expression matrix file with noise reduction that retains only 50% of the genes.

"SpotGF_auto_threshold.gem": Gene expression matrix file after denoising based on an automatic filtering threshold method.

"alpha_shape.png": It shows the tissue contour used by SpotGF during the denoising process, with special attention to the fact that a highly detailed outline is not required; only a rough outline is needed.

"Spatial_automatic.png": Spatial expression levels of data after denoising based on an automatic filtering threshold method.

"Spatial_proportion.png": Spatial expression levels of data after denoising based on a manual threshold method.

"Violinplot_n_genes_by_counts.png": Violin plot distribution of n_genes_by_counts across several datasets.

"Violinplot_total_counts.png": Violin plot distribution of total_counts across several datasets.


### List of Parameters
```
positional arguments:

  -i                Input SRT data files path.

  -o                Outpath for saving results.
    
  -b                Denoising resolution binsize, must int type, default=10.

  -lower            Lower limit for tissue structures capturing optimization, default=0.

  -upper            Upper limit for tissue structures capturing optimization', default=sys.float_info.max.

  -max_iterations   maximum number of iterations when capturing tissue structures', default=10000.

  -p                Proportion of gene numbers, must float type [0,1], default=0.5.

  -auto_threshold   if True, return a denoised gem file using automatic threshold, default=True.

  -v                Visualize SpotGF-denoised data, default=True.

  -s                Spot size used for spatial expression figure, default=5.

  -a                Alpha for tissue boundary detection, default use auto optimizealpha, default=0.
```

### Contact Us
If you have any suggestions/ideas for SpotGF or are having issues trying to use it, please don't hesitate to reach out to us.

Lin Du, dulin@genomics.cn 


### Cite
Du, L., Kang, J., Hou, Y., Sun, H. X., & Zhang, B. (2024). SpotGF: Denoising spatially resolved transcriptomics data using an optimal transport-based gene filtering algorithm. Cell systems, 15(10), 969–981.e6. https://doi.org/10.1016/j.cels.2024.09.005
