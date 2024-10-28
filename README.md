# SpotGF
SpotGF: Denoising Spatially Resolved Transcriptomics Data Using an Optimal Transport-Based Gene Filtering Algorithm
https://doi.org/10.1016/j.cels.2024.09.005

Preprint paper:Optimal Transport Method-Based Gene Filter (GF) Denoising Algorithm for Enhancing Spatially Resolved Transcriptomics Data


![GRAPHICAL ABSTRACT_revised](https://github.com/user-attachments/assets/07108fa9-3784-4aa9-a2dd-1d23041a88dc)



## Installation
We recommend using conda to manage the installation of all dependencies. To do this, simply run:

```
conda create --name SpotGF
conda activate SpotGF
# conda config --add channels conda-forge ##If you do not have this channel added before#
conda install python=3.8 pandas pot=0.8.2 numpy scipy matplotlib descartes
```
Then, download this repo and install it.
```
git clone [repo_path]
cd [path/to/SpotGF-main]
pip install .
```

The total installation time is around 10 mintunes. If error occuors, please upgrade pip and try again.


## Usage
Once the input data have been processed into the supported format, the full SpotGF workflow can be run by calling the `SpotGF.py` script. The input files can include various formats such as `gem`, `txt`, `csv`, `gem.gz`and others, containing the raw SRT data information. These files must contain specific columns including `geneID`, `x`, `y`, `MIDCount`. 

The denoising resolution can be adjusted using the `binsize` parameter. A smaller `binsize` will result in a finer denoising effect but will also increase the processing time.Please note that the `binsize` parameter is only used when calculating SpotGF scores and it does not affect the final denoised data.

The `proportion` parameter determines the proportion of genes to be retained in the final denoised data. For example, setting proportion=0.9 will retain only 90% of the effective genes, resulting in a new denoised SRT dataset. 

If the input file contains columns named `cen_x` and `cen_y`, the denoising process can be performed based on the cell bin, achieving denoising at the single-cell level by setting binsize=1. Please note that the above description provides an overview of the functionality and parameters. Let me know if there is anything else I can help you with.

```
python [path/to/SpotGF.py] [path/to/input.gem] [binsize] [proportion]
```

If you use SpotGF in Jupyter environment, you can choose blow Usage.

```
from SpotGF import SpotGF	

gem_path = './SpotGF/test/demo.gem' 
output = "./SpotGF/test"
binsize =75
lower= 0
upper = 100000
proportion  = 0.1
max_iterations = 10000
auto_threshold = True

#initial class
spotgf = SpotGF.SpotGF(gem_path,binsize,proportion,auto_threshold,lower,upper,max_iterations,output)

#acquire SpotGF socres distribution figure to help choose proportion 
GF_df = spotgf.calculate_GFscore(gem_path,binsize)

#denoised SRT data
new_gem  = spotgf.generate_GFgem(gem_path,GF_df,proportion,auto_threshold)
```


### Data preparation
Sprod workflow requires two mandatory files, a `input.gem` (with "\t" or "," as the delimiter) for gene expression data,

|geneID|x|y|MIDCount|
|-----|-----|-----|-----|
|#gene1|23|36|1|
|#gene1|24|35|1|
|#gene1|23|30|1|
|#gene2|20|31|1|
|#gene2|21|22|1|


### Output files
SpotGF_scores.txt: This file contains the SpotGF scores for each gene. The SpotGF score indicates the degree of clustering or diffusion of a gene. A smaller SpotGF score suggests that the gene is more diffuse, while a larger SpotGF score indicates that the gene is more clustered.


### List of Parameters
```
positional arguments:

  i            Input SRT data files path.

  o            Outpath for saving results.
    
  b            Denoising resolution binsize, must int type, default=5.

  lower            Lower limit for tissue structures capturing optimization, default=0.

  upper            Upper limit for tissue structures capturing optimization', default=sys.float_info.max.

  max_iterations            maximum number of iterations when capturing tissue structures', default=10000.

  p            Proportion of gene numbers, must float type [0,1], default=0.6.

  auto_threshold            if True, return a denoised gem file using automatic threshold, default=True.
```

### Contact Us
If you have any suggestions/ideas for SpotGF or are having issues trying to use it, please don't hesitate to reach out to us.

Lin Du, dulin[dot]@genomics[dot]cn 


### Cite
Du, L., Kang, J., Hou, Y., Sun, H. X., & Zhang, B. (2024). SpotGF: Denoising spatially resolved transcriptomics data using an optimal transport-based gene filtering algorithm. Cell systems, 15(10), 969–981.e6. https://doi.org/10.1016/j.cels.2024.09.005
