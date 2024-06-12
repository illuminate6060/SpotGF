# GeneFilter
Optimal Transport Method-Based Gene Filter (GF) Denoising Algorithm for Enhancing Spatially Resolved Transcriptomics Data

## Installation
We recommend using conda to manage the installation of all dependencies. To do this, simply run:

```
conda create --name GF
conda activate GF
# conda config --add channels conda-forge ##If you do not have this channel added before#
conda install python=3.7 pandas pot=0.8.2 numpy scipy matplotlib descartes
```
Then, download this repo and install it.
```
git clone [repo_path]
cd [path/to/GeneFilter-main]
pip install .
```

The total installation time is around 2 mintunes. If error occuors, please upgrade pip and try again.


## Usage
Once the input data have been processed into the supported format, the full sprod workflow can be run by calling the `GenenFilter.py` script. The input files can include various formats such as `gem`, `txt`, `csv`, and others, containing the raw SRT data information. These files must contain specific columns including `geneID`, `x`, `y`, `MIDCount`. The denoising resolution can be adjusted using the `binsize` parameter. A smaller `binsize` will result in a finer denoising effect but will also increase the processing time. The `proportion` parameter determines the proportion of genes to be retained in the final denoised data. For example, setting proportion=0.9 will retain only 90% of the effective genes, resulting in a new denoised SRT dataset.Set the auto_threshold parameter to True, GF will automatically generate a threshold based on the distribution of GF scores and save the denoised file. If the input file contains columns named `cen_x` and `cen_y`, the denoising process can be performed based on the cell bin, achieving denoising at the single-cell level by setting binsize=1. Please note that the above description provides an overview of the functionality and parameters. Let me know if there is anything else I can help you with.

**Method 1**
```
python [path/to/GenenFilter.py] [path/to/input.gem] [binsize] [proportion] [auto_threshold] [lower] [upper] [max_iterations] 
```

**Method 2**
If you use GeneFilter in Jupyter environment, you can choose blow Usage.

```
import GeneFilter as GF

gem_path = '../arbi-image1-cellbin.gem' 
binsize = 5
proportion = 0.6
auto_threshold = True

GF = GF.GeneFilter(gem_path,binsize,proportion,auto_threshold,lower,upper,max_iterations) #initial class
GF_df = GF.calculate_GFscore(gem_path,binsize)#acquire GF socres distribution figure to help choose proportion parameter
new_gem  = GF.generate_GFgem(gem_path,GF_df,proportion,auto_threshold ) #denoised SRT data
```


### Data preparation
Sprod workflow requires two mandatory files, a `input.gem` (with "\t" as the delimiter) for gene expression data,

|geneID|x|y|MIDCount|
|-----|-----|-----|-----|
|#gene1|23|36|1|
|#gene1|24|35|1|
|#gene1|23|30|1|
|#gene2|20|31|1|
|#...|...|...|...|
|#gene2|21|22|1|


### Output files
**GF_scores.txt**: This file contains the GF scores for each gene. The GF score indicates the degree of clustering or diffusion of a gene. A smaller GF score suggests that the gene is more diffuse, while a larger GF score indicates that the gene is more clustered.

**GF_auto_threshold.gem**: Find the reference threshold according to the automatic threshold method, and generate the gem file after denoising

**GF_proportion.gem**: The gem file generated after denoising according to the custom retained gene ratio


### List of Parameters
```
positional arguments:

  i            type=str, help='input gem file path'
  
  b            type=int, help='Denoising resolution binsize', default = 5

  lower            type=float, help='lower limit for tissue structures capturing optimization', default = 0

  lower            type=float, help='lower limit for tissue structures capturing optimization', default = 0

  upper            type=float, help='upper limit for tissue structures capturing optimization', default=sys.float_info.max

  max_iterations            type=int, help='maximum number of iterations when capturing tissue structures', default== 10000

  p            type=float, help='Proportion of matained genes, must float type [0,1]', default=0.5

  auto_threshold            type=bool, help='Whether generate GF-denoised data based on automatic threshold', default=True





### Contact Us
If you have any suggestions/ideas for GeneFilter or are having issues trying to use it, please don't hesitate to reach out to us.

Lin Du, dulin[dot]@genomics[dot]cn 
