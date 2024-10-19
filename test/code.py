import SpotGF

gem_path = 'path/SpotGF/test/demo.gem'  #please change this file path base on your work directory
output =  'path/SpotGF_test/SpotGF/test'  #please change this output path base on your work directory

binsize =75
lower= 0
upper = 100000
proportion  = 0.1
max_iterations = 10000
auto_threshold = True

spotgf = SpotGF.SpotGF(gem_path,binsize,proportion,auto_threshold,lower,upper,max_iterations,output)
GF_df = spotgf.calculate_GFscore(gem_path,binsize)
new_gem  = spotgf.generate_GFgem(gem_path,GF_df,proportion,auto_threshold)