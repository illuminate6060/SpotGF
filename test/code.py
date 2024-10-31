import SpotGF

gem_path = './demo.gem'  #please change this file path base on your work directory
output =  './result'  #please change this output path base on your work directory


binsize =70
lower= 0
upper = 100000
proportion  = 0.1
max_iterations = 10000
auto_threshold = True
visualize=True
spot_size=5
alpha=0


spotgf =SpotGF.SpotGF(gem_path,binsize,proportion,auto_threshold,lower,upper,max_iterations,output,visualize,spot_size,alpha)
GF_df = spotgf.calculate_GFscore(gem_path,binsize,alpha)
new_gem  = spotgf.generate_GFgem(gem_path,GF_df,proportion,auto_threshold,visualize,spot_size)


"""
other usage
python ../SpotGF.py -i ./demo.gem -o ./result -b 70 -s 5
"""