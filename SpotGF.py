import numpy as np
import ot
import ot.plot
import pandas as pd
import matplotlib.pyplot as plt
import heapq
import math
import alphashape
from descartes import PolygonPatch
from shapely.geometry import Point
from multiprocessing import Pool
import argparse
import csv
import sys
import os
from scipy import sparse
import anndata
import scanpy as sc
import seaborn as sns

class SpotGF():
    def __init__(self,gem_path,binsize,proportion,auto_threshold,lower,upper,max_iterations,outpath,visualize,spot_size,alpha): 
        self.gem_path = gem_path
        self.binsize = binsize
        self.proportion = proportion
        self.auto_threshold = auto_threshold
        self.lower = lower
        self.upper = upper
        self.max_iterations = max_iterations
        self.outpath = outpath
        self.visualize = visualize
        self.spot_size = spot_size
        self.alpha = alpha
        os.chdir(outpath)
        print("Input file:",gem_path)
        print("Output path:",outpath)
        print("Denoising resolution :",binsize)

    def gem2adata(self,result):
        result['x_y'] = result['x'].astype('str')+'_'+result['y'].astype('str')
        result = result[['x_y', 'geneID', 'MIDCount']]
        cell_list = result["x_y"].astype('category')
        gene_list = result["geneID"].astype('category')
        data =  result["MIDCount"].to_numpy()
        row = cell_list.cat.codes.to_numpy()
        col = gene_list.cat.codes.to_numpy()
        obs = pd.DataFrame(index = cell_list.cat.categories)
        var = pd.DataFrame(index = gene_list.cat.categories)
        X = sparse.csr_matrix((data, (row, col)), shape = (len(obs), len(var)))
        adata = anndata.AnnData(X, obs = obs, var = var)
        adata.obsm['spatial'] = pd.Series(adata.obs.index).str.split('_', expand=True).astype(int).values
        adata.obs['spatial1'] = adata.obsm['spatial'][:,1].astype('int64')
        adata.obs['spatial2'] = adata.obsm['spatial'][:,0].astype('int64')
        return adata

    def convert_x_y_to_numeric(self,df):
        """
        This function converts the 'x' and 'y' columns of the DataFrame to numeric types 
        and then casts them to 64-bit integers.
        """
        df['x'] = pd.to_numeric(df['x'])
        df['x'] = df['x'].astype(np.int64)
        df['y'] = pd.to_numeric(df['y'])
        df['y'] = df['y'].astype(np.int64)
        return df

    def open_gem(self,gem_path):
        """
        This function opens a GEM file, reads its content, determines the delimiter using csv Sniffer,
        reads the csv data (handling gzip compression if necessary), renames specific columns if they 
        exist, converts 'x' and 'y' columns to numeric types, and returns the processed DataFrame.
        """
        with open(gem_path, 'r') as file:
            sample_data = file.read(1024)
        dialect = csv.Sniffer().sniff(sample_data)
        if gem_path.endswith('.gz'):
            ex_raw = pd.read_csv(gem_path, delimiter=dialect.delimiter, compression= 'gzip', comment='#')
        else:
            ex_raw = pd.read_csv(gem_path, delimiter=dialect.delimiter, comment='#')
        if "UMICount" in ex_raw.columns:
            ex_raw = ex_raw.rename(columns={'UMICount':'MIDCount'})
        if "UMICounts" in ex_raw.columns:
            ex_raw = ex_raw.rename(columns={'UMICounts':'MIDCount'})
        if "MIDCounts" in ex_raw.columns:
            ex_raw = ex_raw.rename(columns={'MIDCounts':'MIDCount'})    
        ex_raw = self.convert_x_y_to_numeric(ex_raw)
        return ex_raw
    
    def preparedata(self,ex_raw):
        """
        This function processes raw expression data to filter out genes with fewer than 10 occurrences,
        acquires unique gene names and cell positions, sums up 'MIDCount' for each 'geneID'-'x'-'y' group,
        and returns the processed data along with unique cell and gene lists.
        """
        count = ex_raw['geneID'].value_counts()
        gsave = count[count > 10].index.tolist()
        ex_raw = ex_raw[ex_raw.geneID.isin(gsave)].reset_index(drop=True)
        all_gene = ex_raw['geneID'].unique()
        print("Valid genes Number:",len(all_gene))
        ex2 = ex_raw.groupby(['geneID','x','y'], as_index=False)['MIDCount'].sum()
        all_cell = ex_raw[['x', 'y']].drop_duplicates().values
        print("Valid cells number:",len(all_cell)) 
        return ex2,all_cell,all_gene

    def grid_downsample(self,points, num_points, mask_area, alpha_shape):
        """
        This function downsamples a set of points using a grid-based approach. 
        It divides the bounding box of the input points into a grid, counts points in each cell, 
        and selects representative points that lie within a specified alpha shape.
        """
        # Calculate the bounding box of the input points
        min_x, min_y = np.min(points, axis=0)
        max_x, max_y = np.max(points, axis=0)
        width = max_x - min_x
        height = max_y - min_y

        # Determine the bin size based on the mask area and desired number of points
        bin_size = np.sqrt(mask_area / num_points)
        num_cols = math.ceil(width / bin_size)
        num_rows = math.ceil(height / bin_size)

        # Create an output grid to count points in each cell
        grid = np.zeros((num_rows, num_cols), dtype=int)
        
        # Assign input points to grid cells
        cols = ((points[:, 0] - min_x) / bin_size).astype(int)
        rows = ((points[:, 1] - min_y) / bin_size).astype(int)
        for col, row in zip(cols, rows):
            if 0 <= col < num_cols and 0 <= row < num_rows:
                grid[row, col] += 1
        
        # Select points from the grid
        output_points = []
        for row in range(num_rows):
            for col in range(num_cols):
                if grid[row, col] > 0:
                    x = (col + 0.1) * bin_size + min_x
                    y = (row + 0.1) * bin_size + min_y
                    if alpha_shape.contains(Point(x, y)):  # Filter points within the alpha shape
                        output_points.append([x, y])
        
        output_points = np.array(output_points)
        return output_points

    def calculate_ot(self,all_gene,ex2,all_cell,alpha_shape):
        """
        Calculate the optimal transport (OT) for given genes and cells.

        Parameters:
        - all_gene: List of all gene IDs.
        - ex2: DataFrame containing gene expression data with columns 'geneID', 'x', 'y', 'MIDCount'.
        - all_cell: List of cell coordinates.
        - alpha_shape: alpha_shape for tissue counter.

        Returns:
        - result: Array with gene IDs and their corresponding OT values.
        """
        i=0
        emd = []
        gene = []
        mask_area = alpha_shape.area
        grouped_ex2 = ex2.groupby('geneID') 
        for n in  all_gene: 
            gene_c = grouped_ex2.get_group(n) 
            gene.append(n)
            i = i+1
        
            #1.create source distribution
            wnc = gene_c.groupby(['x', 'y']).size().reset_index().rename(columns={0: 'count'})  # Group by coordinates
            source_point = wnc[['x', 'y']].to_numpy()  # Convert to numpy array
            source_w = wnc['count'].to_numpy(dtype='float64')  # Get weights as float

            if len(source_point) > 5000:  #downsample to 5000 points if more than 5000
                source_point = self.grid_downsample(source_point, 5000, mask_area, alpha_shape) 
                source_w = np.ones(len(source_point))  
            
            #2.create target diatribution
            # print(( len(source_point)) )
            target_point = self.grid_downsample(all_cell, len(source_point), mask_area, alpha_shape) #keep same points
            target_w = np.ones(len(target_point), dtype='float64') * (sum(source_w) / len(target_point))  # Normalize weights
            target_w =  target_w.astype('float64')

            #3.calc OT
            M = ot.dist(source_point, target_point, metric='euclidean')
            result2 = ot.emd2(source_w,target_w, M)  / len(source_point)  #nprmalize
            emd.append(result2)
            # print(i,n,result2)
            
        result = np.array([gene,emd])
        return result
        
    def calculate_GFscore(self,gem_path,binsize,alpha,max_iterations,lower,upper):
        """
        Calculate SpotGF scores for genes based on spatial expression data.

        Parameters:
        - gem_path: Path to the GEM file.
        - binsize: Size of the bins for spatial binning.
        - alpha: alpha parameter for detect tissue counter

        Returns:
        - GF_df: DataFrame containing genes and their corresponding SpotGF scores.
        """
        ex_raw = self.open_gem(gem_path) 
        if binsize == 1:
            print('Process with resolution cell bin data')
            ex_raw.rename(columns={'cen_y':'y'},inplace=True) #For cell-bin data
            ex_raw.rename(columns={'cen_x':'x'},inplace=True)
        else:
            print('Process with resolution binsize',str(binsize) ,'data')
            ex_raw['x'] = ex_raw['x'].map(lambda x: int(x/binsize))
            ex_raw['y'] =ex_raw['y'].map(lambda y: int(y/binsize))
        
        ex2,all_cell,all_gene = self.preparedata(ex_raw)

        # decide alpha value
        if alpha == 0:
            alpha_use = alphashape.optimizealpha(all_cell,max_iterations,lower,upper)
        else:
            alpha_use = alpha
        alpha_shape = alphashape.alphashape(all_cell, alpha_use)
        # visualize 
        plt.figure()
        if alpha_shape.geom_type == 'Polygon':
            x, y = alpha_shape.exterior.xy
            plt.fill(x, y, alpha=0.5, fc='#3fc1c9', ec='#393e46')
        elif alpha_shape.geom_type == 'MultiPolygon':
            for polygon in alpha_shape:
                x, y = polygon.exterior.xy
                plt.fill(x, y, alpha=0.5, fc='#3fc1c9', ec='#393e46')
        plt.scatter(all_cell[:, 0], all_cell[:, 1], color='#f38181')
        plt.title('Alpha Shape')
        plt.xlabel('X-axis')
        plt.ylabel('Y-axis')
        plt.savefig('alpha_shape.png', dpi=300)
        plt.close()

        args_list = [(all_gene,ex2,all_cell,alpha_shape)]

        # Data size check
        if len(all_cell) > 50000:  # You can adjust this threshold as needed
            print("Data is too large, we recommend that you appropriately increase the resolution binzise")

        # 使用 multiprocessing 并行计算 OT 分数
        with Pool(processes=20) as pool:
            result = pool.starmap(self.calculate_ot, args_list)

        result = np.array(result).reshape(2, len(all_gene))
        GF_df = pd.DataFrame(result.T,columns=['geneID','SpotGF_score'])
        GF_df.to_csv('SpotGF_scores.txt',sep='\t')
        print("Finished saving SpotGF_scores.txt")

        return GF_df
        
    def cal_threshold(self,emd2):
        """
        Calculate the threshold point based on the change in Earth Mover's Distance (EMD).

        Parameters:
        - emd2: List or array of EMD values.

        Returns:
        - max_point: The point with the maximum change in slope, used as a threshold.
        """

        # 1. Calculate the first derivative to represent the change in EMD
        x = np.arange(len(emd2))  # Construct x-axis data
        emd2 = np.array(emd2, dtype=float)
        dydx = np.gradient(emd2, x)  # Calculate the first derivative

        # 2. Adjust dydx to limit the range of the first derivative
        dydx = np.array(dydx)
        mean_dydx = np.mean(dydx)
        dydx[dydx > mean_dydx] = mean_dydx  # Limit dydx values
        dydx[dydx > mean_dydx] = mean_dydx  # Limit dydx values again for reinforcement

        # 3. Smooth the curve using polynomial fitting
        cur = np.polyfit(x, dydx, 10)  # Polynomial fitting, returns coefficients
        p1 = np.poly1d(cur)  # Generate polynomial expression

        # 4. Find the point with the maximum change in slope (second derivative)
        dy = np.gradient(p1(x), x)  # Calculate the gradient (slope) of the polynomial
        max_idx = np.argmax(dy)  # Find the index of the maximum slope
        max_point = (x[max_idx], emd2[max_idx])  # Get the corresponding (x, y) point

        # Ensure the threshold filters at least 50% of the genes by adjusting the polynomial degree
        for i in [15, 20, 25, 30]:
            if x[max_idx] <= np.median(x):  # Check if less than 50% of the genes are filtered
                cur = np.polyfit(x, dydx, i)  # Polynomial fitting with a higher degree
                p1 = np.poly1d(cur) 
                dy = np.gradient(p1(x), x)  # Recalculate the gradient
                max_idx = np.argmax(dy)  # Find the new maximum slope point
                max_point = (x[max_idx], emd2[max_idx])  # Update the maximum point

        print('thred:', max_point)
        return max_point  # Return the threshold point

    def expression_figure(adata,save_path,spot_size):
        if len(adata.var) > 200:
            adata.var["mt"] = adata.var_names.str.startswith("MT")
            sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"],percent_top=(50,100,200),inplace=True,log1p=True)
            sc.pl.spatial(adata, color = 'total_counts', spot_size=spot_size,title='Total_counts', show=False, return_fig =True,color_map ='Spectral_r') 
            ax = plt.gca()
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            # ax.spines['left'].set_visible(False)
            # ax.spines['bottom'].set_visible(False)
            ax.tick_params(axis='x', labelsize=20)
            ax.tick_params(axis='y', labelsize=20)
            ax.tick_params(width=3) 
            ax.spines['left'].set_linewidth(3) 
            ax.spines['bottom'].set_linewidth(3) 
            plt.savefig(save_path, dpi=300)
            plt.close()
            return adata 
        else:
            print("Warning: Gene numbers below 200 cannot visualize denoised data")

        
    def generate_GFgem(self,gem_path,GF_df,proportion,auto_threshold,visualize,spot_size):
        """
        Generate a gem file filtered by SpotGF scores.

        Parameters:
        - ex_raw: Raw expression data.
        - SpotGF_scores: DataFrame containing genes and their SpotGF scores.
        - proportion: Proportion of top genes to save based on SpotGF scores.
        - auto_threshold: Boolean flag to determine whether to use automatic thresholding.
        - visualize: whether visualize for denoised data.
        - spot_size: spot_size parameterfor draw spatial figure

        Returns:
        - result: Filtered expression data based on SpotGF scores.
        """
        ex_raw = self.open_gem(gem_path) 
        count = ex_raw['geneID'].value_counts()
        gsave_high = count[count <= 10].index.tolist()    

        GF_df['SpotGF_score']  = GF_df['SpotGF_score'].astype(float)
        GF_df = GF_df.sort_values(by='SpotGF_score')
        emd2 = GF_df['SpotGF_score'].tolist()

        if auto_threshold == True:
            print("Generate and visualize SpotGF-denoised data based on automatic threshold")
            thred = self.cal_threshold(emd2)
            save_gene = GF_df[GF_df['SpotGF_score'] >= float(thred[1]) ]['geneID'].tolist()
            save = save_gene + gsave_high
            result = ex_raw[ex_raw.geneID.isin(save)]
            result = result.drop(columns=['Unnamed: 0.1', 'Unnamed: 0'], errors='ignore')
            result.to_csv('SpotGF_auto_threshold.gem','\t')
            if visualize == True:
                # print("Visualize SpotGF-denoised data based on automatic threshold")
                #os.makedirs('automatic', exist_ok=True)
                save_path = './Spatial_automatic.png'
                adata = self.gem2adata(result)
                adata_auto = self.expression_figure(adata,save_path,spot_size)
            
        if proportion is not None:
            print("Generate and visualize SpotGF-denoised data based on proportion")
            drop_pre = int(len(GF_df) / 10  * proportion )
            save_pro = GF_df[GF_df['SpotGF_score'] >= heapq.nlargest(drop_pre,GF_df['SpotGF_score'])[-1]]
            save_gene = list(save_pro.geneID) 
            save = save_gene + gsave_high
            result = ex_raw[ex_raw.geneID.isin(save)]
            result = result.drop(columns=['Unnamed: 0.1', 'Unnamed: 0'], errors='ignore')
            result.to_csv('SpotGF_proportion_'+str(proportion)+'.gem', sep='\t') 
            if visualize == True:
                # print("Visualize SpotGF-denoised data based on proportion")
                #os.makedirs('proportion', exist_ok=True)
                save_path = './Spatial_proportion.png'
                adata = self.gem2adata(result)
                adata_prop = self.expression_figure(adata,save_path,spot_size)
        
        adata_raw = self.gem2adata(ex_raw)
        adata_raw.var["mt"] = adata_raw.var_names.str.startswith("MT")
        sc.pp.calculate_qc_metrics(adata_raw, qc_vars=["mt"],percent_top=(50,100,200),inplace=True,log1p=True)
        df_raw = pd.DataFrame(adata_raw.obs)[['n_genes_by_counts', 'total_counts']]
        df_raw['type'] = ['Raw'] * len(df_raw)

        if auto_threshold == True:
            df_auto = pd.DataFrame(adata_auto.obs)[['n_genes_by_counts', 'total_counts']]
            df_auto['type'] = ['Auto_denoised'] * len(df_auto)
            if proportion is not None:
                df_prop= pd.DataFrame(adata_prop.obs)[['n_genes_by_counts', 'total_counts']]
                df_prop['type'] = ['Prop_denoised'] * len(df_prop)
                df_all = pd.concat([df_raw,df_auto,df_prop])
            else:
                df_all = pd.concat([df_raw,df_auto])
            save_path = "Violinplot_total_counts.png"
            plt.figure(dpi=300, figsize=(5,8))
            sns.violinplot(data=df_all, x='type',y='total_counts',linewidth=1,inner='box', palette = "Set3", hue='type', legend=False)   
            plt.title('Violin Plot of Total Counts by Type', fontsize=14)    
            plt.xticks(rotation=90)
            plt.tight_layout()
            plt.savefig(save_path, dpi=300)
            plt.close()
            save_path = "Violinplot_n_genes_by_counts.png"
            plt.figure(dpi=300, figsize=(5,8))
            g = sns.violinplot(data=df_all, x='type',y='n_genes_by_counts',linewidth=1,inner='box', palette="Set3", hue='type', legend=False)   
            plt.title('Violin Plot of n_genes_by_counts', fontsize=14)    
            plt.xticks(rotation=90)
            plt.tight_layout()  
            plt.savefig(save_path, dpi=300)
            plt.close()  
        return result

if __name__ == '__main__':
    '''
    gem_path (str): Input SRT data files.such as input.gem
    outpath(str): Output path for saving results
    binsize (int): Denoising resolution binsize
    proportion (float): Proportion of matained genes, must float type [0,1]
    auto_threshold(bool): Generate SpotGF-denoised data based on automatic threshold,default = True
    lower (float): lower limit for tissue structures capturing optimization , default = 0
    upper (float): upper limit for tissue structures capturing optimization , default = sys.float_info.max
    max_iterations (bool): maximum number of iterations when capturing finding tissue structures
    visualize (bool): Visualize SpotGF-denoised data
    spot siz (int): Spot size for visualize SpotGF-denoised data
    alpha (float): Alpha for counter detection,default use auto optimizealpha
    '''
    parser = argparse.ArgumentParser(description='SpotGF')
    parser.add_argument('-i', '--arg1', type=str, help='input gem file path')
    parser.add_argument('-o', '--arg2', type=str, help='outpath for saving results')
    parser.add_argument('-b', '--arg3', type=int, help='Denoising resolution binsize', default=5)
    parser.add_argument('-lower', '--arg4', type=float, help='lower limit for tissue structures capturing optimization', default=0)
    parser.add_argument('-upper', '--arg5', type=float, help='upper limit for tissue structures capturing optimization', default=sys.float_info.max)
    parser.add_argument('-max_iterations', '--arg6', type=int, help='maximum number of iterations when capturing tissue structures', default=10000)
    parser.add_argument('-p', '--arg7', type=float, help='Proportion of matained genes, must float type [0,1]', default=0.5)
    parser.add_argument('-auto_threshold', '--arg8', type=bool, help='Automatic threshold', default=True)
    parser.add_argument('-v', '--arg9', type=bool, help='Visualize SpotGF-denoised data', default=True)
    parser.add_argument('-s', '--arg10', type=int, help='Spot size for visualize SpotGF-denoised data', default=5)
    parser.add_argument('-a', '--arg11', type=float, help='Alpha for counter detection,default use auto optimizealpha', default=0 )
    args = parser.parse_args()

    gem_path = args.arg1
    outpath = args.arg2
    binsize = args.arg3
    lower = args.arg4
    upper = args.arg5
    max_iterations = args.arg6
    proportion = args.arg7
    auto_threshold = args.arg8
    visualize = args.arg9
    spot_size = args.arg10
    alpha = args.arg11
    
    os.makedirs(outpath, exist_ok=True)
    spotgf = SpotGF(gem_path,binsize,proportion,auto_threshold,lower,upper,max_iterations,outpath,visualize,spot_size,alpha)
    GF_df = spotgf.calculate_GFscore(gem_path,binsize,alpha,max_iterations,lower,upper)
    new_gem  = spotgf.generate_GFgem(gem_path,GF_df,proportion,auto_threshold,visualize,spot_size)
