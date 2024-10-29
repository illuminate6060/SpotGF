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

class SpotGF():
    def __init__(self,gem_path,binsize,proportion,auto_threshold,lower,upper,max_iterations,outpath): 
        self.gem_path = gem_path
        self.binsize = binsize
        self.proportion = proportion
        self.auto_threshold = auto_threshold
        self.lower = lower
        self.upper = upper
        self.max_iterations = max_iterations
        self.outpath = outpath
        os.chdir(outpath)
        print("input file:",gem_path)
        print("output path:",outpath)

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
        print("Genes Number:",len(all_gene))
        ex2 = ex_raw.groupby(['geneID','x','y'], as_index=False)['MIDCount'].sum()
        all_cell = ex_raw[['x', 'y']].drop_duplicates().values
        print("Cells number:",len(all_cell)) 
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

    def calculate_ot(self,all_gene,ex2,all_cell):
        """
        Calculate the optimal transport (OT) for given genes and cells.

        Parameters:
        - all_gene: List of all gene IDs.
        - ex2: DataFrame containing gene expression data with columns 'geneID', 'x', 'y', 'MIDCount'.
        - all_cell: List of cell coordinates.

        Returns:
        - result: Array with gene IDs and their corresponding OT values.
        """
        i=0
        emd = []
        gene = []
        alpha = alphashape.optimizealpha(all_cell)
        alpha_shape = alphashape.alphashape(all_cell, alpha)
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
        
    def calculate_GFscore(self,gem_path,binsize):
        """
        Calculate SpotGF scores for genes based on spatial expression data.
        
        Parameters:
        - gem_path: Path to the GEM file.
        - binsize: Size of the bins for spatial binning.
        
        Returns:
        - df: DataFrame containing genes and their corresponding SpotGF scores.
        """
        print("Calculate SpotGF scores for genes with binsize ",binsize)

        ex_raw = self.open_gem(gem_path) 
        if binsize == 1:
            print('Process with Cell Bin data')
            ex_raw.rename(columns={'cen_y':'y'},inplace=True) #For cell-bin data
            ex_raw.rename(columns={'cen_x':'x'},inplace=True)
        else:
            print('Process with Bin',str(binsize) ,'data')
            ex_raw['x'] = ex_raw['x'].map(lambda x: int(x/binsize))
            ex_raw['y'] =ex_raw['y'].map(lambda y: int(y/binsize))
        
        ex2,all_cell,all_gene = self.preparedata(ex_raw)
        args_list = [(all_gene,ex2,all_cell)]

        # Data size check
        if len(all_cell) > 100000:  # You can adjust this threshold as needed
            print("Data is too large, we recommend that you appropriately increase the resolution")

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

    def generate_GFgem(self,gem_path,GF_df,proportion,auto_threshold):
        """
        Generate a GEM file filtered by SpotGF scores.
        
        Parameters:
        - ex_raw: Raw expression data.
        - SpotGF_scores: DataFrame containing genes and their SpotGF scores.
        - proportion: Proportion of top genes to save based on SpotGF scores.
        - auto_threshold: Boolean flag to determine whether to use automatic thresholding.
        
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
            print("Generate SpotGF-denoised data based on automatic threshold")
            thred = self.cal_threshold(emd2)
            save_gene = GF_df[GF_df['SpotGF_score'] >= float(thred[1]) ]['geneID'].tolist()
            save = save_gene + gsave_high
            result = ex_raw[ex_raw.geneID.isin(save)]
            result = result.drop(columns=['Unnamed: 0.1', 'Unnamed: 0'], errors='ignore')
            result.to_csv('SpotGF_auto_threshold.gem','\t')
        if proportion is not None:
            print("Generate SpotGF-denoised data based on proportion")
            drop_pre = int(len(GF_df) / 10  * proportion )
            save_pro = GF_df[GF_df['SpotGF_score'] >= heapq.nlargest(drop_pre,GF_df['SpotGF_score'])[-1]]
            save_gene = list(save_pro.geneID) 
            save = save_gene + gsave_high
            result = ex_raw[ex_raw.geneID.isin(save)]
            result = result.drop(columns=['Unnamed: 0.1', 'Unnamed: 0'], errors='ignore')
            result.to_csv('SpotGF_proportion_'+str(proportion)+'.gem', sep='\t') 
            return result



if __name__ == '__main__':
    '''
    gem_path (str): Input SRT data files.such as input.gem
    binsize (int): Denoising resolution binsize
    proportion (float): Proportion of matained genes, must float type [0,1]
    auto_threshold(bool): Generate SpotGF-denoised data based on automatic threshold,default = True
    lower (float): lower limit for tissue structures capturing optimization , default = 0
    upper (float): upper limit for tissue structures capturing optimization , default = sys.float_info.max
    max_iterations (int): maximum number of iterations when capturing finding tissue structures
    outpath(str): Output path for saving results
    '''
    parser = argparse.ArgumentParser(description='SpotGF')
    parser.add_argument('-i', '--arg1', type=str, help='input gem file path')
    parser.add_argument('-b', '--arg2', type=int, help='Denoising resolution binsize', default=5)
    parser.add_argument('-lower', '--arg3', type=float, help='lower limit for tissue structures capturing optimization', default=0)
    parser.add_argument('-upper', '--arg4', type=float, help='upper limit for tissue structures capturing optimization', default=sys.float_info.max)
    parser.add_argument('-max_iterations', '--arg5', type=int, help='maximum number of iterations when capturing tissue structures', default=10000)
    parser.add_argument('-p', '--arg6', type=float, help='Proportion of matained genes, must float type [0,1]', default=0.5)
    parser.add_argument('-auto_threshold', '--arg7', type=bool, help='automatic threshold', default=True)
    parser.add_argument('-o', '--arg8', type=str, help='outpath for saving results')
    args = parser.parse_args()

    gem_path = args.arg1
    binsize = args.arg2
    lower = args.arg3
    upper = args.arg4
    max_iterations = args.arg5
    proportion = args.arg6
    auto_threshold = args.arg7
    outpath = args.arg8

    spotgf = SpotGF(gem_path,binsize,proportion,auto_threshold,lower,upper,max_iterations,outpath)
    GF_df = spotgf.calculate_GFscore(gem_path,binsize)
    new_gem  = spotgf.generate_GFgem(gem_path,GF_df,proportion,auto_threshold )
