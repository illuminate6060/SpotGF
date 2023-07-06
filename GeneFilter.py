import numpy as np
import ot
import ot.plot
import pandas as pd
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull, distance
import heapq
import math
import alphashape
from descartes import PolygonPatch
from shapely.geometry import Point
from multiprocessing import Pool
plt.rcdefaults()
import argparse
import gzip
import csv


class GeneFilter():
    def __init__(self,gem_path,binsize,proportion): #outpath,
        self.gem_path = gem_path
        self.binsize = binsize
        self.proportion = proportion

    def convert_x_y_to_numeric(self,df):
        df['x'] = pd.to_numeric(df['x'])
        df['x'] = df['x'].astype(np.int64)
        df['y'] = pd.to_numeric(df['y'])
        df['y'] = df['y'].astype(np.int64)
        return df

    def point_in_hull(point, hull, tolerance=1e-12):
        return all(
            (np.dot(eq[:-1], point) + eq[-1] <= tolerance)
            for eq in hull.equations)
    
    def open_gem(self,gem_path):
        with open(gem_path, 'r') as file:
            sample_data = file.read(1024)
        dialect = csv.Sniffer().sniff(sample_data)
        if gem_path.endswith('.gz'):
            ex_raw = pd.read_csv(gem_path, delimiter=dialect.delimiter, compression= 'gzip' )
        else:
            ex_raw = pd.read_csv(gem_path, delimiter=dialect.delimiter)
        
        if "UMICount" in ex_raw.columns:
            ex_raw = ex_raw.rename(columns={'UMICount':'MIDCount'})
        if "MIDCounts" in ex_raw.columns:
            ex_raw = ex_raw.rename(columns={'MIDCounts':'MIDCount'})    
        ex_raw = self.convert_x_y_to_numeric(ex_raw)
        return ex_raw
    
    def preparedata(self,ex_raw):
        #Filter genes less than 10
        count = ex_raw['geneID'].value_counts()
        gsave = count[count > 10].index.tolist()
        ex_raw = ex_raw[ex_raw.geneID.isin(gsave)].reset_index() 
        #Aquire all gene name
        all_gene = np.unique(ex_raw.geneID)
        print("Genes Number:",len(all_gene))
        ex2 = ex_raw.groupby(['geneID','x','y'])['MIDCount'].sum().reset_index() 
        #Aquire all cells
        all_c = list(ex_raw.groupby(['x','y']).groups.keys())
        all_cell = np.asarray(all_c)
        print("Spots number:",len(all_cell)) 
        return ex2,all_cell,all_gene

    def grid_downsample(self,points, num_points):
        hull = ConvexHull(points)
        mask_area  = hull.volume
        
        min_x, min_y = np.min(points, axis=0)
        max_x, max_y = np.max(points, axis=0)
        width = max_x - min_x
        height = max_y - min_y
        # print(min_x, min_y ,max_x, max_y)
        # print(width,height)

        # bin_size = math.sqrt((width * height) / num_points)
        bin_size = math.sqrt(mask_area / num_points)   
        # print("bin_size",bin_size)
        
        cols = ((points[:, 0] - min_x) / bin_size).astype(int)
        rows = ((points[:, 1] - min_y) / bin_size).astype(int)
        grid = np.zeros((rows.max()+1 , cols.max()+1), dtype=int) 
        # print('num_cols:',rows.max()+1,'num_rows:',cols.max()+1)
        for col, row in zip(cols, rows):
            grid[row, col] += 1

        output_points = []
        tolerance = 1e-12
        for row in range(rows.max()+1):
            for col in range(cols.max()+1):
                if grid[row, col] > 0:   
                    # print(grid[row, col])
                    x = (col + 0.1) * bin_size + min_x
                    y = (row + 0.1) * bin_size + min_y
                    all((np.dot(eq[:-1], (x,y)) + eq[-1] <= tolerance)
                        for eq in hull.equations)
                    # if  self.point_in_hull((x,y),hull):
                    # if alpha_shape.contains(Point(x,y)): 
                    output_points.append([x, y])
        output_points = np.array(output_points)
        return output_points

    def calculate_ot(self,all_gene,ex2,all_cell):
        i=0
        emd = []
        gene = []
        for n in  all_gene: 
            gene_c= ex2.loc[ex2['geneID'] == n] 
            gene.append(n)
            i = i+1
            print(i,n,len(gene_c))
            #1.create source distribution
            wnc = list(gene_c.groupby(['x', 'y']).groups.keys())
            source_point = np.asarray(wnc)
            source_w = np.asarray(gene_c['MIDCount'])   #/ gene_c['MIDCount'].sum()    #gene
            source_w = source_w.astype('float64')
            # print(source_w,len(source_w))
            if len(source_point) > 5000:  #downsample to 5000 points if more than 5000
                source_point = self.grid_downsample(source_point, 5000) 
                source_w = np.ones(len(source_point))  
            
            #2.create target diatribution
            taget_point = self.grid_downsample(all_cell, len(source_point)) #keep same points
            target_w = np.ones(len(taget_point))  * (sum(source_w) / len(taget_point))
            target_w =  target_w.astype('float64')
            # print(target_w,len(target_w))
            # print(len(taget_point),len(target_w), len(source_point),len(source_w))

            #3.calc OT
            M = ot.dist(source_point, taget_point, metric='euclidean')
            result2 = ot.emd2(source_w,target_w, M)  / len(source_point)  #nprmalize
            emd.append(result2)
        result = np.array([gene,emd])
        return result
    
    def calculate_GFscore(self,gem_path,binsize):
        ex_raw = self.open_gem(gem_path) 
        #include scale
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
        with Pool(processes=20) as pool:
            result = pool.starmap(self.calculate_ot, args_list)

        result = np.array(result).reshape(2, len(all_gene))
        df = pd.DataFrame(result.T,columns=['GeneID','GF_score'])
        df.to_csv('./GF_scores.txt',sep='\t')

        #save gf scores figures
        plt.rcdefaults()
        # plt.figure(dpi=300)
        plt.plot(np.arange(len(df['GF_score'])),df['GF_score'])
        plt.title('GF scores distribution')
        plt.xlabel('different gene')
        plt.ylabel('GF scores')
        plt.savefig('./GF_scores.pdf',dpi=100)

        return df 
    
    def generate_GFgem(self,gem_path,df,proportion):
        ex_raw = self.open_gem(gem_path) 
        drop_pre = int(len(df) / 10  * proportion )
        test = df[df['GF_score'] >= heapq.nlargest(drop_pre,df['GF_score'])[-1]]
        save = list(test['GeneID']) 
        ex_new = ex_raw[ex_raw.geneID.isin(save)]
        
        return ex_new


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='calculate ot')
    parser.add_argument('-i', '--arg1', type=str, help='input gem file path') 
    parser.add_argument('-b', '--arg2', type=int, help='binsize')
    parser.add_argument('-p', '--arg3', type=float, help='proportion save genes')
    args = parser.parse_args()

    gem_path = args.arg1
    binsize = args.arg2
    proportion = args.arg3

    GF = GeneFilter(gem_path,binsize,proportion)
    result = GF.calculate_GFscore(gem_path,binsize)
    new_gem  = GF.generate_GFgem(gem_path,result,proportion)
    new_gem.to_csv('./new_gem.gem',sep='\t')
