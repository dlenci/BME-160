#!/usr/bin/env python3
# Name: David Lenci (dlenci)
# Group Members: none

import pandas as pd
import math
import numpy as np
import feather

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
get_ipython().run_line_magic('matplotlib', 'inline')

from sklearn import linear_model
from sklearn.decomposition import PCA
from sklearn.cluster import AgglomerativeClustering
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import normalize
from sklearn.model_selection import cross_val_score, cross_val_predict
from sklearn import metrics


class findCells():
    '''
    Takes three inputs: single cell data, relevant gene list, and variable
    to model off of. Gene list an be text file, but single cell data and
    variable must be feather files. Gene data must be in the format where
    column names are genes and row names are cells. 
    Parses data by removing cells and genes that got poor coverage, then 
    removing any genes that are not present in provided list. This data
    is then normalized and clustered, and lastly is modeled with input
    variable. Output is R squared score and graphs of clustering and modeled
    data. 
    '''
    
    def __init__(self, dataFrame, genes, stem):
        '''
        Instantiates class and calls cleanUp function. Requires three parameters:
        gene data frame, relevant genes, and variable to model with.
        '''
        cellData = dataFrame
        coolGenes = genes
        self.stemData = stem
        
        self.cleanUp(cellData, coolGenes)
    
    def cleanUp(self, dataFrame, genes): 
        '''
        Removes irrelevant genes (those not in gene list), and removes any genes
        or cells that got poor coverage in sequencing. Requires gene data and 
        relevant gene list as inputs and outputs data frame of clean data.
        '''
        cleanCell = dataFrame
        cleanGenes = genes
        cleanGenes.add('barcode')
        
        for column in cleanCell: #Remove genes with poor coverage.
            sum1 = 0
            if cleanCell[column].dtype == 'int64':
                sum1 = cleanCell[column].sum()
                if sum1 == 0 or sum1 == 0.0:
                    cleanCell = cleanCell.drop(column, 1)
                                
        ind = 0
        for row in cleanCell.iterrows(): #Remove cells with poor coverage - becomes very slow with floating point values
            sum1 = sum(row)
            for j, column in row.iteritems():
                if type(column) == int: or type(column) == float:
                    sum1 += Column
            
            if sum1 == 0:
                cleanCell = cleanCell.drop(cleanCell.index[ind])
                ind = ind - 1
            
            ind += 1 
       
        myGenes = set(cleanCell) 
        toDrop = myGenes.difference(cleanGenes)
        cleanCell.drop(toDrop, axis=1, inplace=True) #Remove genes not in input gene list
        
        
        print('check3')    
        self.normalizeMe(cleanCell)
    
    def normalizeMe(self, dataFrame):
        '''
        Normalizes clean gene data. Input: data frame, output normalized data
        frame. Calls cluster function.
        '''
        normCell = dataFrame
        
        copyCell = normCell.copy()
        for name in normCell.columns:
            if normCell[name].dtype == 'int64' or if normCell[name].dtype == 'float64': #Normalize only columns with numbers
                maxValue = normCell[name].max()
                minValue = normCell[name].min()
                copyCell[name] = (normCell[name] - minValue) / (maxValue - minValue)
        self.cluster(copyCell)
    
    def cluster(self, dataFrame):
        '''
        Uses the sklearn module to perform PCA on normalized data. Input: dataframe,
        output dataframe with PC values. 
        '''
        df = dataFrame
        genes = list(df.columns.values)
        genes = genes[1:]
        
        x = df.loc[:, genes].values
        
        pca = PCA(n_components=2)
        
        components = pca.fit_transform(x) #performs pca on genes in dataset
        
        principalDf = pd.DataFrame(data = components
             , columns = ['pc1', 'pc2'])
        self.label(principalDf)
    
    
    def label(self, dataFrame):
        '''
        Performs agglomerative hierarchal clustering using sklearn package to identify 
        and label clustersin PCA data. Number of clusters found is set to three. 
        Likely will change to be an option. Plots clusters.
        '''
        df = dataFrame
        
        hc = AgglomerativeClustering(n_clusters=3, affinity = 'euclidean', linkage = 'ward')
        
        y_df = hc.fit_predict(df)
        
        df['label'] = pd.Series(y_df, index = df.index)
        df.label = df.label.astype(str)
        
        df['stem'] = self.stemData.iloc[:,1]
            
        self.linear(df)   
        self.plotCluster(df)
        
    def linear(self, dataFrame):
        '''
        Uses sklearn package to perform linear modeling on clustered data vs input variable
        , and prints R squared value and plots results. 
        '''
        data = dataFrame
        
        x = data[['pc1', 'pc2']]
        y = data[['stem']]
        
        X_train, X_test, y_train, y_test = train_test_split(x, y, test_size = 0.2)
        
        lm = linear_model.LinearRegression()
        model = lm.fit(X_train,y_train)
        
        predictions = lm.predict(X_test)
        
        Rsq = lm.score(x,y)
        print(model.score(X_test, y_test))
        self.plotLinear(y_test, predictions, data)
        
    
    def plotCluster(self, dataFrame):
        '''
        Function for plotting clustered data.
        '''
        finalDf = dataFrame
        
        fig = plt.figure(figsize = (8,8)) 
        ax = fig.add_subplot(111, projection='3d')
        ax.set_xlabel('Principal Component 1', fontsize = 15)
        ax.set_ylabel('Principal Component 2', fontsize = 15)
        ax.set_zlabel('Stemness Rating', fontsize = 15)
        ax.set_title('2 component PCA with Agglomerative Hierarchical Clustering', fontsize = 20)
        labels = ['0', '1', '2']
        colors = ['r', 'g', 'b']
        for label, color in zip(labels,colors):
            indicesToKeep = finalDf['label'] == label
            ax.scatter(finalDf.loc[indicesToKeep, 'pc1']
               , finalDf.loc[indicesToKeep, 'pc2']
               , finalDf.loc[indicesToKeep, 'stem']
               , c = color
               , s = 50)
        ax.grid()
        
        fig.savefig(r"C:\Users\dlenc\Downloads\FinalProjectData\pca-Agg2.png")
    
    def plotLinear(self, test, pred, df):
        '''
        Function for plotting modeled data.
        '''
        y_test = test
        finaldf = df
        predictions = pred
        
        fig = plt.figure(figsize = (8,8))
        ax = fig.add_subplot(1,1,1)
        ax.set_xlabel('Predicted Stemness', fontsize = 15)
        ax.set_ylabel('Actual Stemness', fontsize = 15)
        ax.set_title('Predicted Values Using Linear Model', fontsize = 20)
        ax.scatter(y_test
            , predictions
            , c = 'b'
            , s = 50)
        ax.grid()
        
        fig.savefig(r"C:\Users\dlenc\Downloads\FinalProjectData\model.png")
        
def main():
    #Read csv file and transfer into data frame using pandas
    #Call findCells class with data frame
    file1 = input('Gene data:')
    file2 = input('Relevant Genes:')
    file3 = input('Variable:')
    genes =  set(open(file2).read().split())
    cellData=feather.read_dataframe(file1)#,delimiter='\t',encoding='utf-8')
    stemData=feather.read_dataframe(file3)
    findCells(cellData, genes, stemData)
    
if __name__ == "__main__":
    main()  
