import numpy as np
import matplotlib.pyplot as plt  
import pandas as pd
from statistics import median


def read(dataSet): #reads the xyz file, returns lines
    xyz = open(dataSet, "r") #open file 
    lines = xyz.readlines()[1:] #skips first line
    xyz.close()
    return lines #lines is an array

def main():
    
    "--------------------------Histogram CH4 1.0--------------------------"    
    
    data = "ch4_1.0.txt" #big data set
    df = pd.read_csv(data)
    hist = df.hist(bins=30)
    plt.title("Histogram CH4 1.0")
    plt.xlabel('Potentials (cm^-1)')
    plt.ylabel('Frequency')    
    plt.xlim(0, 40000)    
    # plt.ylim(0, 40000)     
    
    df = df.to_numpy()
  
    print("number of energy values: ", 197434 )
    print("maximum:" , max(df))    
    print("minimum:" , min(df))
    print("median:" , median(df))       
     
    "--------------------------Histogram CH4 0.5--------------------------"    
    
    data = "ch4_0.5.txt" #big data set
    df = pd.read_csv(data)
    hist = df.hist(bins=30)
    plt.title("Histogram CH4 0.5")
    plt.xlabel('Potentials (cm^-1)')
    plt.ylabel('Frequency')    
    plt.xlim(0, 40000)    
    # plt.ylim(0, 40000)     
    
    df = df.to_numpy()
    
    print("number of energy values: ", 195578 )    
    print("maximum:" , max(df))    
    print("minimum:" , min(df))
    print("median:" , median(df))     
    
if __name__ == "__main__":
    main()    
 