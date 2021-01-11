import numpy as np
import itertools as it
import csv  
import matplotlib.pyplot as plt  
import pandas as pd

equiCH = 1.09455
equiHH = 1.78739


def file_to_numpy(filename):
    return np.genfromtxt(filename, delimiter=',')

def read(dataSet): #reads the xyz file, returns lines
    xyz = open(dataSet, "r") #open file for glycine data
    lines = xyz.readlines()
    xyz.close()
    return lines #lines is an array

def main():
    
    "--------------------------Normalize--------------------------"   
    
    CH_data = file_to_numpy("CH_SortedDist.csv")
    CH_data = CH_data / equiCH # normalize
    
    #compiling csv file 
    with open("Normalized_CH_SortedDist.csv", 'w') as csvfile:  
        csvwriter = csv.writer(csvfile)   
        csvwriter.writerows(CH_data)   


    HH_data = file_to_numpy("HH_SortedDist.csv")    
    HH_data = HH_data / equiHH # normalize

    #compiling csv file 
    with open("Normalized_HH_SortedDist.csv", 'w') as csvfile:  
        csvwriter = csv.writer(csvfile)   
        csvwriter.writerows(HH_data)   

    "--------------------------Plotting--------------------------"  
    
    xArr, yArr = [], [] #249 holes
    
    for i in range(len(HH_data)): # X Axis: HH shortest
        xArr.append(HH_data[i][0])
    for i in range(len(CH_data)): # Y Axis: CH shortest
        yArr.append(CH_data[i][0])    
    # for i in range(len(CH_data)): # Y Axis: CH longest
    #     yArr.append(CH_data[i][0])           
        
    fig, ax = plt.subplots()
    plt.scatter(xArr, yArr, s = 10, c = 'blue', alpha = 1) 
    
    plt.xlabel('Shortest HH')
    plt.ylabel('Shortest CH')
    # plt.xlim(0.25, 1.5)    
    # plt.ylim(0.25, 1.5)
    plt.axhline(y=1, linewidth=1, color='r')
    plt.title("Shortest CH v. Shortest HH")

    "--------------------------CH distance (Y) > 1--------------------------"  

    
    # print(yArr)
    indexes = [index for index, value in enumerate(yArr) if value > 1]
    print(indexes) # 40 holes with CH > 1.0 
      

     


    

      
     
if __name__ == "__main__":
    main()    
 