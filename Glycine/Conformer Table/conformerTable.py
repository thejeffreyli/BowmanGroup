import numpy as np
import itertools as it
import csv  
import matplotlib.pyplot as plt   
import pandas as pd

def read(dataSet): #reads the xyz file, returns lines
    xyz = open(dataSet, "r") #open file for glycine data
    lines = xyz.readlines()
    xyz.close()
    return lines #lines is an array

def remove(lines): #removes non numerical data lines
    for i in range(len(lines)):  #marks the first two lines in every 12 
        if (i % 12 == 0 or i % 12 == 1):
            lines[i] = "holder"        
    lines[:] = [x for x in lines if "holder" not in x] #removes holders  
    return lines #array of pure numerical data

    
def splitCoord(lines): #split the data accordingly: 4 groups    
    atoms = []
    coordinates = []
    
    for line in lines: 
        atom,x,y,z = line.split() #4 groups
        atoms.append(atom)
        coordinates.append([float(x), float(y), float(z)])
    return atoms, coordinates
       

#calculate distances (coordinates)   
def euclidDist(x1, x2):
       return np.sqrt(np.sum((np.array(x1) - np.array(x2))**2)) 



def compare(distanceA, distanceB):
    outcome = []
    for i in range(len(distanceA)):
        if distanceA[i] < distanceB[i]:
            outcome.append(-1)
        elif distanceA[i] > distanceB[i]:
            outcome.append(1)
        else: 
            outcome.append(0)
    return outcome


#calculating distance between heavy atoms
def distVec(coordinates, lines):
    #arrays for increments
    arrayi, arrayj, arrayk, arrayl, arraym = ([] for i in range(5))

    #increments 
    for i in range(0, len(lines), 10): #atom 1: N
        arrayi.append(i)
    for j in range(3, len(lines), 10): #atom 4: C       
        arrayj.append(j)
    for k in range(6, len(lines), 10): #atom 7: C                
        arrayk.append(k)     
    for l in range(7, len(lines), 10): #atom 8: O               
        arrayl.append(l)   
    for m in range(8, len(lines), 10): #atom 9: O              
        arraym.append(m)          
     
    # 0, 3, 6, 7, 8 
    # N, C, C, O, O       
    incArr = []
    for i in range(len(arrayi)):
        incArr.append(arrayi[i])
        incArr.append(arrayj[i])
        incArr.append(arrayk[i])
        incArr.append(arrayl[i])        
        incArr.append(arraym[i])        

    length = len(incArr) # 52334 arrays with 5 heavy atoms
    length = length/5
    incArr = np.array_split(incArr,(length)) 

    # [r 12 , r 13 , r 14 , r 15 , r 23 , r 24 , r 25 , r 34 , r 35 , r 45 ].
    distArr = []
    for i in range(len(incArr)):
        for a, b in it.combinations(incArr[i], 2):
            distArr.append(euclidDist(coordinates[a], coordinates[b]))
            
    distArr = np.array_split(distArr,(length))
    return distArr

            
def main():

    glyCoord = "GlyCoord.xyz" #small data set
    linesCoord = read(glyCoord)
    linesCoord = remove(linesCoord)
    atomsCoord, coordCoord = splitCoord(linesCoord)

    vecCoord = distVec(coordCoord, linesCoord) #vectors for small data
    vecCoord = np.array(vecCoord)

    # V - U, differences between V and U
    diff = []
    for v in range(len(vecCoord)): #iterate 1 to 8
        for u in range(len(vecCoord)): #iterate 1 to 8
            x = (np.subtract(vecCoord[v],vecCoord[u]))
            x = np.square(x)
            diff.append(x)
            
    # summing up the vector differences
    magn = []
    for i in range(len(diff)):
        x = np.sum(diff[i])
        magn.append(x)      
    #print(magn)
  
    matrix = []
    while magn != []:
        matrix.append(magn[:8])
        magn = magn[8:]
    
    
    #print (matrix)
    df = pd.DataFrame(data=matrix)  
    print(df)
    df.index += 1 
    df.columns += 1     
    # df.to_csv("Conformer Matrix.csv")


if __name__ == "__main__":
    main()      