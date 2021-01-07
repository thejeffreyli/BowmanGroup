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


#calculating distance between all atoms
def distVec(coordinates, lines):
    #arrays for increments
    N, H1, H2, C1, H3, H4, C2, O1, O2, H5 = ([] for i in range(10))

    #increments 
    for i in range(0, len(lines), 10): #atom 1: N
        N.append(i)
    for i in range(1, len(lines), 10): #atom 2: H
        H1.append(i)
    for i in range(2, len(lines), 10): #atom 3: H
        H2.append(i)
    for i in range(3, len(lines), 10): #atom 4: C
        C1.append(i)
    for i in range(4, len(lines), 10): #atom 5: H
        H3.append(i)        
    for i in range(5, len(lines), 10): #atom 6: H
        H4.append(i)        
    for i in range(6, len(lines), 10): #atom 7: C
        C2.append(i)
    for i in range(7, len(lines), 10): #atom 8: O
        O1.append(i)
    for i in range(8, len(lines), 10): #atom 9: O
        O2.append(i)
    for i in range(9, len(lines), 10): #atom 10: H
        H5.append(i)

    # 0, 1, 2, 3, 4, 5, 6, 7, 8, 9
    # N, H, H, C, H, H, C, O, O, H      
    incArr = []
    for i in range(len(N)):
        incArr.append(N[i])
        incArr.append(H1[i])
        incArr.append(H2[i])
        incArr.append(C1[i])        
        incArr.append(H3[i])     
        incArr.append(H4[i])
        incArr.append(C2[i])
        incArr.append(O1[i])
        incArr.append(O2[i])        
        incArr.append(H5[i])           

    length = len(incArr) # 52334 arrays with 10  atoms
    length = length/10
    incArr = np.array_split(incArr, length) 

    distArr = []
    for i in range(len(incArr)):
        for a, b in it.combinations(incArr[i], 2):     # 10C2 = 45 distances between 10 atoms
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
    df.to_csv("Conformer with H Matrix.csv")


if __name__ == "__main__":
    main()      