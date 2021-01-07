import numpy as np
import itertools as it
import matplotlib.pyplot as plt   
from sklearn.preprocessing import MinMaxScaler


def read(dataSet): #reads the xyz file, returns lines
    xyz = open(dataSet, "r") #open file for glycine data
    lines = xyz.readlines()
    xyz.close()
    return lines #lines is an array

def extractEnergy(lines): #removes non numerical data lines
    temp = []
    for i in range(len(lines)): 
        if (i % 12 != 1 ):
            lines[i] = "holder"        
    lines[:] = [x for x in lines if "holder" not in x] #removes holders  
    return lines #array of pure numerical data

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
    linesCoordV = read(glyCoord) #vectors
    linesCoordE = read(glyCoord) # energies
    
    
    linesCoordE = extractEnergy(linesCoordE) # y axis
    
    coordEnergy = []
    for i in range(len(linesCoordE)):
        coordEnergy.append(float(linesCoordE[i]))



    linesCoord = remove(linesCoordV)    
    atomsCoord, coordCoord = splitCoord(linesCoord)
    vecCoord = distVec(coordCoord, linesCoord) #vectors for small data
    

    
    sumDist = []
    for i in range(len(vecCoord)):
        sumDist.append(np.sum(vecCoord[i]))
        
        
    #sumDist = MinMaxScaler().fit_transform(np.asarray(sumDist).reshape(-1, 1))        

 
    plt.scatter(sumDist, coordEnergy)
    plt.xlabel("Sum of All Intermolecular Distances (units)")
    plt.ylabel("Energy (units) ")
    plt.title("Scatterplot of Eight Glycine Conformers")
    plt.ylim(-284.5, -283.5)
    plt.xlim(19, 27)  
    
    for x, y in zip(sumDist, coordEnergy):
        rgb = (np.random.random(), np.random.random(), np.random.random())
        plt.scatter(x, y, c=[rgb])
    
    plt.legend('12345678')
    
    plt.show()    
if __name__ == "__main__":
    main()      
