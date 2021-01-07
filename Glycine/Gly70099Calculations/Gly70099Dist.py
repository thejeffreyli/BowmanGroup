import numpy as np
import itertools as it
import csv  
import matplotlib.pyplot as plt  
import scipy.linalg as la 

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

def splitData(lines): #split the list accordingly: 7 groups
    atoms = []
    coordinates = []
    
    for line in lines: 
        atom,x,y,z,g1,g2,g3 = line.split() #  7 groups
        atoms.append(atom)
        coordinates.append([float(x), float(y), float(z)])
    return atoms, coordinates
    
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

#checking the order of two O atoms
def distanceTest(coordinates, lines):

    #arrays for increments
    arrayi = [] 
    arrayj = []
    arrayk = []

    #increments 
    for i in range(7, len(lines), 10): #atom 8
        arrayi.append(i)
    for j in range(8, len(lines), 10): #atom 9       
        arrayj.append(j)
    for k in range(9, len(lines), 10): #atom 10                
        arrayk.append(k)  

    #distance arrays
    dist810 = []
    dist910 = []
    
    #find distance between atom 8 and 10
    for x in range(len(arrayi)): 
        dist810.append(euclidDist(coordinates[arrayi[x]], coordinates[arrayk[x]]))
    
    #find distance between atom 9 and 10
    for x in range(len(arrayj)): 
        dist910.append(euclidDist(coordinates[arrayj[x]], coordinates[arrayk[x]]))

    return dist810, dist910

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
      
# iterates through the first index of each array, 
# adding 1 to proper conformer count if it matches
def countConf(indexList):
    count1, count2, count3, count4, count5, count6, count7, count8 = (0 for i in range(8))
    for i in range(len(indexList)):
        if indexList[i][0] == 1:
            count1+=1
        elif indexList[i][0] == 2:
            count2+=1
        elif indexList[i][0] == 3:
            count3+=1                    
        elif indexList[i][0] == 4:
            count4+=1
        elif indexList[i][0] == 5:
            count5+=1
        elif indexList[i][0] == 6:
            count6+=1
        elif indexList[i][0] == 7:
            count7+=1     
        else:
            count8+=1                                 
    return count1, count2, count3, count4, count5, count6, count7, count8            
            
            
def main():
    glycData = "Gly70099.xyz" #big data set
    glyCoord = "GlyCoord.xyz" #small data set
    
        
    linesData = read(glycData)
    linesCoord = read(glyCoord)
    #print(len(linesData)) #841188
    
    linesData = remove(linesData)
    linesCoord = remove(linesCoord)
    #print(len(linesData)) #700990
    
    atomsData, coordData = splitData(linesData)
    atomsCoord, coordCoord = splitCoord(linesCoord)
    #print(len(coordData)) #700990  
    
    #checking the two O atoms
    dist810, dist910 = distanceTest(coordData, linesData)
    outcome = compare(dist810, dist910)
    #print(outcome) 
    #print(len(outcome)) #70099 
    ' Output is all 1s, indicating that the distance between '
    ' atoms 8 and 10 are larger than the distance between atoms '
    ' 9 and 10.'

    vecData = distVec(coordData, linesData) #vectors for big data
    print(vecData[0])
    #print(len(vecData)) #70099     
    
    vecCoord = distVec(coordCoord, linesCoord) #vectors for small data
    #print(len(vecCoord)) # 8 

    # V - U, differences between V and U
    diff = []
    for v in range(len(vecData)): #iterate 1 to 70099
        for u in range(len(vecCoord)): #iterate 1 to 8
            x = (np.subtract(vecData[v],vecCoord[u]))
            x = np.square(x)
            diff.append(x)
    #print(len(diff)) # 70099 * 8 = 560792
    
    # summing up the vector differences
    magn = []
    for i in range(len(diff)):
        x = np.sum(diff[i])
        magn.append(x)        
     
    # spitting the array into 70099 (k) groups of 8
    splitMagn = np.array_split(magn, len(vecData))
    
    
    # sorted every groups Vk - Uj 
    sortedMagn = []
    for i in range(len(splitMagn)):
        sortedMagn.append(sorted(splitMagn[i])) 
        
    #compiling csv file      
    with open("sortedResultsGly70099.csv", 'w') as csvfile:  
        #creating a csv writer object  
        csvwriter = csv.writer(csvfile)  
        # writing the data rows  
        csvwriter.writerows(sortedMagn)   
       
    ' Each of the 52334 arrays inside splitMagn is sorted '
    ' and compiled inside sortedMagn. '

    
    # sorted by index
    sortedIdx = []
    for i in range(len(splitMagn)):
        x = np.argsort(splitMagn[i])
        sortedIdx.append(x)   

    sortedIdx = [x+1 for x in sortedIdx]            
        
    #compiling csv file 
    with open("sortedIndexGly70099.csv", 'w') as csvfile:  
        #creating a csv writer object  
        csvwriter = csv.writer(csvfile)  
        # writing the data rows  
        csvwriter.writerows(sortedIdx)      

    ' Each of the 52334 arrays inside splitMagn is sorted by index'
    ' and compiled inside sortedIdx. '

    c1, c2, c3, c4, c5, c6, c7, c8 = countConf(sortedIdx)

    fig = plt.figure()
    ax = fig.add_axes([0,0,1,1])
    confPlot = ['Conf 1', 'Conf 2', 'Conf 3', 'Conf 4', 'Conf 5', 'Conf 6', 'Conf 7', 'Conf 8']
    countPlot = [c1, c2, c3, c4, c5, c6, c7, c8]
    ax.bar(confPlot, countPlot)
    ax.set_ylabel('Number of Configurations')
    ax.set_xlabel('Glycine Conformers')
    ax.set_title('Number of Configurations (Gly70099) for Each of the Eight Glycine Conformers')
    plt.show()
        

if __name__ == "__main__":
    main()      