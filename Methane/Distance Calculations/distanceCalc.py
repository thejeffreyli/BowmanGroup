import numpy as np
import itertools as it
import csv  
import matplotlib.pyplot as plt  


def read(dataSet): #reads the xyz file, returns lines
    xyz = open(dataSet, "r") #open file for glycine data
    lines = xyz.readlines()
    xyz.close()
    return lines #lines is an array

def remove(lines): #removes non numerical data lines
    for i in range(len(lines)):  #marks the first two lines in every 12 
        if (i % 7 == 0 or i % 7 == 1):
            lines[i] = "holder"        
    lines[:] = [x for x in lines if "holder" not in x] #removes holders  
    return lines #array of pure numerical data

def splitData(lines): #split the list accordingly: 4 groups
    atoms = []
    coordinates = []
    
    for line in lines: 
        atom,x,y,z = line.split() #  4 groups
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
    for i in range(0, len(lines), 5): #atom 1: H1
        arrayi.append(i)
    for j in range(1, len(lines), 5): #atom 2: H2      
        arrayj.append(j)
    for k in range(2, len(lines), 5): #atom 3: H3               
        arrayk.append(k)     
    for l in range(3, len(lines), 5): #atom 4: H4              
        arrayl.append(l)   
    for m in range(4, len(lines), 5): #atom 5: C              
        arraym.append(m)          
     
    # 0, 1, 2, 3, 4 
    # H, H, H, H, C      
    incArr = []
    for i in range(len(arrayi)):
        incArr.append(arrayi[i])
        incArr.append(arrayj[i])
        incArr.append(arrayk[i])
        incArr.append(arrayl[i])        
        incArr.append(arraym[i])        

    length = len(incArr) 
    length = length/5
    incArr = np.array_split(incArr,(length)) # number of arrays
    
    #   HH     HH     HH     CH     HH     HH     CH    HH     CH       CH
    # [r 12 , r 13 , r 14 , r 15 , r 23 , r 24 , r 25 , r 34 , r 35 , r 45 ].
    distArr = []
    for i in range(len(incArr)):
        for a, b in it.combinations(incArr[i], 2):
            distArr.append(euclidDist(coordinates[a], coordinates[b]))
            
    distArr = np.array_split(distArr,(length))
    return distArr
      
#removes CH 
def onlyHH(vecData):
    

    for i in range(len(vecData)):
        vecData[i] = np.delete(vecData[i], 9)
        vecData[i] = np.delete(vecData[i], 8)
        vecData[i] = np.delete(vecData[i], 6)
        vecData[i] = np.delete(vecData[i], 3)
    

    return vecData
 
#removes HH            
def onlyCH(vecData):
    
    for i in range(len(vecData)):
        vecData[i] = np.delete(vecData[i], 7)
        vecData[i] = np.delete(vecData[i], 5)
        vecData[i] = np.delete(vecData[i], 4)
        vecData[i] = np.delete(vecData[i], 2)
        vecData[i] = np.delete(vecData[i], 1)
        vecData[i] = np.delete(vecData[i], 0)
        
    return vecData

# seeing which HH is the shortest
def countHH(indexList):
    count1, count2, count3, count4, count5, count6 = (0 for i in range(6))
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
        else:
            count6+=1      
                        
    return count1, count2, count3, count4, count5, count6    

# seeing which HH is the shortest
def countCH(indexList):
    count1, count2, count3, count4 = (0 for i in range(4))
    for i in range(len(indexList)):
        if indexList[i][0] == 1:
            count1+=1
        elif indexList[i][0] == 2:
            count2+=1
        elif indexList[i][0] == 3:
            count3+=1                    
        else:
            count4+=1      
                        
    return count1, count2, count3, count4  


def main():
    
    
    data = "err.xyz" #big data set
       
    linesData = read(data)
    # print(len(linesData)) #1743

    linesData = remove(linesData)
    # print(len(linesData)) #1245

    atomsData, coordData = splitData(linesData)
    # print(len(coordData)) # 1245  
    
    vecData = distVec(coordData, linesData) #vectors 
    # print(len(vecData)) #249 holes   
    vecDataCopy = distVec(coordData, linesData) #vectors (same as original)
    # print(len(vecDataCopy)) #249 holes       

    #compiling csv file 
    with open("distVec.csv", 'w') as csvfile:  
        csvwriter = csv.writer(csvfile)   
        csvwriter.writerows(vecData)   

    
    "--------------------------Determining Shortest HH--------------------------"       

    vecDataHH = onlyHH(vecData)    
    # vecDataHH = vecDataHH / equiHH 
    
    # print(type(vecDataHH))
    
    #compiling csv file 
    with open("HH_DistVec.csv", 'w') as csvfile:  
        csvwriter = csv.writer(csvfile)   
        csvwriter.writerows(vecDataHH)      
    
    sortedMagn = []
    for i in range(len(vecDataHH)):
        sortedMagn.append(sorted(vecDataHH[i]))     

    #compiling csv file 
    with open("HH_SortedDist.csv", 'w') as csvfile:  
        csvwriter = csv.writer(csvfile)  
        csvwriter.writerows(sortedMagn)            
        
    sortedIdx = []
    for i in range(len(vecDataHH)):
        x = np.argsort(vecDataHH[i])
        sortedIdx.append(x)   
         
    sortedIdx = [x+1 for x in sortedIdx]               
        
    #compiling csv file 
    with open("HH_SortedIdx.csv", 'w') as csvfile:  
        csvwriter = csv.writer(csvfile)  
        csvwriter.writerows(sortedIdx) 
        
    c1, c2, c3, c4, c5, c6 = countHH(sortedIdx)        
    
    fig = plt.figure()
    ax = fig.add_axes([0,0,1,1])
    confPlot = ['H1H2', 'H1H3', 'H1H4', 'H2H3', 'H2H4', 'H3H4']
    countPlot = [c1, c2, c3, c4, c5, c6]
    ax.bar(confPlot, countPlot)
    ax.set_ylabel('Number of Shortest Distances')
    ax.set_xlabel('Particular HH')
    ax.set_title('Determination of Shortest HH Distances')
    plt.show()    
        
    "--------------------------Determining Shortest CH--------------------------"         
        
    vecDataCH = onlyCH(vecDataCopy)  
    # vecDataCH = vecDataCH / equiCH 
        
    #compiling csv file 
    with open("CH_DistVec.csv", 'w') as csvfile:  
        csvwriter = csv.writer(csvfile)   
        csvwriter.writerows(vecDataCH)      
    
    sortedMagn = []
    for i in range(len(vecDataCH)):
        sortedMagn.append(sorted(vecDataCH[i]))     

    #compiling csv file 
    with open("CH_SortedDist.csv", 'w') as csvfile:  
        csvwriter = csv.writer(csvfile)  
        csvwriter.writerows(sortedMagn)            
        
    sortedIdx = []
    for i in range(len(vecDataCH)):
        x = np.argsort(vecDataCH[i])
        sortedIdx.append(x)   
         
    sortedIdx = [x+1 for x in sortedIdx]          
        
    #compiling csv file 
    with open("CH_SortedIdx.csv", 'w') as csvfile:  
        csvwriter = csv.writer(csvfile)  
        csvwriter.writerows(sortedIdx)     

    c1, c2, c3, c4  = countCH(sortedIdx)        
    
    fig = plt.figure()
    ax = fig.add_axes([0,0,1,1])
    confPlot = ['CH1', 'CH2', 'CH3', 'CH4']
    countPlot = [c1, c2, c3, c4]
    ax.bar(confPlot, countPlot)
    ax.set_ylabel('Number of Shortest Distances')
    ax.set_xlabel('Particular CH')
    ax.set_title('Determination of Shortest CH Distances')
    plt.show()    
            
    
if __name__ == "__main__":
    main()      