import numpy as np
import pandas as pd
import matplotlib.pyplot as plt   
import argparse

# read an input file and convert it to numpy
def toNumPy(filename):
    df = pd.read_csv(filename)
    return df.to_numpy()

def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument("sortedIndex")
    parser.add_argument("sortedResults")
    args = parser.parse_args()
    
    # load data
    sortedIndex = toNumPy(args.sortedIndex)
    sortedResults = toNumPy(args.sortedResults)
    a12, b12, a13, b13, a14, b14, a15, b15, a16, b16, a17, b17, a18, b18 = ([] for i in range(14))
 
    for i in range(len(sortedIndex)):
        
        # if first index = conformer 1
        if sortedIndex[i][0] == 1: 
            if sortedIndex[i][1] == 2:
                a12.append(sortedResults[i][0]) # 0
                b12.append(sortedResults[i][1])
            elif sortedIndex[i][1] == 3:
                a13.append(sortedResults[i][0]) # 0
                b13.append(sortedResults[i][1])               
            elif sortedIndex[i][1] == 4:
                a14.append(sortedResults[i][0]) # 2935
                b14.append(sortedResults[i][1])               
            elif sortedIndex[i][1] == 5:
                a15.append(sortedResults[i][0]) # 5
                b15.append(sortedResults[i][1])               
            elif sortedIndex[i][1] == 6:
                a16.append(sortedResults[i][0]) # 9833
                b16.append(sortedResults[i][1])              
            elif sortedIndex[i][1] == 7:  
                a17.append(sortedResults[i][0]) # 220
                b17.append(sortedResults[i][1])
            else:
                a18.append(sortedResults[i][0]) # 24
                b18.append(sortedResults[i][1]) 
                
    'for heavy atoms'                
    # fig, ax = plt.subplots()
    # plt.scatter(a16, b16, s = 10, c = 'blue', alpha = 0.7) # Conformer 6
    # plt.scatter(a14, b14, s = 10, c = 'green', alpha = 0.65) # Conformer 4
    # plt.scatter(a17, b17, s = 10, c = 'orange', alpha = 0.60) # Conformer 7
    # plt.scatter(a18, b18, s = 10, c = "purple", alpha = 1) # Conformer 8
    # plt.scatter(a15, b15, s = 10, c = 'red', alpha = 1) # Conformer 5
    # plt.scatter(a12, b12, s = 10, alpha = 0.7) # Conformer 2
    # plt.scatter(a13, b13, s = 10, alpha = 0.7) # Conformer 3
    

    
    # plt.xlabel('Conformer 1')
    # plt.ylabel('Conformer X')
    # plt.xlim(-0.1, 2.5)    
    # plt.ylim(-0.1, 2.5) 

    # count for heavy atoms
    # plt.legend(["Conf. 6: 2935", "Conf. 4: 2935", "Conf. 7: 220", "Conf. 8: 24", "Conf. 5: 5", "Conf. 2: 0", "Conf. 3: 0"])   
    
    
    'with H'         
    fig, ax = plt.subplots()
    plt.scatter(a14, b14, s = 5, c = 'green', alpha = 0.7) # Conformer 4
    plt.scatter(a16, b16, s = 5, c = 'blue', alpha = 0.65) # Conformer 6
    plt.scatter(a13, b13, s = 5, c = 'orange', alpha = 0.6) # Conformer 3
    plt.scatter(a15, b15, s = 5, c = 'red', alpha = 1) # Conformer 5 
    plt.scatter(a18, b18, s = 5, c = "purple", alpha = 1) # Conformer 8    
    plt.scatter(a17, b17, s = 10, alpha = 0.60) # Conformer 7    
    plt.scatter(a12, b12, s = 10, alpha = 0.7) # Conformer 2

    
    plt.xlabel('Conformer 1')
    plt.ylabel('Conformer X')
    # plt.xlim(15, 17.5)    
    # plt.ylim(12, 20) 
    
    
    
    plt.legend(["Conf. 4: 11295", "Conf. 6: 4596", "Conf. 3: 1455", "Conf. 5: 144", "Conf. 8: 2", "Conf. 7: 0", "Conf. 2: 0", ])       
    # plt.legend(["Conf. 4", "Conf. 6", "Conf. 3", "Conf. 5"])      
    #plt.legend(["Conf. 6"])   

    # Bar Graph
    
    # fig = plt.figure()
    # ax = fig.add_axes([0, 0, 1, 1])
    # conformers = ['Conf. 2', 'Conf. 3', 'Conf. 4', 'Conf. 5', 'Conf. 6', 'Conf. 7', 'Conf. 8']
    # count = [0, 0, 2935, 5, 9833, 220, 24]
    # ax.bar(conformers, count)


    plt.show()
    

    
if __name__ == "__main__":
    main()      