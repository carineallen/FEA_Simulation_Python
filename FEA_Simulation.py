import numpy as np 
import math
import matplotlib.pyplot as plt 

def InitializeComponents():
    
    #constructing model
    N_nodes = input("Enter the number of nodes: ")
    Nodes = np.zeros((int(N_nodes),3))
    for i in range(0,int(N_nodes)):
        X_NodePosition = input("Enter the X position of node " + str(i) + ": ")
        Y_NodePosition = input("Enter the Y position of node " + str(i) + ": ")
        IsFixed = input("This node it's fixed?(1 for Fixed, 0 for not fixed): ")
        Nodes[i,:] = [X_NodePosition,Y_NodePosition,IsFixed]
        #Nodes[i,1] = Y_NodePosition
        
    N_elements = input("Enter the number of elements: ")
    Elements = np.zeros((int(N_elements),7))
    for j in range(0,int(N_elements)):
        RigthNode = input("Enter the first Node of element " + str(j) + ": ")
        LeftNode = input("Enter the second Node of element " +str(j) + ": ")
        YoungModulus = float(input("Enter Young Modulus in GPa: ")) * 10**9
        Thickness = float(input("Enter Thickness in meters: "))
        Height = float(input("Enter Height in meters: "))
        Angle = float(input("Enter the angle in degree: "))
        Area = Thickness * Height
        Vector = [Nodes[int(RigthNode),0] - Nodes[int(LeftNode),0], Nodes[int(RigthNode),1] - Nodes[int(LeftNode),1]]
        Length = np.linalg.norm(Vector)
        cosine = math.cos(Angle*np.pi / 180)
        sine = math.sin(Angle*np.pi /180)
        Elements[j,:] = [RigthNode,LeftNode,YoungModulus,Area,Length,cosine,sine]
        
    plot_nodes(Nodes,N_nodes)
    
    for z in range(0,int(N_elements)):
        fromPoint = Elements[z,0]
        toPoint = Elements[z,1]
        area = Elements[z,3]
        LineColor = 'g'
        Draw_elements(fromPoint, toPoint, area,Nodes,LineColor)
        
    return {  'Nodes':Nodes, 'N_nodes':N_nodes,   \
	      	  'Elements':Elements, 'N_elements':N_elements}
         
def plot_nodes(Nodes,N_nodes):
    for i in range(0,int(N_nodes)):
         x = Nodes[i,0]
         y = Nodes[i,1]
         size = 400
         plt.scatter(x, y, c='b', s=size, zorder=5)
        
def Draw_elements(fromPoint, toPoint, area,Nodes,LineColor):
	x1 = Nodes[int(fromPoint),0]
	y1 = Nodes[int(fromPoint),1]
	x2 = Nodes[int(toPoint),0]
	y2 = Nodes[int(toPoint),1]
	plt.plot([x1, x2], [y1, y2], color= LineColor, linestyle='-', linewidth=200*area, zorder=1)
    
def CalculateMatrices(Parameters):
    
    nodes    = Parameters['Nodes']
    elements    = Parameters['Elements']
    N_elements    = Parameters['N_elements']
    N_nodes    = Parameters['N_nodes']
    N_Forces = input("Enter the number of forces: ")
    Forces = np.zeros((int(N_Forces),3))
    for i in range(0,int(N_Forces)):
        F_NodePosition = input("Enter the node where the force " + str(i) + " will be applied: ")
        F_X = float(input("Enter the magnitude in X direction in N: "))
        F_Y = float(input("Enter the magnitude in Y direction in N: "))
        Forces[i,:] = [F_NodePosition,F_X,F_Y]
        
    K = np.zeros((int(N_nodes) * 2,int(N_nodes) * 2))
    for j in range(0,int(N_elements)):
     E = elements[j,2]
     A = elements[j,3]
     L = elements[j,4]
     c = elements[j,5]
     s = elements[j,6]
     K_element = (E*A/L) * np.array([[c*c,c*s,-c*c,-c*s],
                                    [c*s,s*s,-c*s,-s*s],
                                    [-c*c,-c*s,c*c,c*s],
                                    [-c*s,-s*s,c*s,s*s]])
     
     K[int(elements[j,0])*2:int(elements[j,1])*2+2,int(elements[j,0]*2):int(elements[j,1])*2+2] =K[int(elements[j,0])*2:int(elements[j,1])*2+2,int(elements[j,0]*2):int(elements[j,1])*2+2]+ K_element
        
    return { 'K':K, 'N_Forces':N_Forces, 'Forces':Forces}
     
def main():
    
    Parameters = InitializeComponents()
    nodes    = Parameters['Nodes']
    N_nodes    = Parameters['N_nodes']
    elements    = Parameters['Elements']
    N_elements    = Parameters['N_elements']
    Results = CalculateMatrices(Parameters)
    K    = Results['K']
    Forces    = Results['Forces']
    N_Forces    = Results['N_Forces']
    
    for i in range(0,int(N_nodes)):
        if nodes[i,2] == 0:
            ForceVector = np.zeros((2,1))
            for j in range(0,int(N_Forces)):
                ForcePosition = int(Forces[j,0])
                if ForcePosition == i:
                    ForceVector = ForceVector + np.array([[Forces[j,1]],[Forces[j,2]]])
            KC = K[i*2:i*2+2,i*2:i*2+2]
            deformation = np.linalg.inv(KC).dot(ForceVector)
            print(deformation)
            nodes[i,0] = float(nodes[i,0]) + float(deformation[0,0]) * 200
            nodes[i,1] = float(nodes[i,1]) + float(deformation[1,0]) *200
    
    plot_nodes(nodes,N_nodes)
    for z in range(0,int(N_elements)):
        fromPoint = elements[z,0]
        toPoint = elements[z,1]
        area = elements[z,3]
        LineColor = 'r'
        Draw_elements(fromPoint, toPoint, area,nodes, LineColor)
    plt.title('Analysis of Truss Structure')
    plt.show()

if __name__ == '__main__':
	main()
        
        
        
