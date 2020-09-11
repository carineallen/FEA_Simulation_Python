import numpy as np 
import math
import matplotlib.pyplot as plt 
import scipy.linalg as la

def InitializeComponents():
    
    #constructing model
    N_nodes = input("Enter the number of nodes: ")
    Nodes = np.zeros((int(N_nodes),3))
    for i in range(0,int(N_nodes)):
        X_NodePosition = input("Enter the X position of node " + str(i) + ": ")
        Y_NodePosition = input("Enter the Y position of node " + str(i) + ": ")
        IsFixed = input("Enter the degree of freedom of the node: ")
        Nodes[i,:] = [X_NodePosition,Y_NodePosition,IsFixed]
        
    N_elements = input("Enter the number of elements: ")
    Elements = np.zeros((int(N_elements),8))
    for j in range(0,int(N_elements)):
        LeftNode = input("Enter the first Node of element " + str(j) + ": ")
        RigthNode = input("Enter the second Node of element " +str(j) + ": ")
        YoungModulus = float(input("Enter Young Modulus in GPa: ")) * 10**9
        M_Inercia = float(input("Enter the seconds moment of inertia in meters m^4: "))
        Area = float(input("Enter Cross section area in meters: "))
        Angle = float(input("Enter the angle in degree: "))
        Vector = [Nodes[int(RigthNode),0] - Nodes[int(LeftNode),0], Nodes[int(RigthNode),1] - Nodes[int(LeftNode),1]]
        Length = np.linalg.norm(Vector)
        if (Length < 0):
            Length = Length*(-1)
        cosine = math.cos(Angle*np.pi / 180)
        sine = math.sin(Angle*np.pi /180)
        Elements[j,:] = [LeftNode,RigthNode,YoungModulus,Area,Length,cosine,sine,M_Inercia]
        
    #plot model
    plot_nodes(Nodes,N_nodes)
    
    for z in range(0,int(N_elements)):
        fromPoint = Elements[z,0]
        toPoint = Elements[z,1]
        area = Elements[z,3]
        LineColor = 'g'
        Draw_elements(fromPoint, toPoint, area,Nodes,LineColor)
        
    return {  'Nodes':Nodes, 'N_nodes':N_nodes,   \
	      	  'Elements':Elements, 'N_elements':N_elements}
         
#Function to plot the Nodes
def plot_nodes(Nodes,N_nodes):
    for i in range(0,int(N_nodes)):
         x = Nodes[i,0]
         y = Nodes[i,1]
         size = 400
         plt.scatter(x, y, c='b', s=size, zorder=5)
        
#Function to plot the elements
def Draw_elements(fromPoint, toPoint, area,Nodes,LineColor):
	x1 = Nodes[int(fromPoint),0]
	y1 = Nodes[int(fromPoint),1]
	x2 = Nodes[int(toPoint),0]
	y2 = Nodes[int(toPoint),1]
	plt.plot([x1, x2], [y1, y2], color= LineColor, linestyle='-', linewidth=1.5, zorder=1)


def CalculateMatrices(Parameters):
    
    elements    = Parameters['Elements']
    N_elements    = Parameters['N_elements']
    N_nodes    = Parameters['N_nodes']
    N_Forces = input("Enter the number of forces: ")
    Forces = np.zeros((int(N_Forces),4))
    for i in range(0,int(N_Forces)):
        F_NodePosition = input("Enter the node where the force " + str(i) + " will be applied: ")
        F_X = float(input("Enter the magnitude in X direction in N: "))
        F_Y = float(input("Enter the magnitude in Y direction in N: "))
        B_M = float(input("Enter the bending moment in Nm: "))
        Forces[i,:] = [F_NodePosition,F_X,F_Y,B_M]
        
    #construct striffness matrix
    K = np.zeros((int(N_nodes) * 3,int(N_nodes) * 3))
    KG = np.zeros((int(N_nodes) * 3,int(N_nodes) * 3))
    M = np.zeros((int(N_nodes) * 3,int(N_nodes) * 3))
    for j in range(0,int(N_elements)):
     E = elements[j,2]
     A = elements[j,3]
     L = elements[j,4]
     c = elements[j,5]
     s = elements[j,6]
     I = elements[j,7]
     #get matrix for the element
     T = np.array([[c,s,0,0,0,0],
                   [-s,c,0,0,0,0],
                   [0,0,1,0,0,0],
                   [0,0,0,c,s,0],
                   [0,0,0,-s,c,0],
                   [0,0,0,0,0,1]])
     st = 392064.6119
     
     T_Tranpose = T.transpose()
     
     #stiffness matrix
     K_1 = np.array([[E*A/L,0,0,-E*A/L,0,0],
                   [0,12*E*I/L**3,6*E*I/L**2,0,-12*E*I/L**3,6*E*I/L**2],
                   [0,6*E*I/L**2,4*E*I/L,0,-6*E*I/L**2,2*E*I/L],
                   [-E*A/L,0,0,E*A/L,0,0],
                   [0,-12*E*I/L**3,-6*E*I/L**2,0,12*E*I/L**3,-6*E*I/L**2],
                   [0,6*E*I/L**2,2*E*I/L,0,-6*E*I/L**2,4*E*I/L]])
     
     #geometric stiffness matrix
     K_G = np.array([[0,0,0,0,0,0],
                   [0,(36/30)/L,(1/10),0,-(36/30)/L,(1/10)],
                   [0,(1/10),(4/30)*L,0,-(1/10),-(1/30)*L],
                   [0,0,0,0,0,0],
                   [0,-(36/30)/L,-(1/10),0,(36/30)/L,-(1/10)],
                   [0,(1/10),-(1/30)*L,0,-(1/10),(4/30)*L]])
     
     #Mass matrix
     K_M = np.array([[(140/420)*(st*A)*L,0,0,(70/420)*(st*A)*L,0,0],
                    [0,(156/420)*(st*A)*L,(22/420)*(st*A)*L**2,0,(54/420)*(st*A)*L,-(13/420)*(st*A)*L**2],
                    [0,(22/420)*(st*A)*L**2,(4/420)*(st*A)*L**3,0,(13/420)*(st*A)*L**2,-(3/420)*(st*A)*L**3],
                    [(70/420)*(st*A)*L,0,0,(140/420)*(st*A)*L,0,0],
                    [0,(54/420)*(st*A)*L,(13/420)*(st*A)*L**2,0,(156/420)*(st*A)*L,-(22/420)*(st*A)*L**2],
                    [0,-(13/420)*(st*A)*L**2,-(3/420)*(st*A)*L**3,0,-(22/420)*(st*A)*L**2,(4/420)*(st*A)*L**3]])
     
     K_element = T_Tranpose.dot(K_1).dot(T)
     KG_element = T_Tranpose.dot(K_G).dot(T)
     #add to global matrix
     K[int(elements[j,0])*3:int(elements[j,1])*3+3,int(elements[j,0]*3):int(elements[j,1])*3+3] =K[int(elements[j,0])*3:int(elements[j,1])*3+3,int(elements[j,0]*3):int(elements[j,1])*3+3]+ K_element
     KG[int(elements[j,0])*3:int(elements[j,1])*3+3,int(elements[j,0]*3):int(elements[j,1])*3+3] =KG[int(elements[j,0])*3:int(elements[j,1])*3+3,int(elements[j,0]*3):int(elements[j,1])*3+3]+ KG_element 
     M[int(elements[j,0])*3:int(elements[j,1])*3+3,int(elements[j,0]*3):int(elements[j,1])*3+3] =M[int(elements[j,0])*3:int(elements[j,1])*3+3,int(elements[j,0]*3):int(elements[j,1])*3+3]+ K_M
     
    return { 'K':K, 'N_Forces':N_Forces, 'Forces':Forces, 'KG':KG, 'M':M}
     
#main function
def main():
    
    #Initialize the components
    Parameters = InitializeComponents()
    nodes    = Parameters['Nodes']
    N_nodes    = Parameters['N_nodes']
    elements    = Parameters['Elements']
    N_elements    = Parameters['N_elements']
    #Calculete the Stiffness Matrix and Force Vector
    Results = CalculateMatrices(Parameters)
    K    = Results['K']
    KG    = Results['KG']
    M    = Results['M']
    Forces    = Results['Forces']
    N_Forces    = Results['N_Forces']
    
    KC = K
    KG_F = KG
    M_F = M
    
    print("Global Stiffness matrix:")
    print("K = " + str(K) + "\n")
    
    #remove the Fixed Nodes
    N_Nodes_reemove = 0
    for d in range(0,int(N_nodes)):
            if nodes[d,2] == 0:
               N_Nodes_reemove = N_Nodes_reemove + 3
            elif nodes[d,2] == 1:
               N_Nodes_reemove = N_Nodes_reemove + 2
            elif nodes[d,2] == 2:
               N_Nodes_reemove = N_Nodes_reemove + 1
               
    remove = np.zeros((1,N_Nodes_reemove), dtype=int)
    r = 0
    for n in range(0,int(N_nodes)):
            if nodes[n,2] == 0:
               remove[0,r] = int(n*3)
               remove[0,r + 1] = int(n*3+1)
               remove[0,r + 2] = int(n*3+2)
               r = r + 3
            elif nodes[n,2] == 1:
               remove[0,r] = int(n*3)
               remove[0,r + 1] = int(n*3+1)
               r = r + 2
            elif nodes[n,2] == 2:
               remove[0,r] = int(n*3+1)
               r = r + 1
                        
    KC = np.delete(KC,remove,0)
    KC = np.delete(KC,remove,1)
    KG_F = np.delete(KG_F,remove,0)
    KG_F = np.delete(KG_F,remove,1)
    M_F = np.delete(M_F,remove,0)
    M_F = np.delete(M_F,remove,1)
    
    print("Contrained Stiffness matrix")
    print("K: " + str(KC) + "\n")

    #Calculate deformation for each free node
    deformations = np.zeros((int(N_nodes)*3,1))
    F = np.zeros((int(N_nodes)*3,1))
    for i in range(0,int(N_nodes)):
            ForceVector = np.zeros((3,1))
            for j in range(0,int(N_Forces)):
                ForcePosition = int(Forces[j,0])
                if ForcePosition == i:
                    ForceVector = ForceVector + np.array([[Forces[j,1]],[Forces[j,2]],[Forces[j,3]]])

                  
            F[i*3:i*3+3,0] = F[i*3:i*3+3,0] + ForceVector[0:3,0]
            
    F = np.delete(F,remove,0)
    try:
        D =np.linalg.lstsq(KC, F)
    except:
        D = F / KC[0,0]
     
    D = D[0]
    
    r = 0
    h = 0
    for i in range(0,int(N_nodes)):
        if nodes[i,2] == 3:
               nodes[i,0] = float(nodes[i,0]) + float( D[r,0])
               nodes[i,1] = float(nodes[i,1]) + float( D[r+1,0])
               deformations[h,0] = float( D[r,0])
               deformations[h+1,0] = float( D[r+1,0])
               deformations[h+2,0] = float( D[r+2,0])
               r = r + 3
               h = h + 3
        elif nodes[i,2] == 2:
               nodes[i,0] = float(nodes[i,0]) + float( D[r,0])
               deformations[h,0] = float( D[r,0])
               deformations[h+1,0] = 0
               deformations[h+2,0] = float( D[r+1,0])
               r = r + 2
               h = h + 3
        elif nodes[i,2] == 1:
            deformations[h,0] = 0
            deformations[h+1,0] = 0
            deformations[h+2,0] = float( D[r,0])
            r = r + 1
            h = h + 3
        elif nodes[i,2] == 0:
            deformations[h,0] = 0
            deformations[h+1,0] = 0
            deformations[h+2,0] = 0
            h = h + 3

     
    print("Deformation:")
    print(str(D) + "\n")       
    #plot element with new cordenates
    plot_nodes(nodes,N_nodes)
    for g in range(0,int(N_elements)):
        x1 = nodes[int(elements[g,0]),0] - nodes[int(elements[g,0]),0]
        x2 = nodes[int(elements[g,1]),0] - nodes[int(elements[g,0]),0]
        y1 = nodes[int(elements[g,0]),1]
        y2 = nodes[int(elements[g,1]),1]
        
        angle = math.acos(float(elements[g,5]))
        
        try:
            theta1 = angle + deformations[int(elements[g,0])*3+2,0]
        except:
            theta1 = 0
        try:
            theta2 = angle + deformations[int(elements[g,1])*3+2,0]
        except:
            theta2 = 0 
        
        
        CoefMat =np.array([[x1**3,x1**2,x1,1],
                           [x2**3,x2**2,x2,1],
                           [3*x1**2,2*x1,1,0],
                           [3*x2**2,2*x2,1,0]])
                    
        B=np.array([[y1],[y2],[theta1],[theta2]])
        z = np.linalg.lstsq(CoefMat,B)
        z_final = z[0]
        print("y(x) = " + str(float(z_final[0,0])) + "x^3 + (" + str(float(z_final[1,0])) + ")x^2 + (" + str(float(z_final[2,0])) + ")x + (" + str(float(z_final[3,0])) + ")\n" )
        print("M(x) = " + str(float(elements[g,2])*float(elements[g,7]) * 6 * float(z_final[0,0])) + "x + (" + str(float(elements[g,2])*float(elements[g,7]) * 2 * float(z_final[1,0])) + ")\n")
        print("coefficients:")
        print(str(z_final) + "\n")
        distance = x2 - x1
        X_axis = np.zeros((1,20))
        Y_axis = np.zeros((1,20))
        scale_factor = 50
        A = z_final[0,0]
        B = z_final[1,0]
        C = z_final[2,0]
        D = z_final[3,0] 
        for r in range(1,20):
            X_axis[0,r] = (r*(distance/20) + nodes[int(elements[g,0]),0])
        for b in range(1,20):
            Y_axis[0,b] = (A*(X_axis[0,b] - nodes[int(elements[g,0]),0])**3 + B*(X_axis[0,b]- nodes[int(elements[g,0]),0])**2 + C*(X_axis[0,b] - nodes[int(elements[g,0]),0]) +D)
        plt.plot(X_axis, Y_axis,'ro')
        
    
    eigvals, eigvecs = la.eig(KC,KG_F)
    eigvals = eigvals.real
    print("estimated critical loads: ")
    print(str(eigvals) + "\n")
    
    eigvals2, eigvecs2 = la.eig(KC,M_F)
    eigvals2 = eigvals2.real
    print("estimated eigfrequencies: ")
    print(str(eigvals2) + "\n")
    
    
if __name__ == '__main__':
	main()
        
        
        
