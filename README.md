# FEA_Simulation_Python_2

### Example Inputs:
Enter the number of nodes: 2

Enter the X position of node 0: 0

Enter the Y position of node 0: 0

Enter the degree of freedom of the node: 0

Enter the X position of node 1: 2

Enter the Y position of node 1: 0

Enter the degree of freedom of the node: 3

Enter the number of elements: 1

Enter the first Node of element 0: 0

Enter the second Node of element 0: 1

Enter Young Modulus in GPa: 210

Enter the seconds moment of inertia in meters m^4: 0.000002886667

Enter Cross section area in meters: 0.0026

Enter the angle in degree: 0

Enter the number of forces: 1

Enter the node where the force 0 will be applied: 1

Enter the magnitude in X direction in N: 0

Enter the magnitude in Y direction in N: -10000

Enter the bending moment in Nm: 0

### Example Output:
Global Stiffness matrix:
K = [[ 2.73000000e+08  0.00000000e+00  0.00000000e+00 -2.73000000e+08
   0.00000000e+00  0.00000000e+00]
 [ 0.00000000e+00  9.09300105e+05  9.09300105e+05  0.00000000e+00
  -9.09300105e+05  9.09300105e+05]
 [ 0.00000000e+00  9.09300105e+05  1.21240014e+06  0.00000000e+00
  -9.09300105e+05  6.06200070e+05]
 [-2.73000000e+08  0.00000000e+00  0.00000000e+00  2.73000000e+08
   0.00000000e+00  0.00000000e+00]
 [ 0.00000000e+00 -9.09300105e+05 -9.09300105e+05  0.00000000e+00
   9.09300105e+05 -9.09300105e+05]
 [ 0.00000000e+00  9.09300105e+05  6.06200070e+05  0.00000000e+00
  -9.09300105e+05  1.21240014e+06]]

Contrained Stiffness matrix
K: [[ 2.73000000e+08  0.00000000e+00  0.00000000e+00]
 [ 0.00000000e+00  9.09300105e+05 -9.09300105e+05]
 [ 0.00000000e+00 -9.09300105e+05  1.21240014e+06]]

Deformation:
[[ 0.        ]
 [-0.04398988]
 [-0.03299241]]

y(x) = 0.0027493673279626422x^3 + (-0.016496203967775843)x^2 + (1.819632721480302e-18)x + (1.46455005105198e-17)

M(x) = 10000.0x + (-19999.999999999985)

coefficients:
[[ 2.74936733e-03]
 [-1.64962040e-02]
 [ 1.81963272e-18]
 [ 1.46455005e-17]]

estimated critical loads:
[4876986.40099404  376747.53900596              inf]

estimated eigfrequencies: 
[ 45029.33428177    463.85931188 401719.50035667]

![Output result](https://github.com/carineallen/FEA_Simulation_Python/blob/master/exercise%2010.3.3.PNG)

# FEA_Simulation_Python
FEA Simulation for simple Beams 

## Example for 2 degrees of freedom

### Example Inputs:

Enter the degree of freedom of the nodes: 2

Enter the number of nodes: 3

Enter the X position of node 0: 0

Enter the Y position of node 0: 0

This node it's fixed?(1 for Fixed, 0 for not fixed): 1

Enter the X position of node 1: 2

Enter the Y position of node 1: 0

This node it's fixed?(1 for Fixed, 0 for not fixed): 0

Enter the X position of node 2: 0

Enter the Y position of node 2: 3.64

This node it's fixed?(1 for Fixed, 0 for not fixed): 1

Enter the number of elements: 2

Enter the first Node of element 0: 0

Enter the second Node of element 0: 1

Enter Young Modulus in GPa: 210

Enter Thickness in meters: 0.05

Enter Height in meters: 0.1

Enter the angle in degree: 0

Enter the first Node of element 1: 1

Enter the second Node of element 1: 2

Enter Young Modulus in GPa: 210

Enter Thickness in meters: 0.05

Enter Height in meters: 0.1

Enter the angle in degree: -60

Enter the number of forces: 1

Enter the node where the force 0 will be applied: 1

Enter the magnitude in X direction in N: 0

Enter the magnitude in Y direction in N: -100000

### Example Output:

[[ 5.88203306e+08 -1.09471337e+08]
 [-1.09471337e+08  1.89609917e+08]]
 
[[-0.00010997]
 [-0.00059089]] 

![Output result](https://github.com/carineallen/FEA_Simulation_Python/blob/master/FEA_Simulation_python2.png)

## Example for 1 degrees of freedom

### Example Inputs:

Enter the degree of freedom of the nodes: 1

Enter the number of nodes: 3

Enter the X position of node 0: 0

Enter the Y position of node 0: 0

This node it's fixed?(1 for Fixed, 0 for not fixed): 1

Enter the X position of node 1: 0.5

Enter the Y position of node 1: 0

This node it's fixed?(1 for Fixed, 0 for not fixed): 0

Enter the X position of node 2: 1

Enter the Y position of node 2: 0

This node it's fixed?(1 for Fixed, 0 for not fixed): 0

Enter the number of elements: 2

Enter the first Node of element 0: 0

Enter the second Node of element 0: 1

Enter Young Modulus in GPa: 210

Enter Thickness in meters: 0.01131

Enter Height in meters: 1

Enter the angle in degree: 0

Enter the first Node of element 1: 1

Enter the second Node of element 1: 2

Enter Young Modulus in GPa: 210

Enter Thickness in meters: 0.007854

Enter Height in meters: 1

Enter the angle in degree: 0

Enter the number of forces: 2

Enter the node where the force 0 will be applied: 1

Enter the magnitude in X direction in N: 300000

Enter the magnitude in Y direction in N: 0

Enter the node where the force 1 will be applied: 2

Enter the magnitude in X direction in N: 200000

Enter the magnitude in Y direction in N: 0

### Example Output:

[[ 8.04888e+09 -3.29868e+09]
 [-3.29868e+09  3.29868e+09]]

[[0.00010526 0.        ]
 [0.00016589 0.        ]]
 
![Output result](https://github.com/carineallen/FEA_Simulation_Python/blob/master/FEA_Simulation_python10.2.1.png)
