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
