# Explanation

## What you need to install to use these codes?
This repository contains MATLAB files used to generate the results in the paper:

"Simultaneous control and control node selection for nonlinear network  
dynamics"

Independently from this paper and its results, the codes can be used to optimally select control nodes of nonlinear networks whose dynamics is smooth and can be described by 

$$\dot{\mathbf{x}}=\boldsymbol{f}(\mathbf{x})+B\mathbf{u}$$

where $\mathbf{x}$ is the network state, $B$ is the control input matrix, and $\mathbf{u}$ is the control input sequence. The method determines the structure of the matrix $B$ and the control sequence $\mathbf{u}$ that drives the network from an initial to the desired state. 
 
To use these codes, you need to install (and add to the MATLAB path) the following toolboxes (codes):

1.) NOMAD solver.  The easiest way to install and to use the NOMAD solver is through the [OPTI toolbox](https://www.inverseproblem.co.nz/OPTI/index.php/DL/DownloadOPTI). 

2.) [MATLAB BGL]([https://github.com/dgleich/matlab-bgl](https://github.com/dgleich/matlab-bgl)) - this is only necessary if we want to generate the interconnection pattern of Duffing oscillator networks as geometric random graphs. 

3.) CONTEST toolbox- this is only necessary if we want to generate the interconnection pattern of Duffing oscillator networks as geometric random graphs. This toolbox was available [here](http://www.maths.strath.ac.uk/research/groups/numerical_analysis), but last time I checked (March 2020), the link was not working. The toolbox can also be downloaded from [here](https://github.com/jblocher/matlab-network-utilities/tree/master/contest). You only need a single MATLAB file "geo.m"

You only need a single MATLAB file: "geo.m" that generates random geometric graphs. 

## Folder structure and explanation

The following folders contain
1) "duffing_low_dimensonal" - codes for the low-dimensional Duffing oscillator network
2) "duffing_high_dimensonal" - codes for the high-dimensional Duffing oscillator network
3) "final_memory_network" - codes for the assosiative memory network

## Where to start from

1. In every folder there is a file "generate_dynamics.m". This file is used to generate a symbolic description of the function $\mathbf{f}(\mathbf{x})$ and the gradient of this function that is necessary for further computations. Open this file to generate the dynamics. This code generates two files that describe the symbolic dynamics and the symbolic gradients. You can use the same principle to generate the symbolic description of the dynamics of your network.
2. In the folder "duffing_low_dimensonal" and "duffing_high_dimensonal" there is a file "main_file_duffing_implicit_3.m". These are the main files for the Duffing oscillator networks. Similarly, in the folder "final_memory_network", there is a file "main_file_memory.m". This is the main file for the associative memory network. These files call other functions, and they are used to select the control nodes and to design control actions. Additional explanations are given in these main files.
