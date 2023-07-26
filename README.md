# Salt_Inclusion_in_Rocks

Salt is a material characterized by its fluidity, and we model it as an incompressible, nonlinear viscous/viscous creep material. The salt is surrounded by shale/sandstone, and the entire formation is subject to geostress. The purpose of the matlab and julia code is to solve for the stress distribution around the salt. 

Since the results were not published in any papers, some documents have not been published, and the specific equations should refer to the unpublished documents.

# The folder 'literature'
This folder includes references related to natural hydrocarbon seepage, including observation of natural hrdrocarbon seepage, well logging data and seismic imaging, porosity wave model, creep behavior of salt, etc. 

# The folder 'SeepageVideo'
The files in this folder show the natural hydrocarbon seepage process, these videos are recorded by cameras.
Video Source:
Wang, Binbin, Scott A. Socolofsky, John A. Breier, and Jeffrey S. Seewald. "Observations of bubbles in natural seep flares at MC 118 and GC 600 using in situ quantitative imaging." Journal of Geophysical Research: Oceans 121, no. 4 (2016): 2203-2230.

# The folder 'matlab'
Matlab code can be debugged and plotted more conveniently. The codes 'Stokes2D.m' and 'Stokes3D.m' are modified on the basis of FASTICE to simulate the stress distribution of one material under the in-situ geostress. The suffix 2D indicates that the simulation region is two-dimensional, and the suffix 3D indicates that the simulation region is three-dimensional.

The role of the functions 'X_mean.m', 'Z_mean.m', 'XZ_mean.m', 'XY_mean3D.m', 'XZ_mean3D.m', 'YZ_mean3D.m', 'XYZ_mean3D.m' are the same as that of the package 'ParallelStencil.jl' in the julia code. The main functions are 'SaltInclusion2D.m' and 'SaltInclusion3D.m', and these matlab codes are not parallelized, so the computational time is relatively long.

The matlab code can be calculated by the supercomputer Sherlock, 'runcode.sbatch' is a script that uses cpu to run the code. 

'createGUI.m' is a code to create a parameter input interface, which can visualize the input parameters to an input interface.
![](interface.png)


# The folder 'julia'
The Julia code has been parallelized based on the matlab codes. The parallel julia code only considers salt as a viscous material, and does not use the MD model to calculate the change in strain caused by the creep. The change in strain caused by the creep can be calculated by modifying the function in the pseudo-transient loop. When using grid computing, Julia code uses the package code 'ParallelStencil.jl'. You can read the code of the package via this link:
https://github.com/omlins/ParallelStencil.jl

Julia code runs through Jupyter Notebook, which requires the installation of Julia in the environment of Jupyter Notebook, and install related packages accordingly. The julia version I am using is 1.8.5. To install a package in Julia, the command is Pkg.add("PackageName"). To use this command, you need to first import the Pkg package. You can do this by following the steps below:

1. open the Julia command line.
   Enter the following command to import Pkg:
   using Pkg.
2. Use the following command to install the required package. For example, if you want to install a package called "Example", you would enter:
   Pkg.add("Example"). 
   Replace "Example" with the name of the package you want to install.

If you want to install a specific version of a package, you can specify it in Pkg.add(). For example, if you wanted to install version 2.1.3 of the Example package, you would enter: Pkg.add(PackageSpec(name="Example", version="2.1.3")). 




