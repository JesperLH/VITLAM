# VITLAM - Visualization Toolbox for Latent Variable Modeling of fMRI
The Visualization Toolbox for Latent Variable Modeling of fMRI holds a 
collection of visualization methods for latent variable analysis in 
functional magnetic resonance imaging (fMRI). The methods are implemented 
in Matlab™. All code can be used freely in research and other non-profit 
applications. If you publish results obtained with this toolbox we kindly 
ask that our and other relevant sources are properly cited. 

This toolbox has been developed at:

The Technical University of Denmark, 
Department for Applied Mathematics and Computer Science,
Section for Cognitive Systems.

The toolbox was developed in connection with the Brain Connectivity project 
at DTU (https://brainconnectivity.compute.dtu.dk/) .

## Visualizations:

* plotBrain:              
	- Plots 3D view of the brain with the spatial activation of a specified latent component (with intensity and sign).
* plotComponents:         
	- Plots three 3D views of the brain (and the associated temporal activation) with the spatial activation of a specified latent component (positive and/or negative).
* plotBrainBloatedVoxels: 
	- Plots 3D view with bloated indicators for spatial active voxels, can compare multiple latent components and assess their overlap.
* plotBloatedComponents:  
	- Similar to "plotComponent" but with bloated voxels instead.
* plotSpatialSlices:      
	- Plots the spatial activation of a components, as slices of the brain. It is also possible to show the temporal activation (or its power spectrum) and/or histogram of the component.
* plotNoisemap:           
	- Plots heteroscedastic voxel noise as spatial slices (average and standard deviation across subjects)

Common properties

* Highly customizable output, see demostration scripts for further detail.

## Demonstrators:
How to call the functions and illustrating the effect of different parameter settings.

* demoPlotComponents
* demoPlotBloatedComponents
* demoPlotSpatialSlices
* demoPlotNoisemap

## References
For most use cases a link to either this toolbox or the Brain Connectivity project will surfice. The foundations for **plotBrain, plotComponents, plotBloatedVoxels** and **plotBloatedComponents** were introduced in [1]. The foundations for **plotSpatialSlices** is presented in an article currently under review (will be added in the future).

 * [1] Hinrich, J. L., Bardenfleth, S. E., Røge, R. E., Churchill, N. W., Madsen, K. H., & Mørup, M. (2016). Archetypal Analysis for Modeling Multisubject fMRI Data. IEEE Journal of Selected Topics in Signal Processing, 10(7), 1160-1171.
