%% Demostration: Plot component spatial maps and associated temporal activation
clear; close all; clc
load('./demos/example_data')
% 'A' has dimension V x q, where each column is a spatial map
% 'S' has dimension q x T x subjects,
% 'mask' is a 3D binary indicator matrix, of which voxels have been used in
%        the analysis.

%% Plotting all components
plotComponents(A,S,mask)

%% Plotting a subset or using variable arguments
subset = 1:2;
plotComponents(A(:,subset),S(subset,:,:),mask,...,)
            'inConvention','Radiological',... %Input convention (i.e. data)
            'outConvention','Neurological',...%Output convention
            'TR',2.49,... %Time resolution in secounds
            'FontSize',18,... % Self-explanatory. Note titles has Fontsize+2 
            'LineWidth',2,... % Width of the temporal activation
            'save','',... %Saves .png and .eps to the given directory, if '' figures aren't saved
            'Position',[0 0 800 450],... %Figure position lower left and upper right cornor
            'threshold',[-0.9:0.1:-0.3 , 0.3:0.1:0.9],...
            'Scaling','absmax');
%% Plotting only spatial maps
plotComponents(A(:,1),[],mask,'Position',[50 50 800 350])

%% Plotting only positive maps with custome coloring
rng('shuffle')
custom_color = rand(1,3,2);
plotComponents(A(:,2:3),S(2:3,:,:),mask,'threshold',0.1:0.1:0.9,'color',custom_color)

%% Plotting  maps with custome coloring
rng('shuffle')
custom_color = rand(2,3,2);
plotComponents(A(:,2:3),S(2:3,:,:),mask,'color',custom_color)