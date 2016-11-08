function plotBloatedComponents(A,mask,varargin)
% PLOTBLOATEDCOMPONENTS Illustrate the components in A in 3D brain views,
%                       where voxels above a given threshold is shown as
%                       bloated cubes. 
% INPUTS:
% "V" denotes number of voxels, "noc" denotes number of components.
%
%   A:          A matrix of size V x noc. Containing the voxel loadings, 
%               spatial maps, to be visualized. 
%   mask:       A binary voxel space mask (size dimX x dimY x dimZ) with V
%               ones. Mapping the V loadings to the target "brain space"
%               (fx. MNI)
%   varargin:   Passing optional input arguments
%      'inConvention'   
%      'outConvention'  
%      'Fontsize'       Fontsize in figures (default: 18)
%      'linewidth'      
%      'save'           String with the save location of the illustrations,
%                       note figures are saved both as eps and png. Default
%                       is not saving.
%      'Position'       Figure position, units normalized
%      'threshold'      
%      'Scaling'        
%      'color'          
%% Written by Jesper L. Hinrich
% Copyright (C) 2016 Technical University of Denmark - All Rights Reserved
% You may use, distribute and modify this code under the terms of the 
% Visualization Toolbox for Latent Variable Modeling of fMRI license.
% 
% You should have received a copy of the license with this file. 
% If not, please write to: jesper dot hinrich at gmail dot com, 
% or visit : https://brainconnectivity.compute.dtu.dk/ (under software)
% Get input parameters

%% Parse arguments and check if parameter/value pairs are valid
paramNames = {'inConvention','outConvention','FontSize','save',...
    'Position','threshold','color'};
defaults   = {'Radiological', 'Neurological',    18 ,'',...
    [0 0 0.5 0.3],[-0.5,0.5],[]};

[inputConvention, outputConvention, f_size, save_loc,...
    fig_position,threshold,colorlist]...
    = internal.stats.parseArgs(paramNames, defaults, varargin{:});

%TODO customizable??
viewSetting = [140 60; 210 60; 0 0]; % ViewSetting affect zoomLevel
zoomLevel = [2,2,1.3]; % Zoom level affects display of s
s = {'','','Back Side'}; %Note, first and second argument plots weirdly

%% Validate custom input
%Radiological convention (patients left is on the right side)
%Neurological convention (patients left is on the left side)
assert(any(strcmpi(inputConvention,{'Radiological','Neurological'})),'Input convention unsupported')
assert(any(strcmpi(outputConvention,{'Radiological','Neurological'})),'Output convention unsupported')

plotArea = [1 4; 2 5; 3 6];

figure('Units','normalized','Position',fig_position)
%% Plotting S views
for i = 1:3
    subplot(3,3,plotArea(i,:))
    %If input convention is different from output convension
    flip_x = ~strcmpi(inputConvention,outputConvention);
    [clab,cticks] = plotBrainBloatedVoxels(A,mask,'flipx',flip_x,'view',viewSetting(i,:),...
        'zoom',zoomLevel(i),'threshold',threshold,'color',colorlist,...
        'colorbar',false);
    
    set(gca,'Fontsize',f_size)
    title(s{i},'FontSize',f_size+2)
    
    % Offset to the left
    pos = get(gca, 'Position');
    xoffset = -0.06;
    if i ==1
        pos(1) = pos(1) + xoffset;
        pos(3) = pos(3) - 0.5*xoffset;
    elseif i==2
        pos(1) = pos(1) +0.5*xoffset;
        pos(3) = pos(3) - 0.5*xoffset;
    elseif i==3
        pos(3) = pos(3) -.5* xoffset;
    end
    set(gca, 'Position', pos)
end

cbar = colorbar('location','southoutside');
cticks(1:end-1) = cticks(1:end-1)-diff(cticks(1:2))*0.4;
set(cbar,'Ticks',cticks,'TickLabels',clab,'Position',[.05,.2,.9,.05],'fontsize',10)
%Save file
if ~isempty(save_loc)
    results_name = '_bloated';
    print(sprintf('%s%s',save_loc,results_name),'-dpng');
    print(sprintf('%s%s',save_loc,results_name),'-depsc')
end

end