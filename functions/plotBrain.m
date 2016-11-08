function plotBrain(map,mask,varargin)
%% Illustrate a 3D view of the brain
% Inputs:
%   map:    A vector of size V x 1. Containing the voxel loadings, spatial
%           map, to be visualized. 
%   mask:   A binary 3D array of size (dimX x dimY x dimZ) with V ones.
%           Mapping the V loadings to the target "brain space" (fx. MNI)
%   varargin:
%      'threshold'      A vector with values in the interval between -1 and 1. 
%                       Positve values are colored green, negative values
%                       red.
%      'view'           Rotate 3D view, [v,w]
%      'zoom'           Magnification scale
%      'flipx'          Reverses the direction of the first dimension.
%      'save'           Saves the illustration at the specified location. 
%                       Default is '' (empty string) which does not save
%                       the illustration.
%      'anatomical'     Change the underlying image to be different from
%                       the input mask (i.e. the functions second input).
%      'color'          Determines the color used, default is [0 1 0; 1 0 0]
%                       which is green and red. (for positive and negative
%                       values)
%
%% Written by Jesper L. Hinrich
% Copyright (C) 2016 Technical University of Denmark - All Rights Reserved
% You may use, distribute and modify this code under the terms of the 
% Visualization Toolbox for Latent Variable Modeling of fMRI license.
% 
% You should have received a copy of the license with this file. 
% If not, please write to: jesper dot hinrich at gmail dot com, 
% or visit : https://brainconnectivity.compute.dtu.dk/ (under software)
% Get input parameters
paramNames = {'threshold'         ,'view','zoom','flipx','save','anatomical','color'};
defaults   = {[-.9:.1:.1,.1:.1:.9],[50 60],1.3  ,false  , '',[],[0 1 0; 1 0 0]};

[threshold, viewSetting, zoomLevel, flip_x,save_loc, anatom_mask, vcolors]...
    = internal.stats.parseArgs(paramNames, defaults, varargin{:});

%% Plot 3D brain mess
%mask_dim = size(mask);
if isempty(anatom_mask)
    pa = patch(isosurface(mask)); %Show the mask as an underlying image
else
    %Show a user input (anatomical) mask as underlying image
    pa = patch(isosurface(anatom_mask));
end
set(pa,'edgecolor','none','facecolor',[0.5 0.5 0.5],'facealpha',0.2);
hold on;
axis off;
axis equal;
axis tight;
view(viewSetting);
for t = 1:length(threshold)
    m = mask;
    if threshold(t) > 0
        dispColor = vcolors(1,:);% Default 'Green'
        m(m==1)=map>threshold(t);
    else
        dispColor = vcolors(2,:); %Default 'Red'
        m(m==1)=map<threshold(t);
    end
    %If input convention is different from output convension
    if flip_x
        m = m(end:-1:1,:,:); %Flip x-axis
    end
    pb = patch(isosurface(m)); %Add colored brain
    set(pb,'edgecolor','none','facecolor',dispColor...
        ,'facealpha',abs(threshold(t)));
end

camzoom(zoomLevel);
hold off;

%Save file
if ~isempty(save_loc)
    print(save_loc,'-dpng');
    print(save_loc,'-depsc');
end
end

