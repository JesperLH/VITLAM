function [clab,cticks] = plotBrainBloatedVoxels(map,mask,varargin)
%% Illustrate a 3D view of the brain with bloaded voxels
% Inputs:
%   map:    A vector of size V x 1. Containing the voxel loadings, spatial
%           map, to be visualized.
%   mask:   A binary 3D array of size (dimX x dimY x dimZ) with V ones.
%           Mapping the V loadings to the target "brain space" (fx. MNI)
%   varargin:
%      'threshold'      A vector with values in the interval between -1 and 1.
%                       Positve values are colored
%      'view'           Rotate 3D view, [v,w]
%      'zoom'           Magnification scale
%      'flipx'          Reverses the direction of the first dimension.
%      'save'           Saves the illustration at the specified location.
%                       Default is '' (empty string) which does not save
%                       the illustration.
%      'anatomical'     Change the underlying image to be different from
%                       the input mask (i.e. the functions second input).
%
%% Written by Jesper L. Hinrich
% Copyright (C) 2016 Technical University of Denmark - All Rights Reserved
% You may use, distribute and modify this code under the terms of the 
% Visualization Toolbox for Latent Variable Modeling of fMRI license.
% 
% You should have received a copy of the license with this file. 
% If not, please write to: jesper dot hinrich at gmail dot com, 
% or visit : https://brainconnectivity.compute.dtu.dk/ (under software)
%
% Get input parameters
paramNames = {'threshold'         ,'view','zoom','flipx','save','anatomical','color','colorbar'};
defaults   = {[-0.5,0.5],[130 60],1.3  ,false  , '',[],[],true};

[threshold, viewSetting, zoomLevel, flip_x,save_loc, anatom_mask,color, show_cbar]...
    = internal.stats.parseArgs(paramNames, defaults, varargin{:});
noc = size(map,2);
%% Error checking
if all(size(threshold) == [1,2]) || all(size(threshold') == [1,2])
    if min(threshold)<0 && max(threshold)>0
        %Its all good
    else
        error('The two threshold values must be on opposite sides of zero.')
    end
elseif all(squeeze(size(threshold)) == [1,1])
    %Its also all good
else
   error('Only one or two values supported for threshold.')
end

% is there enough space in the colormap
if noc+1>32 && min(map(:))<0
    error('Not possible to visualize more than 31 components at once.')
elseif noc+1>64 && min(map(:))>=0
    error('Not possible to visualize more than 63 components at once.')
elseif min(map(:))<0 && max(map(:))<0
    threshold = abs(threshold);
    map = abs(map);
else
    threshold = sort(threshold,'descend');
end

any_negatives = any(threshold < 0) && any(map(:)<0);
%% Create costume colormap
if isempty(color)
    if any_negatives %Both positive and negative map
        color = distinguishable_colors(noc*2,[1 1 1; 0 0 0; 0.5 0.5 0.5]);
    else %Only positive
        color = distinguishable_colors(noc,[1 1 1; 0 0 0; 0.5 0.5 0.5]);
    end
end

cmap= nan(64,3);
spacing = (64-floor(mod(64,noc)))/noc;
for i = 1:noc
    if any_negatives
        pspace = floor(spacing/2);
        nspace = ceil(spacing/2);
        cmap(((i-1)*spacing+1):((i*spacing)-nspace),:) = repmat(color(1+2*(i-1),:),pspace,1);
        cmap(((i-1)*spacing+1+pspace):((i*spacing)),:) = repmat(color(2*i,:),nspace,1);
    else
        cmap(((i-1)*spacing+1):(i*spacing),:) = repmat(color(i,:),spacing,1);
    end
end
cmap(end-ceil(spacing/2):end,:) = repmat([0,0,0],ceil(spacing/2)+1,1);

c_val = 1:spacing:64-1; c_val = c_val(1:noc);

%% Plot 3D brain mess
if isempty(anatom_mask)%Show the mask as an underlying image
    pa = patch(isosurface(permute(mask,[2 1 3])));
else %Show a user input (anatomical) mask as underlying image
    pa = patch(isosurface(permute(anatom_mask,[2 1 3])));
end

m = zeros(size(mask));
above_threshold = false(size(map));
for d = 1:noc
    voxels = false(size(mask));
    voxels(mask) = map(:,d)>threshold(1);
    m(voxels) = c_val(d);
    above_threshold(map(:,d)>threshold(1),d) = true;
    
    if any_negatives
        voxels = false(size(mask));
        voxels(mask) = map(:,d)<threshold(2);
        m(voxels) = c_val(d)+spacing-1;
        above_threshold(map(:,d)<threshold(2),d) = true;
    end
end

%Color voxels that are part of multiple maps black
voxels = false(size(mask));
voxels(mask) = sum(above_threshold,2)>1;
m(voxels) = 64; 

%Make non active voxels transparrent
m(m == 0) = NaN;

%If input convention is different from output convension
if flip_x
    m = m(end:-1:1,:,:); %Flip x-axis
end
assert(sum(~isnan(m(:)))>0,'There are no voxels active within the given threshold')

%% Show brain
set(pa,'edgecolor','none','facecolor',[0.5 0.5 0.5],'facealpha',0.2,'EdgeLighting','flat');
hold on;
axis off;
axis equal;
axis tight;
view(viewSetting);

%USE!!
PATCH_3Darray(m,'col',cmap,'clim',[1 64]); 
camzoom(zoomLevel)
hold off;

%% Determine ticks and labels for colormap
if any_negatives
    cticks = bsxfun(@plus,repmat(c_val,2,1),[spacing/4;3*spacing/4]);
    cticks = [cticks(:);64];
    clab = cell(length(cticks),1);
    clab(1:2:end-1) = strcat('Comp. ',strsplit(num2str(1:noc)),'(+)');
    clab(2:2:end-1) = strcat('Comp. ',strsplit(num2str(1:noc)),'(-)');
else
    cticks = [c_val+spacing/2,64];
    clab = cell(length(cticks),1);
    clab(1:end-1) = strcat('Comp. ',strsplit(num2str(1:noc)));
end
clab{end} = 'Overlapping';

if show_cbar
    cbar = colorbar;
    set(cbar,'Ticks',cticks,'TickLabels',clab,'Direction','reverse')    
end

%Save file
if ~isempty(save_loc)
    print(save_loc,'-dpng');
    print(save_loc,'-depsc');
end
end