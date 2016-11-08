function plotNoisemap(noise,mask,mask_affine,varargin)
% NOISEMAP Illustrate mean and standard deviation over the columns of noise.
% 
% INPUTS:
% "V" denotes number of voxels, "T" denotes number of time points, 
% "noc" denotes number of components, nSubs denotes number of subjects.
%
%   noise:      A matrix of size V x nSubs. Containing the voxel variance
%               for each subjeect.
%   mask:       A binary voxel space mask (size dimX x dimY x dimZ) with V
%               ones. Mapping the V loadings to the target "brain space"
%               (fx. MNI)
%   mask_affine A 4 x 4 matrix (from SPM) with affine transformation from
%               subject space to MNI space.
%   varargin:   Passing optional input arguments
%      'nslices'        Number of slices to show (default: 9).
%      'slices'         Slices to illustrate, default is 9 equi-distance
%                       slices.
%      'fontsize'       Fontsize in figures (default: 16).
%      'flipx'          Reverses the direction of the first dimension.
%      'save'           String with the save location of the illustrations,
%                       note figures are saved both as eps and png. Default
%                       is not saving.
%      'anatomical'     Show underlying slices from another mask, fx an
%                       anatomical (NOT IMPLEMENTET)
%      'suptitles'       A cell array (size 2), with two strings to be 
%                       displayed as a super titles.
%                       (default:[], which is no super titles)
%      'grid'           Defines a custom grid for the spatial slices,
%                       with a vector with the number of [horizontal,
%                       vertical] slices.
%      'Position'       Figure position, units normalized
%% Written by Jesper L. Hinrich and Sophia E. Bardenfleth
% Copyright (C) 2016 Technical University of Denmark - All Rights Reserved
% You may use, distribute and modify this code under the terms of the 
% Visualization Toolbox for Latent Variable Modeling of fMRI license.
% 
% You should have received a copy of the license with this file. 
% If not, please write to: jesper dot hinrich at gmail dot com, 
% or visit : https://brainconnectivity.compute.dtu.dk/ (under software)

% Get input parameters
paramNames = {'nslices','slices','fontsize','flipx','save','anatomical',...
              'suptitles','grid','Position'};
defaults   = {9,[],16,false, '',mask,[],[],[0 0 0.4, 0.6]};

[nslices, slices, f_size, flip_x,save_loc,anatom_mask,super_titles...
    ,custom_grid,fig_position]=internal.stats.parseArgs(paramNames,...
                                                    defaults, varargin{:});

if isempty(slices)
    slices = floor(linspace(1,size(mask,3),nslices));    
end
%Get NMI slice z-coordinates
if ~isempty(mask_affine)
    slices_MNI = round(slices*abs(mask_affine(1,1)) - abs(mask_affine(3,4)));
end

%% Determine plots required for slices
nslices = length(slices);
if isempty(custom_grid)
    vplots = floor(sqrt(nslices)); %Num. Vertical plots
    hplots = ceil(sqrt(nslices)); %Num. Horizontal plots
    if hplots*vplots < length(slices) %Fixes an off by one situation
        vplots = vplots+1;
    end
else
    hplots = custom_grid(1);
    vplots = custom_grid(2);
    assert(hplots*vplots>=length(slices),['Grid size insufficient!',...
        ' Must allow for atleast %i subplots. Input had %i*%i = %i subplots'],...
        length(slices),hplots,vplots,hplots*vplots)
end

%% Illustrate mean and std noise maps
for i = 1:2 %i=1, avg noise map. i=2 std noise map
    
    %Create brainmap
    noisebrain = nan(size(mask));
    if i == 1
        noisebrain(mask) = mean(noise,2);
        %caxis_limits = [0,2];
        caxis_limits = [min(noisebrain(:)),max(noisebrain(:))];
    elseif i == 2
        noisebrain(mask) = sqrt(var(noise,[],2));
        caxis_limits = [min(noisebrain(:)),max(noisebrain(:))];
    end
    if caxis_limits(1) < 0
        caxis_limits(1) = caxis_limits(1)*1.1;
    else
        caxis_limits(1) = caxis_limits(1)*.9;
    end
    if flip_x %Neurological convention.. Reverse X
        noisebrain = noisebrain(end:-1:1,:,:);
    end
    %scrsz = get(groot,'ScreenSize');
    %figure('Position',[1 scrsz(4)/2 scrsz(3)*2/3 scrsz(4)*2/3])
    figure('Units','Normalized','Position',fig_position)
    for s=1:length(slices)
        subplot(vplots,hplots,s)
        pcolor(squeeze(noisebrain(:,:,slices(s)))');
        shading flat;
        %Make visualization nice
        caxis(caxis_limits);
        axis image;
        set(gca,'XTick',[],...
                'XTickLabel',[],...
                'YTick',[],...
                'YTickLabel',[],...
                'Box','off')
        %Slightly bigger slices (reduces white space between subplots)
        sub_pos = get(gca,'position'); % get subplot axis position
        offset = [-0.075,0,0,0];
        set(gca,'position',sub_pos.*[1 1 1.1 1.1]+offset) % stretch its width and height
        
        if ~isempty(mask_affine)
            text(0,1.05,sprintf('z=%i',slices_MNI(s)),'FontSize',f_size,...
                'Units','Normalized');
        else
            text(0,1.05,sprintf('slice=%i',slices(s)),'FontSize',f_size,...
                'Units','Normalized');
        end
    end

    set(gcf,'defaultaxesfontsize',f_size+2)
    if ~isempty(super_titles)
        if ~isempty(super_titles{i})
            suptitle(super_titles{i});
        end
    elseif i == 1
        suptitle('Mean variance across subject (mean noise)'); 
    elseif i == 2
        suptitle('Standard deviation across subject (s.d. of noise)')
    end
    %
    %colormap jet
    B=colorbar; 
    set(B, 'Position', [.875 .15 .0381 .6150]) 
    set(B,'FontSize',f_size)

    set(gcf,'PaperPositionMode','auto')
    %tightfig()
    if ~isempty(save_loc)
         if i == 1
             filename = [save_loc '_mean'];
         elseif i==2
             filename = [save_loc '_sd'];
         end
         print(filename,'-depsc')
         print(filename,'-dpng')
    end
end
end
