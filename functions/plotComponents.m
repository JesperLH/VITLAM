function plotComponents(A,S,mask,varargin)
% PLOTCOMPONENTS Illustrate 3D spatial maps, along with timeseries.
% INPUTS:
% "V" denotes number of voxels, "T" denotes number of time points, 
% "noc" denotes number of components, nSubs denotes number of subjects.
%
%   A:          A matrix of size V x noc. Containing the voxel loadings, 
%               spatial maps, to be visualized. 
%   S:          A matrix of size noc x T x nSubs. Containing temporal
%               activations.
%   mask:       A binary voxel space mask (size dimX x dimY x dimZ) with V
%               ones. Mapping the V loadings to the target "brain space"
%               (fx. MNI)
%   mask_mat    A 4 x 4 matrix (from SPM) with affine transformation from
%               subject space to MNI space.
%   varargin:   Passing optional input arguments
%      'inConvention'   
%      'outConvention'  
%      'TR'             Time resolution, i.e. seconds between timesteps
%                       (default: [])
%      'Fontsize'       Fontsize in figures (default: 18)
%      'linewidth'      
%      'save'           String with the save location of the illustrations,
%                       note figures are saved both as eps and png. Default
%                       is not saving.
%      'Position'       Figure position, units normalized
%      'threshold'      A range of threshold to apply
%      'Scaling'        
%      'color'          
%
%% Written by Jesper L. Hinrich and Sophia E. Bardenfleth
% Copyright (C) 2016 Technical University of Denmark - All Rights Reserved
% You may use, distribute and modify this code under the terms of the 
% Visualization Toolbox for Latent Variable Modeling of fMRI license.
% 
% You should have received a copy of the license with this file. 
% If not, please write to: jesper dot hinrich at gmail dot com, 
% or visit : https://brainconnectivity.compute.dtu.dk/ (under software)

%% Parse arguments and check if parameter/value pairs are valid 
paramNames = {'inConvention','outConvention','TR','Fontsize','LineWidth','save',...
              'Position','threshold','Scaling','color'};
defaults   = {'Radiological', 'Neurological',2.49,    18    , 2 ,'',...
              [0 0 800 425],[-.9:.1:-.1 , 0.1:0.1:.9],'absmax',...
              repmat([0 1 0; 1 0 0],1,1,size(A,2))};

[inputConvention, outputConvention, TR, f_size, l_size,save_loc,...
    fig_position,threshold,scaling,colors2use]...
    = internal.stats.parseArgs(paramNames, defaults, varargin{:});

%TODO customizable??
viewSetting = [50 60; 130 60; 270 0]; % ViewSetting affect zoomLevel
zoomLevel = [5/3,5/3,1.3]; % Zoom level affects display of s
s = {'','','Back Side'}; %Note, first and second argument plots weirdly

%% Validate custom input
%Radiological convention (patients left is on the right side)
%Neurological convention (patients left is on the left side)
assert(any(strcmpi(inputConvention,{'Radiological','Neurological'})),'Input convention unsupported')
assert(any(strcmpi(outputConvention,{'Radiological','Neurological'})),'Output convention unsupported')

% Scaling
if strcmpi(scaling,'absmax')
    A = A/max(abs(A(:)));
elseif strcmpi(scaling,'off')
    %Do nothing
else
    error('Unknown scaling method.')
    %Furture work could be different scaling techniques
end

if ~isempty(S)
    %% Plotting spatial maps and temporal activation
    Smean = mean(S,3);
    Ssd = sqrt(var(S,0,3));
    T=size(S,2);
    timepoints = (1:T)*TR;
end
plotArea = [1 4; 2 5; 3 6];

for j = 1:size(A,2) %Illustrate every component
    
    figure('Position',fig_position)
    %% Plotting S views
    for i = 1:3
        if ~isempty(S)
            subplot(3,3,plotArea(i,:))
        else
            subplot(2,3,plotArea(i,:))
        end
        %If input convention is different from output convension
        flip_x = ~strcmpi(inputConvention,outputConvention);
        plotBrain(A(:,j),mask,'flipx',flip_x,'view',viewSetting(i,:),...
                                    'zoom',zoomLevel(i),...
                                    'threshold',threshold,...
                                    'color',colors2use(:,:,j));
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

    if ~isempty(S)
        %% Plot temporal activation S
        subplot(3,3,7:9)
        hold on
        plot(timepoints,Smean(j,:),'Color','b','LineWidth',l_size);
        if size(S,3) > 1
            plot(timepoints,Smean(j,:)+2*Ssd(j,:),':','Color','r','LineWidth',l_size-.5);
            plot(timepoints,Smean(j,:)-2*Ssd(j,:),':','Color','r','LineWidth',l_size-.5);
        end
        xlabel('time (s)')
        ylabel('Activation')

        set(gca,'fontsize',f_size)
        if nargin <= 7
            str = sprintf('Latent Component %i',j);
            title(str,'FontSize',f_size+2);
        end
    
        %Move slightly to the left
        pos = get(gca, 'Position');
        xoffset = -0.04;
        pos(1) = pos(1) + xoffset;
        pos(3) = pos(3) -2* xoffset;
        set(gca, 'Position', pos)
    end
    %Save file
    if ~isempty(save_loc)
        results_name = sprintf('_d%i',j);
        print(sprintf('%s%s',save_loc,results_name),'-dpng');
        print(sprintf('%s%s',save_loc,results_name),'-depsc')
    end
    
end
end