%% DEMO: plotSpatialSlices
clear
load('./demos/example_data')

%% Plotting all components
% Shows timeseries and histogram over spatial loadings
plotSpatialSlices(A,S, mask, mask_affine_mat)

%% 
set = [1];%[1,2];
% No timeseries
plotSpatialSlices(A(:,set),[], mask, mask_affine_mat)
% Show power spectrum density
plotSpatialSlices(A(:,set),S(set,:,:), mask, mask_affine_mat,'psd',true)
%% Timeseries below
plotSpatialSlices(A(:,set),S(set,:,:), mask, mask_affine_mat,...
    'placetime','below')
%% Show power spectrum density and histogram instead of last slice.
plotSpatialSlices(A(:,set),S(set,:,:), mask, mask_affine_mat,...
    'placetime','below','psd',true,'histinslices',true)

%% Plotting using variable arguments
% subset = 1:2;
plotSpatialSlices(A(:,set),S(set,:,:),mask,mask_affine_mat,...
            'nslices',12,...    %Number of slices to show (default: 9).
            'slices',6:3:40,... %Slices to illustrate
            'threshold_z',1.5,...%Acceptance threshold for spatial activation
            'fontsize',18,... % Self-explanatory. Note titles has Fontsize+2 
            'save','',... %Saves .png and .eps to the given directory, if '' figures aren't saved
            'suptitle','Super title here',... %A string to be displayed as a super title
            'placetime','right',... %Place time to the 'right' or 'below' slices
            'ratio',1,...           %Display ratio (space) between slices and timeseries
            'TR',2.49,...           %Time resolution in secounds
            'grid',[4 3],...        %Custom grid
            'histogram',true,...    %Show histogram
            'bins',100,...          %Number of histogram bins 
            'histtype','pdf',...    %Type of histogram normalization
            'histinslices',false,...%Show the histogram instead of the last slice
            'psd',true); %Show Power Spectrum Density instead of timeseries.
        
%% Variable input arguments can also be passed as a structure
% Note field names must match the option stated above
opts_struct.nslices = 25;
opts.histinslices = true;
plotSpatialSlices(A(:,set),S(set,:,:),mask,mask_affine_mat,'opts',opts_struct)