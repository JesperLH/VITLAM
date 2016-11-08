%% DEMO: plotNoisemap
clear
load('./demos/example_data')

%% Standard values
plotNoisemap(noisevar, mask, mask_affine_mat)

%% Standard values, but no affine mask => slice numbers are shown
plotNoisemap(noisevar, mask, [])

%% Using variable arguments
plotNoisemap(noisevar, mask, mask_affine_mat,...
            'nslices',12,...    %Number of slices to show (default: 9).
            'slices',6:3:40,... %Slices to illustrate
            'fontsize',18,... % Self-explanatory. Note titles has Fontsize+2 
            'flipx',true,...  %Reverses the direction of the first dimension.
            'save','',... %Saves .png and .eps to the given directory, if '' figures aren't saved
            'suptitles',{'Title mean plot', 'Title sd plot'},... %A string to be displayed as a super title
            'grid',[4 3],...        %Custom grid
            'Position',[0,0,0.3,0.5]); %Figure position