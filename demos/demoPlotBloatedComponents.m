%% DEMO: Plotting Bloated Components
clear
load('./demos/example_data')

%% Plotting all components
plotBloatedComponents(A,mask);

%% Plotting a subset or using variable arguments
subset = 1:5;
plotBloatedComponents(A(:,subset),mask,...
            'inConvention','Radiological',... %Input convention (i.e. data)
            'outConvention','Neurological',...%Output convention
            'FontSize',18,... % Self-explanatory. Note titles has Fontsize+2 
            'save','',... %Saves .png and .eps to the given directory, if '' figures aren't saved
            'Position',[0 0 0.5 0.3],... %Figure position lower left and upper right cornor
            'threshold',[-5, 5]);

%% For more freedom, one can plot a single brain
Amaxabs = A/max(abs(A(:)));
figure('Units','normalized','position',[0.05,0.05,0.7,0.8])
for i = 1:9
    subplot(3,3,i)
    plotBrainBloatedVoxels(Amaxabs,mask,'threshold',[-i/10,i/10]);
    title(sprintf('A_{absmax}, threshold \\pm %1.1f',i/10))
end
