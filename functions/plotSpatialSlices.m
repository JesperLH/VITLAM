function plotSpatialSlices(MAPS, TIMESERIES, mask, mask_mat,varargin)
% PLOTSPATIALSLICES Illustrate active spatial maps, along with timeseries
%                   or power density spectrum, and with an option to show 
%                   component histograms.
% INPUTS:
% "V" denotes number of voxels, "T" denotes number of time points, 
% "noc" denotes number of components, "nSubs" denotes number of subjects.
%
%   MAPS:       A matrix of size V x noc. Containing the voxel loadings, 
%               spatial maps, to be visualized. 
%   TIMESERIES: A matrix of size noc x T x nSubs. Containing temporal
%               activations.
%   mask:       A binary voxel space mask (size dimX x dimY x dimZ) with V
%               ones. Mapping the V loadings to the target "brain space"
%               (fx. MNI)
%   mask_mat    A 4 x 4 matrix (from SPM) with affine transformation from
%               subject space to MNI space.
%   varargin:   Passing optional input arguments
%      'opts'           A structure with field names matching the variable
%                       input argument name. Note if an option is specified
%                       both in "opts" and in varargin, the value in "opts"
%                       will be used.
%      'nslices'        Number of slices to show (default: 9).
%      'slices'         Slices to illustrate, default is 9 equi-distance
%                       slices.
%      'threshold_z'    Acceptance threshold for spatial activation
%                       (default: 1)
%      'fontsize'       Fontsize in figures (default: 12).
%      'save'           String with the save location of the illustrations,
%                       note figures are saved both as eps and png. Default
%                       is not saving.
%      'suptitle'       A string to be displayed as a super title 
%                       (default:[], which is no super title)
%      'placetime'      Place time to the 'right' or 'below' slices 
%                       (default: 'right').
%      'ratio'          Display ratio (space) between slices and timeseries
%                       plot, ex. 0.5 => timeseries fill half as much as
%                       slices, 1 => equal space and 2 => timeseries fills
%                       twice as much as slices. (default: 0.75)
%      'TR'             Time resolution, i.e. seconds between timesteps
%                       (default: [])
%      'grid'           Defines a custom grid for the spatial slices,
%                       with a vector with the number of [horizontal,
%                       vertical] slices.
%      'histogram'      Show histogram (default: false)
%      'bins'           Number of histogram bins (default: 50)
%      'histtype'       Type of histogram normalization(default:'proportion')
%         'frequency'   Show counts in each bin
%         'proportion'  Shows proportion, i.e. counts/(total count)
%         'area'        Shows area normalized proportion, i.e. counts/area.
%                       Trapezoidal numerical integration is used to
%                       approximate the histogram integral(area under the
%                       curve).
%         'fdf'         frequency density estimate, i.e. counts/binwidth .
%         'pdf'         probability density estimate, 
%                       i.e. frequency density/(total counts)
%                       ''
%                       ''
%      'histinslices'   Show the histogram instead of the last slice
%                       (default value: false).
%      'psd'            Show Power Spectrum Density instead of timeseries.
%                       The psd is calculated using "pwelch" for each
%                       subject (default value: false).
%      'time_text'      Cell-array with strings shown above timeseries (if
%                       plotted) for each component
%      'text_below_slices'  Cell-array strings shown below maps for each 
%                           component
%
%
%% Written by Jesper L. Hinrich and S�ren F. V. Nielsen
% Copyright (C) 2016 Technical University of Denmark - All Rights Reserved
% You may use, distribute and modify this code under the terms of the 
% Visualization Toolbox for Latent Variable Modeling of fMRI license.
% 
% You should have received a copy of the license with this file. 
% If not, please write to: jesper dot hinrich at gmail dot com, 
% or visit : https://brainconnectivity.compute.dtu.dk/ (under software)

%% Get parameters
% Parse arguments and check if parameter/value pairs are valid 
paramNames = {'nslices','slices','threshold_z','fontsize','save',...
              'suptitle','placetime','ratio','TR','grid','histogram',...
              'bins','histtype','histinslices','psd','opts', 'time_text',...
              'text_below_slices'};
defaults = {9, [], 1 , 12 , [],...
            [],'right', 0.75, [], [], false,...
            50, 'proportion', false, false,[],[],[]};

[numslices, slices, thresh, fs, save_loc,...
    super_title, place_time, slice_time_ratio, TR, custom_grid, show_hist,...
    nbins,hist_method,histogramInSlice,show_psd, opts, time_text, ...
    text_below_slices]...
    = internal.stats.parseArgs(paramNames, defaults, varargin{:});
%Supporting passing a struct "opts" with optional parameters. In order to
%allow simultaneous an "opts" structure and variable input arguments(varargin) 
% the default values need to be whatever they became from the above line.
% Note that if an option is specified both in "opts" and in varargin, the
% value in "opts" will be used.
if ~isempty(opts)
    numslices = mgetopt(opts,'nslices',numslices);
    slices = mgetopt(opts,'slices',slices);
    thresh = mgetopt(opts, 'threshold_z', thresh);
    fs = mgetopt(opts,'fontsize',fs);
    save_loc = mgetopt(opts,'save',save_loc);
    super_title = mgetopt(opts,'suptitle',super_title);
    place_time = mgetopt(opts,'placetime',place_time);
    time_text = mgetopt(opts,'time_text',time_text); % text above time
    text_below_slices = mgetopt(opts,'text_below_slices',text_below_slices); 
    slice_time_ratio = mgetopt(opts,'ratio',slice_time_ratio);
    TR = mgetopt(opts,'TR',TR);
    custom_grid = mgetopt(opts,'grid',custom_grid);
    show_hist = mgetopt(opts,'histogram',show_hist);
    nbins = mgetopt(opts,'bins',nbins);
    hist_method = mgetopt(opts,'histtype',hist_method);
    histogramInSlice = mgetopt(opts,'histinslices',histogramInSlice);
    show_psd = mgetopt(opts,'psd',show_psd); %Power Spectrum Density
end

% Make sure slices are initialized
if isempty(slices)
    slices = floor(linspace(1,size(mask,3),numslices));
end
numslices = length(slices);

%% Figure out grid dimensions
%Determine plots required for slices
if isempty(custom_grid)
    vplots = floor(sqrt(numslices)); %Num. Vertical plots
    hplots = ceil(sqrt(numslices)); %Num. Horizontal plots
    if hplots*vplots < length(slices) %Fixes an off by one situation
        vplots = vplots+1;
    end
    hslices = ceil(sqrt(numslices));
else
    hplots = custom_grid(1);
    vplots = custom_grid(2);
    assert(hplots*vplots>=length(slices),['Grid size insufficient!',...
        ' Must allow for atleast %i subplots. Input had %i*%i = %i subplots'],...
        length(slices),hplots,vplots,hplots*vplots)
    hslices = hplots;
end
%Determine additional plots needed for time series
if strcmp(place_time,'below')
    vplots = vplots+~isempty(TIMESERIES);
elseif strcmp(place_time,'right')
    hplots = ceil(hplots*(1+slice_time_ratio*(~isempty(TIMESERIES) | show_hist)));
end

%% 
%Determine "statistical" significans, by z-scoring maps and thresholding
Az = zscore(MAPS')';
Az(abs(Az)<thresh) = NaN;

if ~isempty(mask_mat) % define slices in z-coordinates
    slices_MNI = round(slices*abs(mask_mat(1,1)) - abs(mask_mat(3,4)));
end

%% Illustrating the components
noc = size(MAPS,2); % number of components
for nc = 1:noc
    if strcmp(place_time,'below') %? A gull ?? TODO: Smarter figure init.!
        h=figure('Position',[0,0,1000,1000]);
    elseif strcmp(place_time,'right') && isempty(custom_grid)
        h=figure('Position',[0,0,500*(2+slice_time_ratio),600]);
    elseif ~isempty(custom_grid)
        h=figure('Position',[0,0,500*(2+slice_time_ratio),450]);
    end
    
    brain = nan(size(mask));
    brain(mask) = Az(:,nc);
    %brain = brain(end:-1:1,:,:); %Neurological convention.. Reverse X
    for s = 1:length(slices);
        % Determine subplot location
        if strcmp(place_time,'below')
            subplot(vplots,hplots,s)
        elseif strcmp(place_time,'right')
            subplot(vplots,hplots,1+mod(s-1,hslices)+...
                (ceil(s/hslices)-1)*hplots)
        end
        
        IMG = convertSliceToImage(brain(:,:,slices(s)), ...
            mask(:,:,slices(s)) );
        image( IMG )
        axis image
        hold on
        if ~isempty(mask_mat)
            text(0,0,sprintf('z=%i',slices_MNI(s)),'FontSize',fs);
        else
            text(0,0,sprintf('slice=%i',slices(s)),'FontSize',fs);
        end
        hold off
        set(gca,'XTick',[],...
            'XTickLabel',[],...
            'YTick',[],...
            'YTickLabel',[],...
            'Box','off')
        %Slightly bigger slices (reduces white space between subplots)
        % get subplot axis position, then stretch its width and height.
        sub_pos = get(gca,'position'); 
        set(gca,'position',sub_pos.*[1 1 1.15 1.15])
        
        % Show histogram over spatial values
        if show_hist && s == length(slices) && histogramInSlice
            illustrateHistogram(MAPS(:,nc),nbins,hist_method,fs);
            sub_pos = get(gca,'position'); % get subplot axis position
            set(gca,'position',sub_pos.*[1 1 0.8 1]+[sub_pos(3)*0.15,0,0,0])
        end
        % get lower left corner position and lower right corner position
        if s == length(slices)
            pos = get(gca,'Position');
            pos = [pos(1)-pos(3)*(hplots-1)+1/2*pos(3) pos(2)-1/10 (hplots-1)*pos(3) 1/10];
            text_ax = axes('Position',pos,'Visible','off');
        end
        
    end
    
    % Draw text
    axes(text_ax);
    text(0.5,0.5,text_below_slices{nc}, 'FontSize',fs)
    
    if ~isempty(TIMESERIES)% plot time series
        %% Determine subplot location
        if strcmp(place_time,'below')
            subplot(vplots+1,hplots,((vplots*hplots) +1):((vplots+1)*hplots))
        elseif strcmp(place_time,'right')
            time_idx = bsxfun(@plus,(hslices+1):(hplots),...
                    bsxfun(@times,ones(vplots,ceil(hslices*slice_time_ratio))...
                    ,(0:(vplots-1))')*hplots);
            if ~histogramInSlice %Plot histogram below timeseries
                height = size(time_idx,1);
                if height > 1
                    %Illustrate histrogram
                    hist_start = double(height==2)+floor(height/2)+mod(height,2);
                    t_idx = time_idx(hist_start:end,:);
                    subplot(vplots,hplots,t_idx(:))
                    illustrateHistogram(MAPS(:,nc),nbins,hist_method,fs);
                    set(gca,'YAxisLocation','right')
                    %Switch to timeseries area..
                    t_idx = time_idx(1:(hist_start-1),:);
                    subplot(vplots,hplots,t_idx(:))
                else
                    warning(['Histogram not shown! Too few vertical '...
                             'subplots/slices (minimum is 2)'])
                end
            else 
                subplot(vplots,hplots,time_idx(:))
            end
        end
        
        if ~show_psd 
            %% Plot timeseries
            if isempty(TR)
                xtime = 1:size(TIMESERIES,2);
            else
                xtime = (1:size(TIMESERIES,2))*TR;
            end
            
            if size(TIMESERIES,3) == 1;
                plot(xtime,TIMESERIES(nc,:),'k','Linewidth',2)
            else
                %Plot individual timeseries for each subject/session
                for s = 1:size(TIMESERIES,3)
                    plot(xtime,TIMESERIES(nc,:,s),'g','Linewidth',.75); hold on
                end
                %Plot mean
                plot(xtime,mean(TIMESERIES(nc,:,:),3),'k','Linewidth',2)
                
            end
            yax = max(abs(min(min(TIMESERIES(nc,:,:)))),max(max(TIMESERIES(nc,:,:))));
            axis([1,xtime(end),-yax*1.001,yax*1.001])
            
            set(gca,'Fontsize',fs)
            xlabel('Time'); ylabel('Activation')
            if isempty(time_text)
                title( sprintf('Component %d' , nc) , 'FontSize',fs+2)
            else
                title( time_text{nc} , 'FontSize',fs+2)
            end
        elseif ~isempty(TIMESERIES) && show_psd 
            %% Plot power spectrum density
            pwelch(squeeze(TIMESERIES(nc,:,:)),[],[],[],TR^(-1));
            set(gca,'Fontsize',fs*.9)
        end
        
        
    elseif show_hist && ~histogramInSlice 
        %% Only show histogram, not timeseries
        if strcmp(place_time,'below')
            subplot(vplots+1,hplots,((vplots*hplots) +1):((vplots+1)*hplots))
        elseif strcmp(place_time,'right')
            time_idx = bsxfun(@plus,(hslices+1):(hplots),...
                    bsxfun(@times,ones(vplots,ceil(hslices*slice_time_ratio))...
                    ,(0:(vplots-1))')*hplots);
            subplot(vplots,hplots,time_idx(:))
        end
        illustrateHistogram(MAPS(:,nc),nbins,hist_method,fs);
    end
    
    %Show Y-axis is to the right, avoids overlab with plotted slices
    if strcmp(place_time,'right') && (show_hist || ~isempty(TIMESERIES))
            set(gca,'YAxisLocation','right')
    end
    
    if ~isempty(super_title)
        suptitle(super_title);
    end
    
    if ~isempty(save_loc)
        h.PaperPositionMode = 'auto';
        print(sprintf('%s_d%i',save_loc,nc),'-dpng')
        print(sprintf('%s_d%i',save_loc,nc),'-depsc')
        if mod(nc,25) == 0
            warning('Closing all figures to avoid heap space error')
            close all
            %h.delete % removes figure (to avoid heap space exceeded)
        end
    end
end

%eof
end

function illustrateHistogram(MAP,nbins,hist_method,fontsize)
%# create histogram from a normal distribution.
[counts,binlimits]=hist(MAP,nbins);

%Illustrate the histogram under the following constraint
if strcmpi(hist_method,'frequency')
    % frequencies, or bin counts
    bar(binlimits,counts);
    ylabel('Frequency');
elseif strcmpi(hist_method,'proportion')
    bar(binlimits,counts/sum(counts)); 
    ylabel('Proportions');
elseif strcmpi(hist_method,'area')
    bar(binlimits,counts/trapz(binlimits,counts));
    ylabel('Proportions (area)');
elseif strcmpi(hist_method,'fdf')
    %A frequency density estimate which is distribution independent
    bar(binlimits,counts/abs(diff(bin(limits(1:2)))))
elseif strcmpi(hist_method,'pdf')
    % probability density, i.e. frequency density/total count
    bar(binlimits,(counts/abs(diff(binlimits(1:2))))/sum(counts));
    ylabel('Prob. density')
end

xlim = max(abs(get(gca,'Xtick')));
ylim = max(abs(get(gca,'Ytick')));
axis([-xlim,xlim,0,ylim]);
set(gca,'Fontsize',ceil(fontsize))
xlabel('Value');

text(-xlim*0.9,ylim*0.8,...
        sprintf('k=%2.2f',kurtosis(MAP)-3),'fontsize',fontsize*1.1)

end
