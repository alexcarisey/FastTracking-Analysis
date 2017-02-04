%% Custom script to calculate the speed, velocities, distances, displacement, RoG and MSD of tracked objects in Imaris

% This script takes the Imaris export csv file where all the measurement
% are listed time point after time point and re-shuffle them into collated
% columns with one track per column

% This code allows the user to select a folder of csv files and then the
% column inside the csv file that will be collated into a single
% spreadsheet with one file per column.

%% Cleanup

close all
clearvars
home

%% User variables

frame_rate = 0.011; % Time interval betweent the frames in second
SPACE_UNITS = 'µm'; % Space unit in which the measurements are imported
TIME_UNITS = 's'; % Time unit in which the measurements are imported
minimum_track_length = 100; % Threshold for including a track in the analysis (this is a frame number, not a duration!)
every_n_frame = 1; % To use if you want to change the sampling rate
value_for_fitLogLogMSD = 0.25; % This is the portion of the MSD curve used for the fit
R2LIMIT = 0.75; % Threshold value for what is considered a good fit
boundary_alpha_confined = 0.25; % Limit for alpha value between confined and diffusion
boundary_alpha_directed = 1.5; % Limit for alpha value between diffusion and directed

%% Loading and various parameters

% Header for the POSITION import table is expected as follows:
% Position X	Position Y	Position Z	Unit	Category	Collection	Birth [s]	Death [s]	TrackID	ID	OriginalID	Original Component Name	Original Component ID	Original Image Name	Original Image ID
% Convert as:
% {'PositionX','PositionY','PositionZ','Unit','Category','Collection','Birth','Death','TrackID','ID','OriginalID','OriginalComponentName','OriginalComponentID','OriginalImageName','OriginalImageID'}

%POOL = parpool('local',8);
tic

disp('FastTracker - v1.0 - 2017-02-02');
disp(' ');
disp('Loading source csv file');

[filename, path, ~] = uigetfile('.csv');
delimiter = ',';
startRow = 4; % note that sometimes, Imaris adds a padding line on the top!? So start at 5 if textscan encounters an error
formatSpec = '%f%f%f%s%s%s%f%f%s%s%s%s%s%s%s%[^\n\r]'; % Important: 'PositionX','PositionY','PositionZ','Birth','Death' are imported as a number, rest as strings
fileID = fopen([path filename],'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);
fclose(fileID);
Extracted_data = table(dataArray{1:end-1}, 'VariableNames', {'PositionX','PositionY','PositionZ','Unit','Category','Collection','Birth','Death','TrackID','ID','OriginalID','OriginalComponentName','OriginalComponentID','OriginalImageName','OriginalImageID'});
disp(['File: ', filename]);
disp(['Total number of points: ', mat2str(height(Extracted_data))]);
clearvars delimiter startRow formatSpec fileID dataArray ans;

[~,filename,~] = fileparts(filename);

%% Location to save the files
%pathname = uigetdir('~/Desktop/', 'Directory where csv files will be stored');

%% Foolproofing the TrackID in case they are not unique by merging OriginalImageName and TrackID into TruelyUniqueName (16th column)

C_single = cell(size(Extracted_data,1),1); % Pre-allocation for speed
C_multiple = [Extracted_data.OriginalImageName Extracted_data.TrackID];
textprogressbar('Foolproofing the TrackIDs:   ');
indicator_progress = 0;
for hh = 1:size(C_multiple,1)
    C_single{hh,1} = strjoin(C_multiple(hh,:));
    indicator_progress = indicator_progress + (100/size(C_multiple,1));
    textprogressbar(indicator_progress);
end
textprogressbar(' Done');
C_single_table = cell2table(C_single,'VariableNames',{'TruelyUniqueName'});
Extracted_data = [Extracted_data C_single_table];
clear C_* hh indicator_progress

%% Creating a smaller table with only a few columns: {'PositionX','PositionY','PositionZ'}

disp('Extracting XYZ positions array');
Extracted_limited_data = table(Extracted_data.PositionX,Extracted_data.PositionY,Extracted_data.PositionZ,'VariableNames',{'X','Y','Z'});

%% Separate the main table into cells using their TruelyUniqueName

[TUN, TUN_idx_last, TUN_idx] = unique(Extracted_data(:,16),'stable');
disp(['Number of tracks: ', mat2str(height(TUN))]);
unique_idx = accumarray(TUN_idx(:),(1:length(TUN_idx))',[],@(x) {sort(x)});
textprogressbar('Re-constructing the tracks:   ');
indicator_progress = 0;
all_the_tracks_all_values = cell(height(TUN),1); % Pre-allocation for speed
all_the_tracks_time_position = cell(height(TUN),1); % Pre-allocation for speed
for jj = 1:height(TUN)
    temporary_table_1 = Extracted_data(unique_idx{jj},:);
    timeline_1 = array2table(linspace(0,(height(temporary_table_1)-1)*frame_rate,height(temporary_table_1))','VariableNames',{'Time'});
    temporary_table_2 = [timeline_1 temporary_table_1];
    all_the_tracks_all_values{jj} = temporary_table_2;
    temporary_table_3 = table2array(Extracted_limited_data(unique_idx{jj},:));
    timeline_2 = linspace(0,(height(temporary_table_2)-1)*frame_rate,height(temporary_table_2))';
    temporary_table_4 = [timeline_2 temporary_table_3];
    all_the_tracks_time_position{jj} = temporary_table_4;
    clear temporary_table_* timeline_*
    indicator_progress = indicator_progress + (100/height(TUN));
    textprogressbar(indicator_progress);
end
textprogressbar(' Done');
clear jj indicator_progress unique_idx TUN*

%% Adding a threshold parameter for the minimum number of frames per track and plot the distribution

disp(['Minimum track duration: ', mat2str(minimum_track_length), ' frames (', mat2str(minimum_track_length*frame_rate),' ',TIME_UNITS,').']);
index_minimum_length = false(length(all_the_tracks_all_values),1); % Pre-allocation for speed
track_duration = zeros(length(all_the_tracks_all_values),1); % Pre-allocation for speed

textprogressbar('Removing all the short tracks:   ');
indicator_progress = 0;
for i=1:length(all_the_tracks_all_values)
    index_minimum_length(i,1) = (height(all_the_tracks_all_values{i,1}) > minimum_track_length);
    track_duration(i,1) = height(all_the_tracks_all_values{i,1});
    indicator_progress = indicator_progress + (100/length(all_the_tracks_all_values));
    textprogressbar(indicator_progress);
end
textprogressbar(' Done');

mini_length_tracks_all_values = all_the_tracks_all_values(index_minimum_length);
mini_length_tracks_time_position = all_the_tracks_time_position(index_minimum_length);

disp([mat2str((length(all_the_tracks_all_values)-length(mini_length_tracks_all_values))),' tracks out of ', mat2str(length(all_the_tracks_all_values)),' were deleted. ',mat2str(length(mini_length_tracks_all_values)),' tracks left.']);
disp(['Longest track is ', mat2str(max(track_duration)),' frames long (', mat2str(max(track_duration)*frame_rate),' ',TIME_UNITS,').']);

clear i indicator_progress index_minimum_length

%% Re-sampling the acquisition frequency

skipping_tracks_all_values = cell(length(mini_length_tracks_all_values),1); % Pre-allocation for speed
skipping_tracks_time_position = cell(length(mini_length_tracks_time_position),1); % Pre-allocation for speed
for track = 1:length(mini_length_tracks_all_values)
    skipping_tracks_all_values{track,1} = mini_length_tracks_all_values{track,1}(1:every_n_frame:end,:);
    skipping_tracks_time_position{track,1} = mini_length_tracks_time_position{track,1}(1:every_n_frame:end,:);
end

clear mini_length_tracks_*
mini_length_tracks_all_values = skipping_tracks_all_values;
mini_length_tracks_time_position = skipping_tracks_time_position;

%% Calculating additional parameters using only the X,Y,Z positions and time interval

disp('Pre-allocation');
textprogressbar('Progress:   ');
indicator_progress = 0;
for track = 1:length(mini_length_tracks_all_values)
    nrow = size(mini_length_tracks_all_values{track},1);
    mini_length_tracks_all_values{track}.PTPDistance = zeros(nrow, 1);
    mini_length_tracks_all_values{track}.CumulativeDistance = zeros(nrow, 1);
    mini_length_tracks_all_values{track}.Displacement = zeros(nrow, 1);
    mini_length_tracks_all_values{track}.InstantVelocity = zeros(nrow, 1);
    mini_length_tracks_all_values{track}.CumulativeSpeed = zeros(nrow, 1);
    mini_length_tracks_all_values{track}.AverageSpeed = zeros(nrow, 1);
    mini_length_tracks_all_values{track}.RoG_ComponentX = zeros(nrow, 1);
    mini_length_tracks_all_values{track}.RoG_ComponentY = zeros(nrow, 1);
    mini_length_tracks_all_values{track}.RoG_ComponentZ = zeros(nrow, 1);
    mini_length_tracks_all_values{track}.RadiusOfGyration2 = zeros(nrow, 1);
    indicator_progress = indicator_progress + (100/length(mini_length_tracks_all_values));
    textprogressbar(indicator_progress);
end
textprogressbar(' Done');
clear indicator_progress track nrow

disp('Computing distances, displacement and velocities...');
textprogressbar('Progress:   ');
indicator_progress = 0;
for track = 1:length(mini_length_tracks_all_values)
    for position = 2:(height(mini_length_tracks_all_values{track}))
        mini_length_tracks_all_values{track}.PTPDistance(position) = sqrt((mini_length_tracks_all_values{track}.PositionX(position)-mini_length_tracks_all_values{track}.PositionX(position-1))^2 + (mini_length_tracks_all_values{track}.PositionY(position)-mini_length_tracks_all_values{track}.PositionY(position-1))^2 + (mini_length_tracks_all_values{track}.PositionZ(position)-mini_length_tracks_all_values{track}.PositionZ(position-1))^2);
        mini_length_tracks_all_values{track}.CumulativeDistance(position) = mini_length_tracks_all_values{track}.PTPDistance(position) + mini_length_tracks_all_values{track}.CumulativeDistance(position-1);
        mini_length_tracks_all_values{track}.Displacement(position) = sqrt((mini_length_tracks_all_values{track}.PositionX(end)-mini_length_tracks_all_values{track}.PositionX(1))^2 + (mini_length_tracks_all_values{track}.PositionY(end)-mini_length_tracks_all_values{track}.PositionY(1))^2 + (mini_length_tracks_all_values{track}.PositionZ(end)-mini_length_tracks_all_values{track}.PositionZ(1))^2);
        mini_length_tracks_all_values{track}.InstantVelocity(position) = mini_length_tracks_all_values{track}.PTPDistance(position) / frame_rate;
        mini_length_tracks_all_values{track}.CumulativeSpeed(position) = mini_length_tracks_all_values{track}.CumulativeDistance(position) / (frame_rate * (position-1));
    end
    indicator_progress = indicator_progress + (100/length(mini_length_tracks_all_values));
    textprogressbar(indicator_progress);
end
textprogressbar(' Done');
clear indicator_progress track position nrow

disp('Computing components for radius of gyration...');
textprogressbar('Progress:   ');
indicator_progress = 0;
for track = 1:length(mini_length_tracks_all_values)
    for position = 1:(height(mini_length_tracks_all_values{track}))
        mini_length_tracks_all_values{track}.RoG_ComponentX(position) = (mini_length_tracks_all_values{track}.PositionX(position) - mean(mini_length_tracks_all_values{track}.PositionX(:)))^2;
        mini_length_tracks_all_values{track}.RoG_ComponentY(position) = (mini_length_tracks_all_values{track}.PositionY(position) - mean(mini_length_tracks_all_values{track}.PositionY(:)))^2;
        mini_length_tracks_all_values{track}.RoG_ComponentZ(position) = (mini_length_tracks_all_values{track}.PositionZ(position) - mean(mini_length_tracks_all_values{track}.PositionZ(:)))^2;
    end
    indicator_progress = indicator_progress + (100/length(mini_length_tracks_all_values));
    textprogressbar(indicator_progress);
end
textprogressbar(' Done');
clear indicator_progress track position nrow

disp('Computing radius of gyration^2 and average speed...');
textprogressbar('Progress:   ');
indicator_progress = 0;
for track = 1:length(mini_length_tracks_all_values)
    mini_length_tracks_all_values{track}.AverageSpeed(:) = mean(mini_length_tracks_all_values{track}.InstantVelocity(:));
    mini_length_tracks_all_values{track}.RadiusOfGyration2(:) = ( sum(mini_length_tracks_all_values{track}.RoG_ComponentX(:)) + sum(mini_length_tracks_all_values{track}.RoG_ComponentY(:)) + sum(mini_length_tracks_all_values{track}.RoG_ComponentZ(:))) / height(mini_length_tracks_all_values{track});
    indicator_progress = indicator_progress + (100/length(mini_length_tracks_all_values));
    textprogressbar(indicator_progress);
end
textprogressbar(' Done');
clear indicator_progress track position nrow

%% @msdanalyzer class loading and computing

ma = msdanalyzer(3, 'µm', 's'); % Build a new MSD analyzer object.
ma = ma.addAll(mini_length_tracks_time_position); % Load tracks to msdanalyzer.
ma = ma.computeMSD; % Compute MSD
disp(['Fitting of MSD is currently set to be done with a factor of ',mat2str(((max(track_duration)*frame_rate/(1*frame_rate))*value_for_fitLogLogMSD)/every_n_frame),' between smallest delay and the largest one.']);
if ((((max(track_duration)*frame_rate/(1*frame_rate))*value_for_fitLogLogMSD)/every_n_frame)<=1000)
    disp('Power laws can only be determined reliably if they are sampled over at least 3 orders of magnitude.');
    disp('This means that there must be at least a factor of 1000 between the smallest delay and the largest one.');
end
ma = ma.fitLogLogMSD(value_for_fitLogLogMSD); % Fit the individual MSD and store the alpha (slope), gamma (value at origin) and r2 (fitness) into the object

%% Alpha values and statistical analysis

disp('Extracting all the alpha/gamma/r2fit values from the MSD curves fit.');
alphas_MSD = ma.loglogfit.alpha;
gammas_MSD = ma.loglogfit.gamma;
r2fits_MSD = ma.loglogfit.r2fit;
bad_fits = r2fits_MSD < R2LIMIT; % Remove bad fits
fprintf('Keeping %d fits out of %d (R2 > %.2f).\n', sum(~bad_fits), length(r2fits_MSD), R2LIMIT);
alphas_MSD(bad_fits) = [];
[htest, pval] = ttest(alphas_MSD, 1, 0.05, 'left'); % T-test
if ~htest
    [htest, pval] = ttest(alphas_MSD, 1, 0.05);
end

% Prepare string
str = { [ '\alpha = ' sprintf('%.2f ± %.2f (mean ± std, N = %d)', mean(alphas_MSD), std(alphas_MSD), numel(alphas_MSD)) ] };
if htest
    str{2} = sprintf('Significantly below 1, with p = %.2g', pval);
else
    str{2} = sprintf('Not significantly differend from 1, with p = %.2g', pval);
end

% Creating indexes according to the alpha value and the goodness of the fit of the MSD
index_confined_MSD = ma.loglogfit.alpha < boundary_alpha_confined & ma.loglogfit.r2fit > R2LIMIT; % Confined
index_diffusive_MSD = boundary_alpha_confined <= ma.loglogfit.alpha & ma.loglogfit.alpha <= boundary_alpha_directed & ma.loglogfit.r2fit > R2LIMIT; % Diffusive
index_directed_MSD = ma.loglogfit.alpha > boundary_alpha_directed & ma.loglogfit.r2fit > R2LIMIT; % Directed
index_bad_fits = ma.loglogfit.r2fit <= R2LIMIT; % Bad fits

% Data preparation for the plots
array_categories_plot = [sum(index_confined_MSD) sum(index_diffusive_MSD) sum(index_directed_MSD)]; % For the "Mode of Motion" plot
average_speeds = zeros(length(mini_length_tracks_all_values),1); % Pre-allocation for the Granule velocity (mean) plot
radius_of_gyration2 = zeros(length(mini_length_tracks_all_values),1); % Pre-allocation for the Radius of gyration plot
for track = 1:length(mini_length_tracks_all_values)
    average_speeds(track,1) = mini_length_tracks_all_values{track}.AverageSpeed(1); % For the Granule velocity (mean) plot
    radius_of_gyration2(track,1) = mini_length_tracks_all_values{track}.RadiusOfGyration2(1); % For the Radius of gyration plot
end

%% Compute velocity autocorrelations

% ma = ma.computeVCorr; % Compute velocity autocorrelation

%% MEGA PLOT

figure1 = figure('rend','painters','pos',[800 300 1800 1000]);

subplot(2,7,[1 2])
hold on
for track = 1:length(mini_length_tracks_all_values)
    plot(mini_length_tracks_all_values{track}.PositionX,mini_length_tracks_all_values{track}.PositionY)
end
hold off
box on
axis equal
axis square
title('Trajectories','FontSize', 16,'FontName','Arial')
xlabel(SPACE_UNITS,'FontName','Arial')
ylabel(SPACE_UNITS,'FontName','Arial')

subplot(2,7,3)
hold on
plotSpread(average_speeds);
boxplot(average_speeds);
hold off
title('Granule velocity (mean)','FontSize', 16,'FontName','Arial')
ylabel([SPACE_UNITS,'/',TIME_UNITS],'FontName','Arial')
set(gca,'XTickLabel','')

subplot(2,7,4);
hold on
plotSpread(radius_of_gyration2);
boxplot(radius_of_gyration2);
hold off
title('Radius of gyration','FontSize', 16,'FontName','Arial')
ylabel([SPACE_UNITS,'^2'],'FontName','Arial')
set(gca,'XTickLabel','')

subplot(2,7,[5 7])
hold on
rectangle('position',[0 0 0.25 0.14],'facecolor',[0.729411780834198 0.831372559070587 0.95686274766922]);
rectangle('position',[0.25 0 1.25 0.14],'facecolor',[0.756862759590149 0.866666674613953 0.776470601558685]);
rectangle('position',[1.5 0 0.5 0.14],'facecolor',[0.925490200519562 0.839215695858002 0.839215695858002]);
histogram(alphas_MSD,40,'Normalization','probability');
xlabel('\alpha','FontName','Arial')
ylabel('Frequency','FontName','Arial')
text(0.8, 0.06, str, ...
    'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'bottom', ...
    'FontSize', 14)
title('\alpha values distribution','FontSize', 16,'FontName','Arial')
hold off
box on;

subplot(2,7,[8 9])
hold on
rectangle('position',[0 0.00001 minimum_track_length 1],'facecolor',[0.9 0.8 0.8]);
histogram(track_duration,100,'Normalization','probability');
box on;
axes1=gca;
set(axes1,'YMinorTick','on','YScale','log');
line([minimum_track_length minimum_track_length],get(axes1,'YLim'),'Color',[1 0 0],'LineWidth',2)
hold off
title('Distribution of track lengths','FontSize', 16,'FontName','Arial')
xlabel('Duration (Number of frames)','FontName','Arial')
ylabel('Frequency (log scale)','FontName','Arial')

subplot(2,7,[10 11]);
histogram(average_speeds,40,'Normalization','probability')
title('Distribution of mean granule velocity','FontSize', 16,'FontName','Arial')
xlabel(['Speed (',SPACE_UNITS,'/',TIME_UNITS,')'],'FontName','Arial');
ylabel('Frequency','FontName','Arial');

subplot(2,7,[12 13]);
histogram(radius_of_gyration2,40,'Normalization','probability')
title('Distribution of radius of gyration','FontSize', 16,'FontName','Arial')
xlabel(['Radius of gyration (',SPACE_UNITS,'^2)'],'FontName','Arial');
ylabel('Frequency','FontName','Arial')

subplot(2,7,14)
bar(array_categories_plot,'FaceColor',[0.4 0.67 0.84]);
box on
label_modes = {'Confined','Diffusive','Directed'};
set( gca(), 'XTickLabel', label_modes, 'FontName','Arial')
rotateXLabels( gca(), 45 )
ylabel('Count','FontName','Arial')
title('Modes of motion','FontSize', 16,'FontName','Arial')

suptitle(filename);

save = [filename '_fig1.eps'];
print('-painters','-depsc','-loose',save);
clear track
close all

%%

% meanmsd = ma.getMeanMSD; % Compute the weighted mean of all MSD curves. 
% meanvcorr = ma.getMeanVCorr; % Compute the weighted mean of velocity autocorrelation. 
% instant_velocities = ma.getVelocities; % Generate and return the instantaneous velocities. 
% ma = ma.fitMSD; % Fit all MSD curves by a linear function.

%% Plotting part 2

% Plotting and curve fitting

figure2 = figure('rend','painters','pos',[800 300 2000 900]);

subplot(2,5,1)
[hps1, ha1] =  ma.plotMSD; % Plot the mean square displacement curves.
ma.labelPlotMSD(ha1); % A convenience method to set the axes labels.
xlim(ha1,[0 10]); ylim(ha1,[0 0.1]); title('MSD for all tracks','FontSize', 16,'FontName','Arial'); box on; axis square

if sum(index_bad_fits) ~= 0 % Required because if the index is empty, plotMSD plot ALL the MSD values by default!
    subplot(2,5,2)
    [hps2, ha2] =  ma.plotMSD(gca,find(index_bad_fits)); % Plot the mean square displacement curves.
    ma.labelPlotMSD(ha2); % A convenience method to set the axes labels.
    xlim(ha2,[0 10]); ylim(ha2,[0 0.1]); title(['MSD for tracks with r2<',mat2str(R2LIMIT)],'FontSize', 16,'FontName','Arial'); box on; axis square
end

if sum(index_confined_MSD) ~= 0 % Required because if the index is empty, plotMSD plot ALL the MSD values by default!
    subplot(2,5,3)
    [hps3, ha3] =  ma.plotMSD(gca,find(index_confined_MSD)); % Plot the mean square displacement curves.
    ma.labelPlotMSD(ha3); % A convenience method to set the axes labels.
    xlim(ha3,[0 10]); ylim(ha3,[0 0.1]); title('MSD for tracks \alpha<0.25','FontSize', 16,'FontName','Arial'); box on; axis square
end

if sum(index_diffusive_MSD) ~= 0 % Required because if the index is empty, plotMSD plot ALL the MSD values by default!
    subplot(2,5,4)
    [hps4, ha4] =  ma.plotMSD(gca,find(index_diffusive_MSD)); % Plot the mean square displacement curves.
    ma.labelPlotMSD(ha4); % A convenience method to set the axes labels.
    xlim(ha4,[0 10]); ylim(ha4,[0 0.1]); title('MSD for tracks 0.25<\alpha<1.5','FontSize', 16,'FontName','Arial'); box on; axis square
end

if sum(index_directed_MSD) ~= 0 % Required because if the index is empty, plotMSD plot ALL the MSD values by default!
    subplot(2,5,5)
    [hps5, ha5] =  ma.plotMSD(gca,find(index_directed_MSD)); % Plot the mean square displacement curves.
    ma.labelPlotMSD(ha5); % A convenience method to set the axes labels.
    xlim(ha5,[0 10]); ylim(ha5,[0 0.1]); title('MSD for tracks \alpha>1.5','FontSize', 16,'FontName','Arial'); box on; axis square
end

subplot(2,5,6)
[hps6, ha6] =  ma.plotMeanMSD(gca,true); % Plot the MEAN mean square displacement curves.
ma.labelPlotMSD(ha6); % A convenience method to set the axes labels.
xlim(ha6,[0 10]); ylim(ha6,[0 0.1]); title('Mean MSD for all tracks','FontSize', 16,'FontName','Arial'); box on; axis square

if sum(index_bad_fits) ~= 0 % Required because if the index is empty, plotMSD plot ALL the MSD values by default!
    subplot(2,5,7)
    [hps7, ha7] =  ma.plotMeanMSD(gca,true,find(index_bad_fits)); % Plot the mean square displacement curves.
    ma.labelPlotMSD(ha7); % A convenience method to set the axes labels.
    xlim(ha7,[0 10]); ylim(ha7,[0 0.1]); title(['MSD for tracks with r2<',mat2str(R2LIMIT)],'FontSize', 16,'FontName','Arial'); box on; axis square
end

if sum(index_confined_MSD) ~= 0 % Required because if the index is empty, plotMSD plot ALL the MSD values by default!
    subplot(2,5,8)
    [hps8, ha8] =  ma.plotMeanMSD(gca,true,find(index_confined_MSD)); % Plot the MEAN mean square displacement curves.
    ma.labelPlotMSD(ha8); % A convenience method to set the axes labels.
    xlim(ha8,[0 10]); ylim(ha8,[0 0.1]); title('Mean MSD for tracks \alpha<0.25','FontSize', 16,'FontName','Arial'); box on; axis square
end

if sum(index_diffusive_MSD) ~= 0 % Required because if the index is empty, plotMSD plot ALL the MSD values by default!
    subplot(2,5,9)
    [hps9, ha9] =  ma.plotMeanMSD(gca,true,find(index_diffusive_MSD)); % Plot the MEAN mean square displacement curves.
    ma.labelPlotMSD(ha9); % A convenience method to set the axes labels.
    xlim(ha9,[0 10]); ylim(ha9,[0 0.1]); title('Mean MSD for tracks 0.25<\alpha<1.5','FontSize', 16,'FontName','Arial'); box on; axis square
end

if sum(index_directed_MSD) ~= 0 % Required because if the index is empty, plotMSD plot ALL the MSD values by default!
    subplot(2,5,10)
    [hps10, ha10] =  ma.plotMeanMSD(gca,true,find(index_directed_MSD)); % Plot the MEAN mean square displacement curves.
    ma.labelPlotMSD(ha10); % A convenience method to set the axes labels.
    xlim(ha10,[0 10]); ylim(ha10,[0 0.1]); title('Mean MSD for tracks \alpha>1.5','FontSize', 16,'FontName','Arial'); box on; axis square
end

suptitle(filename);

save = [filename '_fig2.eps'];
print('-painters','-depsc','-loose',save);

clear track
close all

%% Clean up

clear ha* hps* hmsd fo gof ans

toc