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
clc
home

%% Loading and various parameters

% Header for the POSITION import table is as follows:
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
clearvars filename path delimiter startRow formatSpec fileID dataArray ans;

frame_rate = 0.011;
SPACE_UNITS = 'µm';
TIME_UNITS = 's';
minimum_track_length = 100;

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
textprogressbar('Re-construction the tracks:   ');
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
clear jj indicator_progress pathname unique_idx TUN*

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
value_for_fitLogLogMSD = 0.50; % This is the portion of the MSD curve used for the fit
disp(['Fitting of MSD is currently set to be done with a factor of ',mat2str((max(track_duration)*frame_rate/(1*frame_rate))*value_for_fitLogLogMSD),' between smallest delay and the largest one.']);
if (((max(track_duration)*frame_rate/(1*frame_rate))*value_for_fitLogLogMSD)<=1000)
    disp('Power laws can only be determined reliably if they are sampled over at least 3 orders of magnitude.');
    disp('This means that there must be at least a factor of 1000 between the smallest delay and the largest one.');
end
ma = ma.fitLogLogMSD(value_for_fitLogLogMSD); % Fit the individual MSD and store the alpha (slope), gamma (value at origin) and r2 (fitness) into the object

%% Alpha values and statistical analysis

disp('Extracting all the alpha/gamma/r2fit values from the MSD curves fit.');
alphas_MSD = ma.loglogfit.alpha;
gammas_MSD = ma.loglogfit.gamma;
r2fits_MSD = ma.loglogfit.r2fit;
R2LIMIT = 0.8;
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
% Split the data into the 3 categories
boundary_confined = 0.25;
boundary_directed = 1.5;
lower_cat = alphas_MSD(alphas_MSD<boundary_confined);
middle_cat = alphas_MSD(alphas_MSD>=boundary_confined & alphas_MSD<=boundary_directed);
top_cat = alphas_MSD(alphas_MSD>boundary_directed);
array_categories_plot = [length(lower_cat) length(middle_cat) length(top_cat)];

% Data preparation

average_speeds = zeros(length(mini_length_tracks_all_values),1);
radius_of_gyration2 = zeros(length(mini_length_tracks_all_values),1);
for track = 1:length(mini_length_tracks_all_values)
    average_speeds(track,1) = mini_length_tracks_all_values{track}.AverageSpeed(1);
    radius_of_gyration2(track,1) = mini_length_tracks_all_values{track}.RadiusOfGyration2(1);
end

%% MEGA PLOT

figure1 = figure('rend','painters','pos',[800 300 1800 950]);

subplot(2,7,[1 2])
hold on
for track = 1:length(mini_length_tracks_all_values)
    plot(mini_length_tracks_all_values{track}.PositionX,mini_length_tracks_all_values{track}.PositionY)
end
hold off
box on
axis equal
axis square
title('Trajectories','FontSize', 16)
xlabel(SPACE_UNITS)
ylabel(SPACE_UNITS)

subplot(2,7,3)
hold on
plotSpread(average_speeds);
boxplot(average_speeds);
hold off
title('Average speed','FontSize', 16)
ylabel([SPACE_UNITS,'/',TIME_UNITS])
set(gca,'XTickLabel','')

subplot(2,7,4);
hold on
plotSpread(radius_of_gyration2);
boxplot(radius_of_gyration2);
hold off
title('Radius of gyration','FontSize', 16)
ylabel([SPACE_UNITS,'^2'])
set(gca,'XTickLabel','')

subplot(2,7,[5 7])
hold on
histogram(alphas_MSD,40,'Normalization','probability');
xlabel('\alpha')
ylabel('Probability')
text(0.8, 0.06, str, ...
    'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'bottom', ...
    'FontSize', 14)
title('\alpha values distribution','FontSize', 16)
hold off
box on;

subplot(2,7,[8 9])
hold on
histogram(track_duration,100,'Normalization','probability');
box on;
axes1=gca;
set(axes1,'YMinorTick','on','YScale','log');
line([minimum_track_length minimum_track_length],get(axes1,'YLim'),'Color',[1 0 0],'LineWidth',2)
hold off
title('Distribution of track lengths','FontSize', 16)
xlabel('Duration (Number of frames)')
ylabel('Frequency (log scale)')

subplot(2,7,[10 11]);
histogram(average_speeds,40,'Normalization','probability')
title('Distribution of average speed','FontSize', 16)
xlabel(['Speed (',SPACE_UNITS,'/',TIME_UNITS,')']);
ylabel('Frequency');

subplot(2,7,[12 13]);
histogram(radius_of_gyration2,40,'Normalization','probability')
title('Distribution of radius of gyration','FontSize', 16)
xlabel(['Radius of gyration (',SPACE_UNITS,'^2)']);
ylabel('Frequency')

subplot(2,7,14)
bar(array_categories_plot,'FaceColor',[0.4 0.67 0.84]);
box on

label_modes = {'Confined','Diffusive','Directed'};
set( gca(), 'XTickLabel', label_modes )
rotateXLabels( gca(), 45 )
ylabel('Count')
title('Modes of motion','FontSize', 16)

print -painters -depsc figure1.eps
clear track
close all

%%

% ma = ma.computeVCorr;% Compute velocity autocorrelation

%%

meanmsd = ma.getMeanMSD; % Compute the weighted mean of all MSD curves. 
meanvcorr = ma.getMeanVCorr; % Compute the weighted mean of velocity autocorrelation. 
instant_velocities = ma.getVelocities; % Generate and return the instantaneous velocities. 

%% Plotting part 2

% Plot trajectories

figure;
[hps1, ha1] = ma.plotTracks; % Plot the tracks stored in this object.
ma.labelPlotTracks(ha1); % A convenience method to set the axes labels.
print -painters -depsc figure2.eps
close all

% Plotting and curve fitting

figure;
[hps2, ha2] =  ma.plotMSD; % Plot the mean square displacement curves. 
ma.labelPlotMSD(ha2); % A convenience method to set the axes labels. 
ma = ma.fitMSD( 0.5 ); % Fit all MSD curves by a linear function.
print -painters -depsc figure3.eps
close all

% Plot Mean MSD

figure;
hmsd =  ma.plotMeanMSD(gca, true); % Plot the weighted mean of the MSD curves.
[fo, gof] = ma.fitMeanMSD( 0.1 ); % Fit the weighted averaged MSD by a linear function.
plot(fo)
legend off
ma.labelPlotMSD
print -painters -depsc figure4.eps
close all

%%

ma = ma.fitLogLogMSD(0.5); % Fit the log-log MSD to determine behavior.
ma.loglogfit

mean(ma.loglogfit.alpha)

%% MeanVCorr

figure;
[hps4, ha4] =  ma.plotMeanVCorr; % Plot the weighted mean of the velocity autocorrelation curves. 
ma.labelPlotVCorr(ha4); % A convenience method to set the axes labels.
print -painters -depsc figure5.eps
close all

%% Clean up

clear ha* hps* hmsd fo gof ans

toc