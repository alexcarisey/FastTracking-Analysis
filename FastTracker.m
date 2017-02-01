%% Custom script to calculate the speed and MSD of tracked objects in Imaris

% This script takes the Imaris export csv file where all the measurement
% are listed time point after time point and re-shuffle them into collated
% columns with one track per column

% This code allows the user to select a folder of csv files and then the
% column inside the csv file that will be collated into a single
% spreadsheet with one file per column.

%% Loading and various parameters

% Header for the POSITION import table is as follows:
% Position X	Position Y	Position Z	Unit	Category	Collection	Birth [s]	Death [s]	TrackID	ID	OriginalID	Original Component Name	Original Component ID	Original Image Name	Original Image ID
% Convert as:
% {'PositionX','PositionY','PositionZ','Unit','Category','Collection','Birth','Death','TrackID','ID','OriginalID','OriginalComponentName','OriginalComponentID','OriginalImageName','OriginalImageID'}

[filename, path, ~] = uigetfile('.csv');
delimiter = ',';
startRow = 4; % note that sometimes, Imaris adds a padding line on the top!? So start at 5 if textscan encounters an error
formatSpec = '%f%f%f%s%s%s%f%f%s%s%s%s%s%s%s%[^\n\r]'; % Important: 'PositionX','PositionY','PositionZ','Birth','Death' are imported as a number, rest as strings
fileID = fopen([path filename],'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);
fclose(fileID);
Extracted_data = table(dataArray{1:end-1}, 'VariableNames', {'PositionX','PositionY','PositionZ','Unit','Category','Collection','Birth','Death','TrackID','ID','OriginalID','OriginalComponentName','OriginalComponentID','OriginalImageName','OriginalImageID'});

clearvars filename path delimiter startRow formatSpec fileID dataArray ans;

frame_rate = 0.011;
SPACE_UNITS = 'µm';
TIME_UNITS = 's';

%% Location to save the files

%pathname = uigetdir('~/Desktop/', 'Directory where csv files will be stored');

%% Foolproofing the TrackID in case they are not unique by merging OriginalImageName and TrackID into TruelyUniqueName (16th column)

C_single = cell(size(Extracted_data,1),1); % Pre-allocation for speed
C_multiple = [Extracted_data.OriginalImageName Extracted_data.TrackID];
for hh = 1:size(C_multiple,1)
    C_single{hh,1} = strjoin(C_multiple(hh,:));
end
C_single_table = cell2table(C_single,'VariableNames',{'TruelyUniqueName'});
Extracted_data = [Extracted_data C_single_table];
clear C_* hh

%% Creating a smaller table with only a few columns: {'PositionX','PositionY','PositionZ'}

Extracted_limited_data = table(Extracted_data.PositionX,Extracted_data.PositionY,Extracted_data.PositionZ,'VariableNames',{'X','Y','Z'});

%% Separate the main table into cells using their TruelyUniqueName

[TUN, TUN_idx_last, TUN_idx] = unique(Extracted_data(:,16),'stable');
disp(['Number of trajectories: ', mat2str(height(TUN))]);
unique_idx = accumarray(TUN_idx(:),(1:length(TUN_idx))',[],@(x) {sort(x)});
textprogressbar('Progress:   ');
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
indicator_progress = 100;
textprogressbar(' Done');
clear jj indicator_progress pathname unique_idx TUN*

%% Calculating additional parameters using only the X,Y,Z positions and time interval

% General pre-allocation for the additional columns

disp('Pre-allocation');
textprogressbar('Progress:   ');
indicator_progress = 0;
for track = 1:length(all_the_tracks_all_values)
    nrow = size(all_the_tracks_all_values{track},1);
    all_the_tracks_all_values{track}.PTPDistance = zeros(nrow, 1);
    all_the_tracks_all_values{track}.CumulativeDistance = zeros(nrow, 1);
    all_the_tracks_all_values{track}.Displacement = zeros(nrow, 1);
    all_the_tracks_all_values{track}.InstantVelocity = zeros(nrow, 1);
    all_the_tracks_all_values{track}.CumulativeSpeed = zeros(nrow, 1);
    all_the_tracks_all_values{track}.AverageSpeed = zeros(nrow, 1);
    indicator_progress = indicator_progress + (100/length(all_the_tracks_all_values));
    textprogressbar(indicator_progress);
end
textprogressbar(' Done');

clear indicator_progress track nrow

% Populate the additional fields

disp('Computing distances, displacement, speed and velocities...');
textprogressbar('Progress:   ');
indicator_progress = 0;
for track = 1:length(all_the_tracks_all_values)
    for position = 2:(height(all_the_tracks_all_values{track}))
        all_the_tracks_all_values{track}.PTPDistance(position) = sqrt((all_the_tracks_all_values{track}.PositionX(position)-all_the_tracks_all_values{track}.PositionX(position-1))^2 + (all_the_tracks_all_values{track}.PositionY(position)-all_the_tracks_all_values{track}.PositionY(position-1))^2 + (all_the_tracks_all_values{track}.PositionZ(position)-all_the_tracks_all_values{track}.PositionZ(position-1))^2);
        all_the_tracks_all_values{track}.CumulativeDistance(position) = all_the_tracks_all_values{track}.PTPDistance(position) + all_the_tracks_all_values{track}.CumulativeDistance(position-1);
        all_the_tracks_all_values{track}.Displacement(position) = sqrt((all_the_tracks_all_values{track}.PositionX(end)-all_the_tracks_all_values{track}.PositionX(1))^2 + (all_the_tracks_all_values{track}.PositionY(end)-all_the_tracks_all_values{track}.PositionY(1))^2 + (all_the_tracks_all_values{track}.PositionZ(end)-all_the_tracks_all_values{track}.PositionZ(1))^2);
        all_the_tracks_all_values{track}.InstantVelocity(position) = all_the_tracks_all_values{track}.PTPDistance(position) / frame_rate;
        all_the_tracks_all_values{track}.CumulativeSpeed(position) = all_the_tracks_all_values{track}.CumulativeDistance(position) / (frame_rate * (position-1));
        all_the_tracks_all_values{track}.AverageSpeed(:) = mean(all_the_tracks_all_values{track}.InstantVelocity(:));
    end
    indicator_progress = indicator_progress + (100/length(all_the_tracks_all_values));
    textprogressbar(indicator_progress);
end
indicator_progress = 100;
textprogressbar(' Done');
clear indicator_progress track position nrow

%% @msdanalyzer 

% Initilization

ma = msdanalyzer(3, 'µm', 's'); % Builds a new MSD analyzer object.
ma = ma.addAll(all_the_tracks_time_position); % Add specified trajectories to msdanalyzer.

% Compute

ma = ma.computeMSD; % Compute the mean-squared-displacement for this object.
ma = ma.computeVCorr; % Compute velocity autocorrelation.

meanmsd = ma.getMeanMSD; % Compute the weighted mean of all MSD curves. 
meanvcorr = ma.getMeanVCorr; % Compute the weighted mean of velocity autocorrelation. 
instant_velocities = ma.getVelocities; % Generate and return the instantaneous velocities. 

% Plot trajectories

figure
[hps1, ha1] = ma.plotTracks; % Plot the tracks stored in this object.
ma.labelPlotTracks(ha1); % A convenience method to set the axes labels.

% Plotting and curve fitting

figure
[hps2, ha2] =  ma.plotMSD; % Plot the mean square displacement curves. 
ma.labelPlotMSD(ha2); % A convenience method to set the axes labels. 

ma = ma.fitMSD( 0.5 ); % Fit all MSD curves by a linear function.

% Plot Mean MSD

figure
hmsd =  ma.plotMeanMSD(gca, true); % Plot the weighted mean of the MSD curves.
[fo, gof] = ma.fitMeanMSD( 0.1 ); % Fit the weighted averaged MSD by a linear function.
plot(fo)
legend off
ma.labelPlotMSD

ma = ma.fitLogLogMSD(0.5); % Fit the log-log MSD to determine behavior.
ma.loglogfit

mean(ma.loglogfit.alpha)

%% MeanVCorr

figure
[hps4, ha4] =  ma.plotMeanVCorr; % Plot the weighted mean of the velocity autocorrelation curves. 
ma.labelPlotVCorr(ha4); % A convenience method to set the axes labels.
