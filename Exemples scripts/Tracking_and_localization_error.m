%%Impact of tracking and localization error.

clear all
close all

% Typical bad localization error. Large compared to typical displacement.
BAD_XY_TYPICAL_OFFSET = 0.2; % µm

SPACE_UNITS = 'µm';
TIME_UNITS = 's';
N_PARTICLES = 100;
N_TIME_STEPS = 500;
N_DIM = 2; % 2D
D  = 1e-3; % µm^2/s - diffusion coefficient
dT = 0.2; % s,
SIZE = 8; % µm

k = sqrt(2 * D * dT);
tracks = cell(N_PARTICLES, 1);

for i_spot = 1 : N_PARTICLES

    % Time
    time = (0 : N_TIME_STEPS-1)' * dT;

    % Initial position
    X0 = SIZE .* rand(1, N_DIM);

    % Integrate uncorrelated displacement
    dX = k * randn(N_TIME_STEPS, N_DIM);
    dX(1, :) = X0;
    X = cumsum(dX, 1);

    % Deal with incorrect detection
    bad_dx = BAD_XY_TYPICAL_OFFSET * randn(N_TIME_STEPS, N_DIM);
    X = X + bad_dx;

    % Store
    tracks{i_spot} = [ time X];

end

fprintf('Generated %d tracks over %d time steps.\n', N_PARTICLES, N_TIME_STEPS)

ma = msdanalyzer(2, SPACE_UNITS, TIME_UNITS);
ma = ma.addAll(tracks);
ma = ma.computeMSD;
ma.plotMeanMSD(gca, true);

ma = ma.fitMSD;

good_enough_fit = ma.lfit.r2fit > 0.8;
Dmean = mean( ma.lfit.a(good_enough_fit) ) / 2 / ma.n_dim;
Dstd  =  std( ma.lfit.a(good_enough_fit) ) / 2 / ma.n_dim;

fprintf('Found D = %.3e ± %.3e (mean ± std, N = %d)\n', ...
    Dmean, Dstd, sum(good_enough_fit));

fprintf('Simulated localization error magnitude: %.3e %s.\n', ...
    BAD_XY_TYPICAL_OFFSET, SPACE_UNITS);
fprintf('Typical brownian displacement magnitude: %.3e %s.\n', ...
    k, SPACE_UNITS);

bmean = mean( ma.lfit.b(good_enough_fit) );
sigma_locmean = 0.5 * sqrt(bmean);

% Standard deviation derived by variance composition
sigma_sigma_locmean = 0.5 * std( ma.lfit.b(good_enough_fit) ) / sigma_locmean;

fprintf('Localization error estimated to be s = %.3e ± %.3e (mean ± std, N = %d),\n', ...
    sigma_locmean, sigma_sigma_locmean, sum(good_enough_fit));
fprintf('to compare to the simulated value: %.3e.\n', BAD_XY_TYPICAL_OFFSET);

close all
clear all

% Number of particles to simulate
N_PARTICLES = 100;

% Probability to miss a detection.
P_GAPS = 0.3;

% Probability that a XY position is incorrect.
P_BAD_XY = 0.2;

% Typical bad localization error.
BAD_XY_TYPICAL_OFFSET = 0.2; % µm

% Probability that a track actually follows two particles.
P_BAD_TRACK = 0.2;

% Typical distance between the 2 particles erroneously tracked together
BAD_TRACK_DISTANCE = 0.5; %

N_TIME_STEPS = 500;
N_DIM = 2; % 2D

kT = 4.2821e-21; % 37ºC
D  = 1e-3; % µm^2/s
dT = 0.2; % s

SIZE = 10; % µm

tracks = cell(N_PARTICLES, 1);

k = sqrt(N_DIM * D * dT);

n_missed = 0;
n_bad_xy = 0;
n_bad_track = 0;

for i_spot = 1 : N_PARTICLES

    % Time
    time_steps = max(1, round(N_TIME_STEPS + N_TIME_STEPS/4*randn));
    time = (0 : time_steps-1)' * dT + dT * floor(N_TIME_STEPS / 4  * rand);


    % Initial position
    X0 = SIZE .* rand(1, N_DIM);

    % Integrate uncorrelated displacement
    dX = k * randn(time_steps, N_DIM);
    dX(1, :) = X0;
    X = cumsum(dX, 1);

    % First deal with incorrect detection

    incorrect_detection = rand(time_steps, 1) < P_BAD_XY;
    n_incorrect_detection = sum(incorrect_detection);
    bad_dx = BAD_XY_TYPICAL_OFFSET * randn(n_incorrect_detection, N_DIM);
    X(incorrect_detection, :) = X(incorrect_detection, :) + bad_dx;

    n_bad_xy = n_bad_xy + n_incorrect_detection;


    % Deal with two particle confused as one track.

    bad_track = rand < P_BAD_TRACK;
    if bad_track
        % It is a bad track. So at a random time, all the X coordinates
        % will actually follow another particle, which is off by a certain
        % distance:
        switch_time = 1 + floor(rand * (time_steps-1));
        dx_other_particle = BAD_TRACK_DISTANCE * randn(1, N_DIM);
        dx_other_particle = repmat(dx_other_particle, [(time_steps-switch_time+1) 1]);
        X(switch_time:end, :) = X(switch_time:end, :) + dx_other_particle;

        n_bad_track = n_bad_track + 1;
    end

    % Deal with missing frames

    missing_frames = rand(time_steps, 1) < P_GAPS;
    X(missing_frames, :) = [];
    time(missing_frames) = [];

    n_missed = n_missed + sum(missing_frames);


    % Store
    tracks{i_spot} = [ time X];

end

fprintf('Generated %d tracks, with:\n', N_PARTICLES)
fprintf(' - %d missed detections\n', n_missed)
fprintf(' - %d bad detections\n', n_bad_xy)
fprintf(' - %d bad tracks\n', n_bad_track)

ma = msdanalyzer(2, 'µm', 's');
ma = ma.addAll(tracks);
ma.plotTracks;
ma.labelPlotTracks;

ma = ma.computeMSD;
figure
ma.plotMeanMSD(gca, true)
ma = ma.fitMSD;

good_enough_fit = ma.lfit.r2fit > 0.8;
Dmean = mean( ma.lfit.a(good_enough_fit) ) / 2 / ma.n_dim;
Dstd  =  std( ma.lfit.a(good_enough_fit) ) / 2 / ma.n_dim;

fprintf('Found D = %.3e ± %.3e (mean ± std, N = %d)\n', ...
    Dmean, Dstd, sum(good_enough_fit));
