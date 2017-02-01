%% Velocities and velocity auto-correlation.

close all
clear all

SPACE_UNITS = 'µm';
TIME_UNITS = 's';
SIZE = 2; % µm
N_PARTICLES = 100;
N_TIME_STEPS = 100;
N_DIM = 2; % 2D
D  = 1e-3; % µm^2/s
dT = 0.05; % s
k = sqrt(2 * D * dT);

tracks = cell(N_PARTICLES, 1);
for i = 1 : N_PARTICLES

    time = (0 : N_TIME_STEPS-1)' * dT;
    X0 = SIZE .* rand(1, N_DIM);

    % Integrate uncorrelated displacement
    dX = k * randn(N_TIME_STEPS, N_DIM);
    dX(1, :) = X0;
    X = cumsum(dX, 1);

    % Store
    tracks{i} = [time X];

end
clear i X dX time X0

ma = msdanalyzer(2, SPACE_UNITS, TIME_UNITS);
ma = ma.addAll(tracks);

v = ma.getVelocities %#ok<NOPTS>

V = vertcat( v{:} );

figure

hist(V(:, 2:end), 50) % we don't want to include the time in the histogram
box off
xlabel([ 'Velocity (' SPACE_UNITS '/' TIME_UNITS ')' ])
ylabel('#')

mean(V(:,2:end))

std(V(:,2:end))

k/dT %#ok<NOPTS>

% we use a syntax close to the one we used for MSD curves
ma = ma.computeVCorr;
ma.vcorr

figure

ma.plotMeanVCorr