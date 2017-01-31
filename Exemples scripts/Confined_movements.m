%% Confined movements.

clear all
close all
clc


SPACE_UNITS = 'µm';
TIME_UNITS = 's';

N_PARTICLES = 100;
N_TIME_STEPS = 200;
N_DIM = 2; % 2D
SIZE = 5; % µm

kT = 4.2821e-21; % kBoltzman x T @ 37ºC
D  = 1e-3; % µm^2/s
dT = 0.05; % s short frame interval - important for later

k = sqrt(2 * D * dT);

% Particle in a potential: settings the 'stiffness' of the energy potential
% Typical diameter of the trap (still in micron)
Ltrap = 0.05; % µm
Ktrap = kT / Ltrap^2; % = thermal energy / trap size ^ 2

tracks = cell(N_PARTICLES, 1);

for i_spot = 1 : N_PARTICLES

    % Time
    time = (0 : N_TIME_STEPS-1)' * dT;

    % Initial position
    X0 = SIZE .* rand(1, N_DIM);

    % Energy potential:
    V = @(x) 0.5 * Ktrap * sum (x .^ 2); % Unused, just to show
    Fx = @(x) - Ktrap * (x - X0); % Is a vector

    % Position
    X = zeros(N_TIME_STEPS, N_DIM);

    % Init first step
    X(1, :) = X0;

    % Iterate
    for j = 2 : N_TIME_STEPS

        dxtrap = D/kT * Fx(X(j-1,:)) * dT; % ad hoc displacement
        dxbrownian = k * randn(1, N_DIM);

        X(j,:) = X(j-1,:) + dxtrap + dxbrownian;

    end

    % Store
    tracks{i_spot} = [ time X];

end

ma = msdanalyzer(N_DIM, 'µm', 's');
ma = ma.addAll(tracks);

% Plot trajectories
[hps, ha] = ma.plotTracks;
ma.labelPlotTracks(ha);

ma = ma.computeMSD;
figure
hmsd = ma.plotMeanMSD(gca, true);

fprintf(['Time threshold for confined motion: %.1f ' TIME_UNITS '.\n'], ...
    kT / (2*Ktrap*D) )

[fo, gof] = ma.fitMeanMSD( 0.1 );
plot(fo)
legend off
ma.labelPlotMSD

ma = ma.fitLogLogMSD(0.5);
ma.loglogfit

mean(ma.loglogfit.alpha)

r2fits = ma.loglogfit.r2fit;
alphas = ma.loglogfit.alpha;

R2LIMIT = 0.8;

% Remove bad fits
bad_fits = r2fits < R2LIMIT;
fprintf('Keeping %d fits (R2 > %.2f).\n', sum(~bad_fits), R2LIMIT);
alphas(bad_fits) = [];

% T-test
[htest, pval] = ttest(alphas, 1, 0.05, 'left');

if ~htest
    [htest, pval] = ttest(alphas, 1, 0.05);
end

% Prepare string
str = { [ '\alpha = ' sprintf('%.2f ± %.2f (mean ± std, N = %d)', mean(alphas), std(alphas), numel(alphas)) ] };

if htest
    str{2} = sprintf('Significantly below 1, with p = %.2g', pval);
else
    str{2} = sprintf('Not significantly differend from 1, with p = %.2g', pval);
end

figure
hist(alphas);
box off
xlabel('\alpha')
ylabel('#')

yl = ylim(gca);
xl = xlim(gca);
text(xl(2), yl(2)+2, str, ...
    'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'top', ...
    'FontSize', 16)
title('\alpha values distribution', ...
    'FontSize', 20)
ylim([0 yl(2)+2])

gammas = ma.loglogfit.gamma;
gammas(bad_fits) = []; % discard bad fits, like for alpha

Dmean = mean( gammas ) / 2 / ma.n_dim;
Dstd  =  std( gammas ) / 2 / ma.n_dim;

fprintf('Estimation of the diffusion coefficient from log-log fit of the MSD curves:\n')
fprintf('D = %.2e ± %.2e (mean ± std, N = %d)\n', ...
    Dmean, Dstd, numel(gammas));
