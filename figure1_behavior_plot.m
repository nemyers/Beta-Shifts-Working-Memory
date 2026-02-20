% figure1_behavior_plot.m - Plot behavioral switch costs and psychometric curves
%
% This script loads pre-calculated behavioral data and generates Figure 1
% from the paper, showing:
%   - Psychometric curves with swap error modeling
%   - Precision, guess rate, and swap rate by trial type
%   - Regression coefficients for item distance and congruency
%   - Item distance effects on accuracy and RT
%
% KEY FINDING:
%   Switch costs in accuracy and RT scale with item distance.
%   Swap errors are more frequent on switch trials.
%
% INPUTS:
%   - results/Study12_TFPow_Beh_Decoding.mat
%     (contains behavioral data for all subjects)
%
% OUTPUTS:
%   - Figure 1 panels (saved as PNG, EPS, and PDF)
%   - Statistical results printed to console
%
% Author: Nicholas E. Myers
% Cleaned for GitHub: 2026
% For: Myers, Stokes, & Muhle-Karbe (2026) J. Neuroscience

clearvars;
close all;
clc;

%% 1. SETUP AND LOAD DATA

% Load configuration
paths = switchbeta_config();

% Subject information
sublist = [2 4 5 7 11 14 15 16 17 18 19 21 22 23 27 29 31 32 35 36 2:31];
explist = [ones(1,20)*1 ones(1,30)*2];
nsubs   = length(sublist);

% Subject pairs (participated in both experiments)
subpairs = [05 10 11 12 13 14 15; 
            28 23 22 24 34 42 41];

% Load pre-calculated behavioral data
fprintf('Loading behavioral data...\n');
data_file = fullfile(paths.results, 'Study12_TFPow_Beh_Decoding.mat');
if ~exist(data_file, 'file')
    error('Behavioral data file not found: %s\nRun figure1_behavior_calculate.m first.', data_file);
end
load(data_file);
nsub = length(tfdata);

%% 2. EXTRACT AND ORGANIZE BEHAVIORAL DATA

fprintf('Processing behavioral data...\n');

% Initialize output arrays
outdata = [];
dvdata = [];
repeffect = [];

% Progress indicator
fprintf('Analyzing: ');
backspace_string = '';

for isub = 1:nsub
    % Progress
    progress_msg = sprintf('S%02d/%02d', isub, nsub);
    fprintf([backspace_string progress_msg]);
    backspace_string = repmat('\b', [1 length(progress_msg)]);
    
    % Extract behavioral data for this subject
    beh = tfdata(isub).behaviour;
    experiment = tfdata(isub).experiment;
    explist(isub, 1) = experiment;
    
    % Response, RT, accuracy
    rsp = 2 - beh.resp;  % flip so positive DV = response 1
    rt  = beh.rt;
    ac  = beh.acc;
    
    % Decision variables (probe offset from cued and uncued items)
    dv = round(beh.dv * 1e5) / 1e5;
    
    % Repetition number (1=switch, 2=first repeat, 3=second+ repeat)
    rep = beh.repnum;
    maxrep = 3;
    rep(rep > maxrep) = maxrep;
    
    % Trial info
    bt = beh.blocktrial;
    bl = beh.block;
    
    % Template distance (angular distance between two memory items)
    td = circ_dist(beh.template_angles(:,1)*2*pi/180, ...
                   beh.template_angles(:,2)*2*pi/180) / pi * 90;
    td = round(td * 1e4) / 1e4;
    td = abs(td);
    
    utd = unique(td);
    ntd = length(utd);
    
    % Wrap probe offsets to -45 to +45 degrees
    d = dv(:, 1);
    idx_wrap = abs(d) > 45;
    d(idx_wrap) = (90 - abs(d(idx_wrap))) .* sign(d(idx_wrap));
    dwrap = d;
    d = abs(d);
    
    d2 = dv(:, 2);
    idx_wrap = abs(d2) > 45;
    d2(idx_wrap) = (90 - abs(d2(idx_wrap))) .* sign(d2(idx_wrap));
    dwrap(:, 2) = d2;
    d2 = abs(d2);
    
    udv = unique(d);
    ndv = length(udv);
    
    % Congruent responses (same response for both items)
    congresp = diff(sign(beh.dv), [], 2) == 0;
    
    % Minimum block trial (exclude first trial in block)
    blocktrialmin = 2;
    
    %% 2.1 CALCULATE SUMMARY STATISTICS
    
    % Accuracy and RT by trial type
    for irep = 1:maxrep
        itrl = bt >= blocktrialmin & rep == irep;
        repeffect(isub, irep, 1) = mean(ac(itrl));
        itrl = itrl & ac == 1;
        repeffect(isub, irep, 2) = median(rt(itrl));
    end
    
    % Organize data by template distance, trial type, congruency
    nd2 = length(unique(dwrap(:, 1)));
    udv2 = unique(dwrap(:, 1));
    udv2_abs = unique(abs(dwrap(:, 1)));
    nd2_abs = length(udv2_abs);
    
    for irep = 1:3
        for idv = 1:2  % cued vs. uncued item
            for ind2 = 1:nd2
                itrl = bt >= blocktrialmin & rep == irep & dwrap(:, idv) == udv2(ind2);
                if nnz(itrl) > 0
                    misbnd_resp(isub, ind2, irep, idv) = mean(rsp(itrl));
                else
                    misbnd_resp(isub, ind2, irep, idv) = NaN;
                end
                itrl = itrl & ac == 1;
                if nnz(itrl) > 0
                    misbnd_rt(isub, ind2, irep, idv) = median(rt(itrl));
                else
                    misbnd_rt(isub, ind2, irep, idv) = NaN;
                end
            end
            
            for ind2_abs = 1:nd2_abs
                itrl = bt >= blocktrialmin & rep == irep & abs(dwrap(:, idv)) == udv2_abs(ind2_abs);
                if nnz(itrl) > 0
                    misbnd_ac(isub, ind2_abs, irep, idv) = mean(ac(itrl));
                else
                    misbnd_ac(isub, ind2_abs, irep, idv) = NaN;
                end
            end
        end
    end
    
    % Organize by template distance and trial type
    for itd = 1:ntd
        for irep = 1:2
            for icong = 1:2
                for id = 1:ndv
                    itrl = bt >= blocktrialmin & rep == irep & td == utd(itd) & ...
                           congresp == (2-icong) & d == udv(id);
                    outdata(isub, itd, irep, 1, icong, id) = mean(ac(itrl));
                    itrl = itrl & ac == 1;
                    outdata(isub, itd, irep, 2, icong, id) = median(rt(itrl));
                end
            end
        end
    end
    
    %% 2.2 FIT PSYCHOMETRIC MODELS
    
    % Fit models for each trial type
    itrl = bt >= blocktrialmin;
    dmin = -45;
    dmax = +45;
    minrt = 0.1;
    maxrt = 2;
    
    for it = 1:3
        itrl1 = rep == it & d > dmin & d < dmax & itrl;
        itrl2 = itrl1 & (ac == 1) & rt > minrt & rt < maxrt;
        
        tac(isub, it, 1) = nanmean(ac(itrl1));
        trt(isub, it, 1) = nanmedian(rt(itrl2));
        
        % 2-parameter model: slope (precision) and lapse rate
        design = dwrap(itrl1, 1);
        data = (rsp(itrl1) > 0);
        
        start_params = [10.0, 0.05];  % [sigma, lapse_rate]
        lower_bounds = [0, 0];
        upper_bounds = [50, 1];
        
        options = optimset('fmincon');
        options = optimset(options, 'Algorithm', 'interior-point', 'Display', 'off');
        
        [params, fval] = fmincon(@(P) costaccfun(P, design, data, false), ...
                                 start_params, [], [], [], [], ...
                                 lower_bounds, upper_bounds, [], options);
        mle_vals(isub, it, 1, :) = [params, fval];
        
        % 3-parameter model: add swap rate
        design = dwrap(itrl1, :);
        start_params = [10.0, 0.05, 0.05];  % [sigma, lapse_rate, swap_rate]
        lower_bounds = [0, 0, 0];
        upper_bounds = [50, 1, 1];
        
        [params, fval] = fmincon(@(P) costaccfun_misbind(P, design, data, false), ...
                                 start_params, [], [], [], [], ...
                                 lower_bounds, upper_bounds, [], options);
        mle_vals_misbind(isub, it, 1, :) = [params, fval];
        
        % Predict responses with fitted model
        respfun = @(P, x) ((1 - P(2)/2 - P(3)) .* normcdf(x(:,1), 0, P(1)) + ...
                           P(3) .* normcdf(x(:,2), 0, P(1)) + P(2)/2);
        yhat = respfun(params, design);
        
        for idv = 1:2
            for ind2 = 1:nd2
                misbnd_resp_hat(isub, ind2, it, idv) = mean(yhat(design(:, idv) == udv2(ind2)));
            end
        end
        
        %% 2.3 FIT REGRESSION MODELS
        
        % Logistic regression for accuracy
        y = rsp(itrl1) > 0;
        X = dwrap(itrl1, :);
        betas = glmfit(X, y, 'binomial', 'link', 'probit');
        betamat_misbnd(isub, it, :) = betas;
        
        % Logistic regression with multiple predictors
        y = ac(itrl1) > 0;
        X = [td(itrl1), congresp(itrl1), d(itrl1), d2(itrl1), bt(itrl1), bl(itrl1)];
        betas = glmfit(X, y, 'binomial', 'link', 'probit');
        betamat(isub, it, :) = betas;
        
        % Linear regression for RT
        y = rt(itrl1);
        X = [td(itrl1), congresp(itrl1), d(itrl1), d2(itrl1), bt(itrl1), bl(itrl1)];
        betas = glmfit(X, y, 'normal', 'link', 'identity');
        betamat_rt(isub, it, :) = betas;
    end
end

fprintf(' - done!\n');

%% 3. MERGE DATA FROM SUBJECTS IN BOTH EXPERIMENTS

fprintf('Merging data from subjects in both experiments...\n');

misbnd_acs = Switch_MergeExperiments(misbnd_ac, subpairs, explist, false, +Inf);
misbnd_resps = Switch_MergeExperiments(misbnd_resp, subpairs, explist, false, +Inf);
misbnd_rts = Switch_MergeExperiments(misbnd_rt, subpairs, explist, false, +Inf);
bmat_misbnd = Switch_MergeExperiments(betamat_misbnd, subpairs, explist, false, +Inf);
bmat_rt = Switch_MergeExperiments(betamat_rt, subpairs, explist, false, +Inf);
bmat = Switch_MergeExperiments(betamat, subpairs, explist, false, +Inf);
tacs = Switch_MergeExperiments(tac, subpairs, explist, false, +Inf);
trts = Switch_MergeExperiments(trt, subpairs, explist, false, +Inf);
repfx = Switch_MergeExperiments(repeffect, subpairs, explist, false, +Inf);
outdat = Switch_MergeExperiments(outdata, subpairs, explist, false, +Inf);
mles = Switch_MergeExperiments(mle_vals, subpairs, explist, false, +Inf);
mles_misbind = Switch_MergeExperiments(mle_vals_misbind, subpairs, explist, false, +Inf);
misbnd_resps_hat = Switch_MergeExperiments(misbnd_resp_hat, subpairs, explist, false, +Inf);

% Extract model parameters
mles = squeeze(mles(:, :, :, 1:2));
mles_misbind = squeeze(mles_misbind(:, :, :, 1:4));

%% 4. FIGURE: PSYCHOMETRIC CURVES WITH SWAP ERRORS

fprintf('Creating Figure: Psychometric curves...\n');

close all;
SetupPlotting();
figure('Position', [50, 50, 1500, 600]);

% Define colors
plotcol = [0.2000, 0.3500, 0.8000;   % blue (switch)
           0.3000, 0.3000, 0.3000;    % dark grey (repeat 1)
           0.7000, 0.7000, 0.7000];   % light grey (repeat 2+)

% Panel A & B: Psychometric curves relative to cued and uncued items
for ipanel = 1:2
    subplot(3, 9, [1 2 3 10 11 12 19 20 21] + 3*(ipanel-1));
    
    for icond = 1:3
        hold on;
        
        % Reference lines
        plot([-45 45], [0.5 0.5], 'k--', 'Color', ones(1,3)*0.75, 'LineWidth', 1);
        plot([0 0], [0 1], 'k--', 'Color', ones(1,3)*0.75, 'LineWidth', 1);
        
        % Plot model fit
        cfg = struct();
        cfg.alpha = 0.25;
        cfg.linewidth = 0.5;
        hp(icond) = plotpatch(misbnd_resps_hat(:, :, icond, ipanel), udv2', plotcol(icond, :), cfg);
        
        % Plot data points
        cfg = struct();
        cfg.dotsize = 25;
        plotDotsWithErrors(misbnd_resps(:, :, icond, ipanel), udv2', plotcol(icond, :), cfg);
    end
    
    % Format axes
    set(gca, 'XTick', round(udv2));
    set(gca, 'YTick', 0.0:0.1:1, 'YTickLabel', ...
        {'0', '10', '20', '30', '40', '50', '60', '70', '80', '90', '100'});
    
    if ipanel == 1
        legend(hp(1:3), {'Switch', 'Repeat1', 'Repeat2+'}, 'Location', 'NorthWest');
        legend boxoff;
        ylabel('% clockwise');
        title(sprintf('Resp relative to\nCued Item'));
    else
        title(sprintf('Resp relative to\nUncued Item'));
    end
    
    ylim([0 1]);
    xlabel('Probe Offset (°)');
end

% Panels C-E: Model parameters (precision, guess rate, swap rate)
varnames = {sprintf('%c', 963), 'Guess Rate', 'Swap Rate'};  % sigma symbol

for ivar = 1:3
    subplot(3, 9, [7 8 9] + 9*(ivar-1));
    
    varlabel = varnames{ivar};
    plotdat = mles_misbind(:, :, ivar);
    title(varlabel);
    
    if ivar == 1
        ylims = [0.0, 50.0];
    else
        ylims = [0.0, 1.0];
    end
    xlims = [0.5, 3.5];
    
    hold on;
    
    nx = size(plotdat, 1);
    nt = size(plotdat, 2);
    
    mubar = mean(plotdat, 1);
    sdbar = std(plotdat, [], 1) / sqrt(nx);
    
    % Plot individual subjects and group statistics
    for ibar = 1:nt
        xloc = ibar;
        boxchart(xloc*ones(nx, 1), plotdat(:, ibar), ...
                 'BoxFaceColor', plotcol(ibar, :), 'Notch', 'off', 'MarkerColor', [0 0 0]);
        plot(xloc + [-0.25 +0.25], ones(1, 2)*mubar(ibar), 'k-', 'LineWidth', 0.5);
        
        for isub = 1:nx
            scatter(xloc + (rand(1) - 0.5)*0.2, plotdat(isub, ibar), 25, ...
                   'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', plotcol(ibar, :), ...
                   'MarkerFaceAlpha', 0.5);
        end
    end
    
    offsetAxes();
    maxdat = max(abs(plotdat(:))) + range(plotdat(:))*0.05;
    ylabel(varlabel);
    xlabel('Trial Type');
    set(gca, 'XTick', 1:nt, 'XTickLabel', {'Switch', 'Repeat 1', 'Repeat 2+'});
    xlim(xlims);
    
    % Print statistics
    fprintf('\nEstimate for %s:\n', varlabel);
    fprintf('Switch:   %3.3f ± %3.3f (s.e.m.)\n', mubar(1), sdbar(1));
    fprintf('Repeat 1: %3.3f ± %3.3f (s.e.m.)\n', mubar(2), sdbar(2));
    fprintf('Repeat 2+:%3.3f ± %3.3f (s.e.m.)\n', mubar(3), sdbar(3));
    
    % Statistical tests
    [~, pval, ~, stat] = ttest(diff(plotdat(:, [1 2]), [], 2));
    d = computeCohen_d(plotdat(:, 1), plotdat(:, 2), 'paired');
    mysigstar(gca, [1.1 1.9], maxdat + range(plotdat(:))*0.1, pval, 'k', 0);
    fprintf('Switch vs. Rep1: t(%d)=%.4f, p=%g, Cohen''s d=%.2f\n', ...
            stat.df, stat.tstat, pval, d);
    
    [~, pval, ~, stat] = ttest(diff(plotdat(:, [2 3]), [], 2));
    d = computeCohen_d(plotdat(:, 2), plotdat(:, 3), 'paired');
    mysigstar(gca, [2.1 2.9], maxdat + range(plotdat(:))*0.1, pval, 'k', 0);
    fprintf('Rep1 vs Rep2+: t(%d)=%.4f, p=%g, Cohen''s d=%.2f\n', ...
            stat.df, stat.tstat, pval, d);
    
    [~, pval, ~, stat] = ttest(diff(plotdat(:, [1 3]), [], 2));
    d = computeCohen_d(plotdat(:, 1), plotdat(:, 3), 'paired');
    mysigstar(gca, [1.1 2.9], maxdat + range(plotdat(:))*0.15, pval, 'k', 0);
    fprintf('Switch vs Rep2+: t(%d)=%.4f, p=%g, Cohen''s d=%.2f\n\n', ...
            stat.df, stat.tstat, pval, d);
    
    ylim([0, maxdat + range(plotdat(:))*0.20]);
end

% Save figure
do_print = false;
if do_print
    figpath = fullfile(paths.results, 'figures', 'Figure1_Supp1_Misbinding');
    if ~exist(figpath, 'dir')
        mkdir(figpath);
    end
    fname = fullfile(figpath, sprintf('Psychometric-with-Misbinding_%s', ...
                     datestr(now, 'yyyymmdd-HHMM')));
    print(gcf, '-dpng', fname);
    print(gcf, '-depsc2', fname, '-painters');
    print(gcf, '-dpdf', fname, '-painters');
    fprintf('Saved: %s\n', fname);
end

%% 5. FIGURE: ITEM DISTANCE EFFECTS

fprintf('Creating Figure: Item distance effects...\n');

close all;
SetupPlotting();
figure('Position', [50, 50, 700, 300]);

% Panel A: Accuracy by item distance (congruent trials)
subplot(2, 2, 1);
congruent_index = 1;
hold on;
hpatch = [];

for icond = 1:2
    % Fit linear regression
    betas = permute(massGLM(utd, permute(outdat(:, :, icond, 1, congruent_index), [2 1 3 4])), [2 1]);
    
    % Plot regression line
    cfg = struct();
    cfg.linewidth = 0.5;
    cfg.alpha = 0.25;
    hpatch(icond) = plotpatch(betas(:, 1) + betas(:, 2)*utd', utd', plotcol(icond, :), cfg);
    
    % Plot data points
    cfg = struct();
    cfg.dotsize = 25;
    plotDotsWithErrors(outdat(:, :, icond, 1, congruent_index), utd', plotcol(icond, :), cfg);
end

legend(hpatch, {'Switch', 'Repeat'}, 'Location', 'SouthWest');
legend boxoff;
offsetAxes();
ylim([0.50, 0.90]);
set(gca, 'YTick', 0.50:0.05:0.90, ...
    'YTickLabel', {'50', '', '60', '', '70', '', '80', '', '90'});
set(gca, 'XTick', round(utd));
ylabel('% correct');
xlim([0, 90]);
title('Congruent');

% Statistical tests for accuracy
betas = permute(massGLM(utd, permute(outdat(:, :, :, 1, congruent_index), [2 1 3])), [2 1 3]);
mubeta = mean(squeeze(betas(:, 2, :)), 1);
sdbeta = std(squeeze(betas(:, 2, :)), [], 1) / sqrt(size(betas, 1));
[~, p, ~, s] = ttest(squeeze(betas(:, 2, :)));
d1 = computeCohen_d(betas(:, 2, 1), betas(:, 2, 1)*0, 'paired');
d2 = computeCohen_d(betas(:, 2, 2), betas(:, 2, 2)*0, 'paired');

fprintf('\n\nCONGRUENT TRIALS - Effect of template distance (Accuracy):\n');
fprintf('Switch: %.2f ± %.2f%% per 10 deg, t(%d)=%.2f, p=%g, d=%.2f\n', ...
        mubeta(1)*1000, sdbeta(1)*1000, s.df(1), s.tstat(1), p(1), d1);
fprintf('Repeat: %.2f ± %.2f%% per 10 deg, t(%d)=%.2f, p=%g, d=%.2f\n', ...
        mubeta(2)*1000, sdbeta(2)*1000, s.df(2), s.tstat(2), p(2), d2);

[~, p, ~, s] = ttest(-diff(betas(:, 2, :), [], 3));
d = computeCohen_d(betas(:, 2, 1), betas(:, 2, 2), 'paired');
fprintf('Switch vs Repeat: t(%d)=%.2f, p=%g, d=%.2f\n', s.df, s.tstat, p, d);

% Panel B: RT by item distance (congruent trials)
subplot(2, 2, 3);
hold on;

for icond = 1:2
    % Fit linear regression
    betas = permute(massGLM(utd, permute(outdat(:, :, icond, 2, congruent_index), [2 1 3 4])), [2 1]);
    
    % Plot regression line
    cfg = struct();
    cfg.linewidth = 0.5;
    cfg.alpha = 0.25;
    plotpatch(betas(:, 1) + betas(:, 2)*utd', utd', plotcol(icond, :), cfg);
    
    % Plot data points
    cfg = struct();
    cfg.dotsize = 25;
    plotDotsWithErrors(outdat(:, :, icond, 2, congruent_index), utd', plotcol(icond, :), cfg);
end

offsetAxes();
ylim([0.55, 0.95]);
set(gca, 'YTick', 0.6:0.1:0.9);
xlim([0, 90]);
set(gca, 'XTick', round(utd));
ylabel('RT (s)');
xlabel('Template Distance (°)');

% Statistical tests for RT
betas = permute(massGLM(utd, permute(outdat(:, :, :, 2, congruent_index), [2 1 3])), [2 1 3]);
mubeta = mean(squeeze(betas(:, 2, :)), 1);
sdbeta = std(squeeze(betas(:, 2, :)), [], 1) / sqrt(size(betas, 1));
[~, p, ~, s] = ttest(squeeze(betas(:, 2, :)));
d1 = computeCohen_d(betas(:, 2, 1), betas(:, 2, 1)*0, 'paired');
d2 = computeCohen_d(betas(:, 2, 2), betas(:, 2, 2)*0, 'paired');

fprintf('\nEffect of template distance (RT):\n');
fprintf('Switch: %.2f ± %.2f ms per 10 deg, t(%d)=%.2f, p=%g, d=%.2f\n', ...
        mubeta(1)*10000, sdbeta(1)*10000, s.df(1), s.tstat(1), p(1), d1);
fprintf('Repeat: %.2f ± %.2f ms per 10 deg, t(%d)=%.2f, p=%g, d=%.2f\n', ...
        mubeta(2)*10000, sdbeta(2)*10000, s.df(2), s.tstat(2), p(2), d2);

[~, p, ~, s] = ttest(-diff(betas(:, 2, :), [], 3));
d = computeCohen_d(betas(:, 2, 1), betas(:, 2, 2), 'paired');
fprintf('Switch vs Repeat slope: t(%d)=%.2f, p=%g, d=%.2f\n', s.df, s.tstat, p, d);

% Save figure
do_print = 0;
if do_print
    figpath = fullfile(paths.results, 'figures', 'Figure1');
    if ~exist(figpath, 'dir')
        mkdir(figpath);
    end
    fname = fullfile(figpath, sprintf('ItemDistance-Effects_%s', ...
                     datestr(now, 'yyyymmdd-HHMM')));
    print(gcf, '-dpng', fname);
    print(gcf, '-depsc2', fname, '-painters');
    print(gcf, '-dpdf', fname, '-painters');
    fprintf('Saved: %s\n', fname);
end

fprintf('\n=== Plotting Complete! ===\n');

%% HELPER FUNCTIONS

function cost = costaccfun(P, design, data, verbose)
    % COSTACCFUN - Cost function for 2-parameter psychometric fit
    % P = [sigma, lapse_rate]
    sigma = P(1);
    lapse = P(2);
    
    % Predicted probability of clockwise response
    pred = (1 - lapse) .* normcdf(design, 0, sigma) + lapse/2;
    
    % Negative log likelihood
    cost = -sum(log(pred(data == 1))) - sum(log(1 - pred(data == 0)));
    
    if verbose
        fprintf('sigma=%.2f, lapse=%.3f, cost=%.2f\n', sigma, lapse, cost);
    end
end

function cost = costaccfun_misbind(P, design, data, verbose)
    % COSTACCFUN_MISBIND - Cost function for 3-parameter psychometric fit with swap errors
    % P = [sigma, lapse_rate, swap_rate]
    sigma = P(1);
    lapse = P(2);
    swap = P(3);
    
    % Predicted probability accounting for both items
    pred = (1 - lapse/2 - swap) .* normcdf(design(:, 1), 0, sigma) + ...
           swap .* normcdf(design(:, 2), 0, sigma) + lapse/2;
    
    % Negative log likelihood
    cost = -sum(log(pred(data == 1))) - sum(log(1 - pred(data == 0)));
    
    if verbose
        fprintf('sigma=%.2f, lapse=%.3f, swap=%.3f, cost=%.2f\n', ...
                sigma, lapse, swap, cost);
    end
end
