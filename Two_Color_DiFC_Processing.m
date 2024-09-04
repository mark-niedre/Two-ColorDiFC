function [out_dat] = Two_Color_DiFC_Processing(saveDir, rel_thresh_color1, rel_thresh_color2, limLast60Min, saveOrPlot,varargin)
% Two-Color DiFC Processing Code
% Processes Two-Color DiFC data
% saveDir: string of the directory name e.g. 'C:Files\Foo\FileName1\'
% rel_thresh_color1: the relative threshold for color 1 (GFP), usually 5
% rel_thresh_color2: the relative threshold for color 2 (tdTomato), usually 5
% limLast60Min: input 1 if you want to trim the data down to the last 60
%   minutes of data, else input 0 to process the entire set of data
% saveOrPlot: input 1 if you want the processed data plotted, input 2 if
%   you want to save the data to a folder and input the folder as an extra
%   input, input 0 if you don't want either of these things
%
% Make sure to change line 16 for your proccodes directory

%% Proccodes library (where processing functions are stored)
proccodes_library = "D:\Lenovo\OneDrive - Northeastern University\Niedre_Lab\DiFC Files\codes\proccodes";
addpath(genpath(proccodes_library));

%% Loading data, processing functions
%Getting file name
if ispc
    slash = '\';
else
    slash = '/';
end
index = find(saveDir == slash,1,'last');
stem = saveDir(index+1:end);
fname = saveDir;
    
if saveOrPlot == 2
    plot_flag = 1;
else
    plot_flag = 0;
end

% Set pre - processing parameters
color1_name = 'GFP';
color2_name = 'tdTomato';
background_sliding_window = 2.5;                  % length (in seconds) of median filter used for background subtraction
smoothing_siding_window = 0.001;                  % length (in seconds) of mean(?) filter used for data\ smoothing
snr_sliding_window = 1;                           % length (in seconds) of window for sliding estimate of background standard deviation
mov_thresh_window = 60;                           % length (in seconds) of segments used for moving threshold analysis
noise_alpha = 0.5;                                % smoothing factor for standard-deviation estimation (exponential smoothing)
std_thresh = 1.05;                                % multiplicative factor for excluding the current std:
prominence_factor = 1;                            % minimum peak prominence (versus adjacent peaks) is the threshold/prominence_factor
probe_distance = 3;                               % Distance between probes 1 and 2 (mm)
coinc_window = 0.06;                              % Peak coincident window
max_speed = 300;                                  % Maximum speed (mm/s) permitted for peak matching
max_score = 4;

%-----------------------------------------------------------------------------------------------------------------%
%% Running proccessing code
fprintf('Running Red_preProc...\n\n')
fprintf('Pre-processing %s...\n',stem)

% Load data from both probes
load([fname '_F1.mat'], 'time', 'data', 'params')
data_probe1 = data;
params1 = params;
load([fname '_F2.mat'], 'data', 'params')
data_probe2 = data;
params2 = params;

% Format data from both probes
data_color1 = [data_probe1(:,1) data_probe2(:,1)];
data_color2 = [data_probe1(:,2) data_probe2(:,2)];

if limLast60Min == 1
    data_color1 = data_color1(end+1-60*60*2000:end,:);
    data_color2 = data_color2(end+1-60*60*2000:end,:);
    time = time(1:size(data_color1,1));
end

params_color1 = [params1(1) params2(1)];
params_color1(1).name = [stem ' Probe 1 ' color1_name];
params_color1(2).name = [stem ' Probe 2 ' color1_name];
params_color1(1).color = [0 1 0];
params_color1(2).color = [0 1 0];
params_color2 = [params1(1) params2(1)];
params_color2(1).name = [stem ' Probe 1 ' color2_name];
params_color2(2).name = [stem ' Probe 2 ' color2_name];
params_color2(1).color = [1 0 0];
params_color2(2).color = [1 0 0];
clear data_probe1 data_probe2 params1 params2 data

% Pre-proc background subtracts data, calculates the noise/ moving peak
% threshhold, and identifies peak candidates. Processes all data in the
% data array
%%%%%%%%%%%%%%%%%% Color 1 preprocessing and coincident peak removal
[data_bs_color1, noise_color1, peaks_color1, thresh_curve_color1, ~] = preProc(data_color1, time, 'RelativeThresh', rel_thresh_color1,...
    'StdThresh', std_thresh, 'BackgroundWindow', background_sliding_window, 'SmoothWindow', smoothing_siding_window,...
    'RemoveBunchedPeaks', false, 'CoincidenceWindow', coinc_window, 'ProminenceFactor', prominence_factor,...
    'NoiseAlpha', noise_alpha, 'DynamicWindow', mov_thresh_window, 'StdWindow', snr_sliding_window);

%%%%%%%%%%%%%%%%%% Color 2 preprocessing and coincident peak removal
[data_bs_color2, noise_color2, peaks_color2, thresh_curve_color2, ~] = preProc(data_color2, time, 'RelativeThresh', rel_thresh_color2,...
    'StdThresh', std_thresh, 'BackgroundWindow', background_sliding_window, 'SmoothWindow', smoothing_siding_window,...
    'RemoveBunchedPeaks', false, 'CoincidenceWindow', coinc_window, 'ProminenceFactor', prominence_factor,...
    'NoiseAlpha', noise_alpha, 'DynamicWindow', mov_thresh_window, 'StdWindow', snr_sliding_window);

all_peaks(1) = peaks_color1(1);
all_peaks(2) = peaks_color2(1);
all_peaks(3) = peaks_color1(2);
all_peaks(4) = peaks_color2(2);

% Remove coincident peaks that occur on all channels (aka motion artifact)
disp('Removing 4-channel Coincident Peaks');
[all_peaks(1), all_peaks(2), all_peaks(3), all_peaks(4),...
    motion(1), motion(2), motion(3), motion(4), motion_count] = ...
    remove4ChannelCoincPeaks(all_peaks(1), all_peaks(2), all_peaks(3), all_peaks(4), time, 'CoincidenceWindow', coinc_window);

%-----------------------------------------------------------------------------------------------------------------%
% MATCH PEAKS IN FORWARD DIRECTION
%-----------------------------------------------------------------------------------------------------------------%
disp('Matching peaks in the forward direction')
%%%%%%%%%%%%%%%%%% Color 1 forward matches
[fwd_peaks_color1(1), fwd_peaks_color1(2), fwd_speed_color1, fwd_score_color1] =...
    matchDirectionalPeaks(all_peaks(1), all_peaks(3), time, probe_distance, 'MaxSpeed', max_speed, 'MaxScore', max_score);
%%%%%%%%%%%%%%%%%% Color 2 forward matches
[fwd_peaks_color2(1), fwd_peaks_color2(2), fwd_speed_color2, fwd_score_color2] =...
    matchDirectionalPeaks(all_peaks(2), all_peaks(4), time, probe_distance, 'MaxSpeed', max_speed, 'MaxScore', max_score);
%%%%%%%%%%%%%%%%%% Color 1 reverse matches
[rev_peaks_color1(1), rev_peaks_color1(2), rev_speed_color1, rev_score_color1] =...
    matchDirectionalPeaks(all_peaks(1), all_peaks(3), time, probe_distance, 'MaxSpeed', max_speed, 'MaxScore', max_score,'Direction','rev');
%%%%%%%%%%%%%%%%%% Color 2 reverse matches
[rev_peaks_color2(1), rev_peaks_color2(2), rev_speed_color2, rev_score_color2] =...
    matchDirectionalPeaks(all_peaks(2), all_peaks(4), time, probe_distance, 'MaxSpeed', max_speed, 'MaxScore', max_score,'Direction','rev');

%%%%%%%%%%%%%%%%%% Color 1 removing double matches
%Looking at potential double matches in Probe 1
[fwd_peaks_color1, rev_peaks_color1, fwd_score_color1, rev_score_color1, ...
    fwd_speed_color1,rev_speed_color1, doubleMatchCountsP1_color1] = ...
    determineBestDirectionalMatch(fwd_peaks_color1, rev_peaks_color1, fwd_score_color1, rev_score_color1,fwd_speed_color1,rev_speed_color1,1);
%Looking at potential double matches in Probe 2
[fwd_peaks_color1, rev_peaks_color1, fwd_score_color1, rev_score_color1, ...
    fwd_speed_color1,rev_speed_color1, doubleMatchCountsP2_color1] = ...
    determineBestDirectionalMatch(fwd_peaks_color1, rev_peaks_color1, fwd_score_color1, rev_score_color1,fwd_speed_color1,rev_speed_color1,2);

%%%%%%%%%%%%%%%%%% Color 2 removing double matches
%Looking at potential double matches in Probe 1
[fwd_peaks_color2, rev_peaks_color2, fwd_score_color2, rev_score_color2, ...
    fwd_speed_color2,rev_speed_color2, doubleMatchCountsP1_color2] = ...
    determineBestDirectionalMatch(fwd_peaks_color2, rev_peaks_color2, fwd_score_color2, rev_score_color2, fwd_speed_color2, rev_speed_color2,1);
%Looking at potential double matches in Probe 2
[fwd_peaks_color2, rev_peaks_color2, fwd_score_color2, rev_score_color2, ...
    fwd_speed_color2,rev_speed_color2, doubleMatchCountsP2_color2] = ...
    determineBestDirectionalMatch(fwd_peaks_color2, rev_peaks_color2, fwd_score_color2, rev_score_color2, fwd_speed_color2, rev_speed_color2,2);

%-----------------------------------------------------------------------------------------------------------------%
% FIND POTENTIAL CLUSTERS (Coincident peaks in both colors in each probe, PMTS 1/2 and PMTS 3/4)
%-----------------------------------------------------------------------------------------------------------------%

params_p1 = [params_color1(1) params_color2(1)];
params_p1(1).name = [stem ' Probe 1 ' color1_name];
params_p1(2).name = [stem ' Probe 1 ' color2_name];
params_p2 = [params_color1(2) params_color2(2)];
params_p2(1).name = [stem ' Probe 2 ' color1_name];
params_p2(2).name = [stem ' Probe 2 ' color2_name];

disp('Matching Coincident Peaks');
[TwoLambda_peaks_probe1(1), TwoLambda_peaks_probe1(2), ~] = match2ColorPeaks(all_peaks(1), all_peaks(2), time, 'CoincidenceWindow', coinc_window);
[TwoLambda_peaks_probe2(1), TwoLambda_peaks_probe2(2), ~] = match2ColorPeaks(all_peaks(3), all_peaks(4), time, 'CoincidenceWindow', coinc_window);

%-----------------------------------------------------------------------------------------------------------------%
% IDENTIFY COINCIDENT PEAKS THAT ARE HIGH OR LOW RATIO (AKA 2FP or 1FP)
%-----------------------------------------------------------------------------------------------------------------%

% Probe 1
if TwoLambda_peaks_probe1(1).count > 0
    [probe1_2FP, probe1_1FP] = TwoFP_vs_OneFP_Peaks(TwoLambda_peaks_probe1);
else
    probe1_2FP = struct; probe1_2FP(1).pks = []; probe1_2FP(1).locs = [];
    probe1_2FP(1).widths = []; probe1_2FP(1).proms = []; probe1_2FP(1).count = 0;
    probe1_2FP(2) = probe1_2FP(1); probe1_1FP = probe1_2FP;
end
% Probe 2
if TwoLambda_peaks_probe2(1).count > 0
    [probe2_2FP, probe2_1FP] = TwoFP_vs_OneFP_Peaks(TwoLambda_peaks_probe2);
else
    probe2_2FP = struct; probe2_2FP(1).pks = []; probe2_2FP(1).locs = [];
    probe2_2FP(1).widths = []; probe2_2FP(1).proms = []; probe2_2FP(1).count = 0;
    probe2_2FP(2) = probe2_2FP(1); probe2_1FP = probe2_2FP;
end

% For ease later, 2lambda peak amplitudes in matrix form:
OneFP_ampl = [probe1_1FP(1).pks probe1_1FP(2).pks;...
        probe2_1FP(1).pks probe2_1FP(2).pks];
[~, inds] = sort(OneFP_ampl, 2);
if ~isempty(inds) 
    OneFP_color = inds(:,2); 
else
    OneFP_color = [];
end

TwoFP_ampl = [probe1_2FP(1).pks probe1_2FP(2).pks;...
    probe2_2FP(1).pks probe2_2FP(2).pks];
[~, inds] = sort(TwoFP_ampl, 2);
if ~isempty(inds) 
    TwoFP_color = inds(:,2); 
else
    TwoFP_color = [];
end

%-----------------------------------------------------------------------------------------------------------------%
% PLOT FWD PEAKS AND COINCIDENT PEAKS (SINGLE COLOR DOM AND CLUSTER CANDIDATES
%-----------------------------------------------------------------------------------------------------------------%

if plot_flag == 1
    figure;
    t = tiledlayout(4,1);
    nexttile();
    plot_4channels_clusters(data_bs_color1(:,1), time, all_peaks(1), probe1_2FP(1), probe1_1FP(1), fwd_peaks_color1(1), rev_peaks_color1(1), motion(1), [0 1 0], thresh_curve_color1(:,1))
    title('Channel 1 (GFP Probe 1)');
    nexttile();
    plot_4channels_clusters(data_bs_color2(:,1), time, all_peaks(2), probe1_2FP(2), probe1_1FP(2), fwd_peaks_color2(1), rev_peaks_color2(1), motion(2), [1 0 0], thresh_curve_color2(:,1))
    title('Channel 2 (tdT Probe 1)');
    nexttile();
    plot_4channels_clusters(data_bs_color1(:,2), time, all_peaks(3), probe2_2FP(1), probe2_1FP(1), fwd_peaks_color1(2), rev_peaks_color1(2), motion(3), [0 1 0], thresh_curve_color1(:,2))
    title('Channel 3 (GFP Probe 2)');
    nexttile();
    plot_4channels_clusters(data_bs_color2(:,2), time, all_peaks(4), probe2_2FP(2), probe2_1FP(2), fwd_peaks_color2(2), rev_peaks_color2(2), motion(4), [1 0 0], thresh_curve_color2(:,2))
    title('Channel 4 (tdT Probe 2)');
    title(t,stem, 'Interpreter', 'none')
    xlabel(t,'Time (s)')
    ylabel(t,['Signal (' params(1).units ')'])
    t.Padding = 'compact';
end

%-----------------------------------------------------------------------------------------------------------------%
% Display useful results
%-----------------------------------------------------------------------------------------------------------------%

fprintf('\n--------------------------GFP---------------------------\n');
fprintf('Identified %g and %g peak candidates in GFP channels\n', peaks_color1(1).count, peaks_color1(2).count);
fprintf('Identified %g forward GFP matches\n', fwd_peaks_color1(1).count);
fprintf('Identified %g reverse GFP matches\n', rev_peaks_color1(1).count);
fprintf('Average relative threshold for GFP: %.3f nA\t %.3f nA\n', mean(thresh_curve_color1(:,1)), mean(thresh_curve_color1(:,2)));
fprintf('Average peak matched peak amplitude for GFP: %.3f nA\t %.3f nA\n', mean(fwd_peaks_color1(1).pks), mean(fwd_peaks_color1(2).pks));
fprintf('Calculated noise for GFP: %.3f nA\t %.3f nA\n', noise_color1(1), noise_color1(2));
fprintf('GFP SNR: %.3f dB\t %.3f dB\n', mean(20*log10(peaks_color1(1).pks./noise_color1(1))), mean(20*log10(peaks_color1(2).pks./noise_color1(2))));
fprintf('GFP matched Peaks SNR: %.3f dB\n', mean(20*log10(fwd_peaks_color1(1).pks./mean(noise_color1))));

fprintf('\n------------------------tdTomato------------------------\n');
fprintf('Identified %g and %g peak candidates in tdTomato channels\n', peaks_color2(1).count, peaks_color2(2).count);
fprintf('Identified %g forward tdTomato matches\n', fwd_peaks_color2(1).count);
fprintf('Identified %g reverse tdTomato matches\n', rev_peaks_color2(1).count);
fprintf('Average relative threshold for tdTomato: %.3f nA\t %.3f nA\n', mean(thresh_curve_color2(:,1)), mean(thresh_curve_color2(:,2)));
fprintf('Average peak matched peak amplitude for tdTomato: %.3f nA\t %.3f nA\n', mean(fwd_peaks_color2(1).pks), mean(fwd_peaks_color2(2).pks));
fprintf('Calculated noise for GFP: %.3f nA\t %.3f nA\n', noise_color2(1), noise_color2(2));
fprintf('tdTomato SNR: %.3f dB\t %.3f dB\n', mean(20*log10(peaks_color2(1).pks./noise_color2(1))), mean(20*log10(peaks_color2(2).pks./noise_color2(2))));
fprintf('tdTomato matched Peaks SNR: %.3f dB\n', mean(20*log10(fwd_peaks_color2(1).pks./mean(noise_color2))));
fprintf('\n\n');

out_dat.scan_length_samples = length(time);
out_dat.scan_length_minutes = time(end)/60;
out_dat.fwd_peaks_color1 = fwd_peaks_color1;
out_dat.fwd_peaks_color2 = fwd_peaks_color2;
out_dat.rev_peaks_color1 = rev_peaks_color1;
out_dat.rev_peaks_color2 = rev_peaks_color2;
out_dat.probe1_2FP = probe1_2FP;
out_dat.probe2_2FP = probe2_2FP;
out_dat.probe1_1FP = probe1_1FP;
out_dat.probe2_1FP = probe2_1FP;
out_dat.noise_color1 = noise_color1;
out_dat.noise_color2 = noise_color2;

out_dat.OneFP_ampl = OneFP_ampl;
out_dat.OneFP_color = OneFP_color;
out_dat.TwoFP_ampl = TwoFP_ampl;
out_dat.TwoFP_color = TwoFP_color;

out_dat.num_peaks = [all_peaks(1).count; all_peaks(2).count; all_peaks(3).count; all_peaks(4).count];
out_dat.all_peaks = all_peaks;
out_dat.motion = motion;
out_dat.data_bs_color1 = data_bs_color1;
out_dat.data_bs_color2 = data_bs_color2;
out_dat.thresh_curve_color1 = thresh_curve_color1;
out_dat.thresh_curve_color2 = thresh_curve_color2;


if saveOrPlot == 1
    saveToDir = varargin{1};
    save([saveToDir stem], 'out_dat');
end

end

function [TwoFP, OneFP] = TwoFP_vs_OneFP_Peaks(coincidentPeaks)
TR_GFP = 0.057;
TR_TDT = 0.089;

coincPks = [coincidentPeaks(1).pks coincidentPeaks(2).pks];% All coincident peaks as an array of peak amplitudes (column 1 - GFP, column 2 - tdTomato)
coincPksSort = sort(coincPks,2); % Sort each pair by amplitude (largest in the second column)
TwoLambdaRatio = coincPksSort(:,1) ./ coincPksSort(:,2); % The ratio of the smaller amplitudes against the larger amplitudes

% Identify GFP-dominant peaks
isGFPDom = (coincPks(:,1) >= coincPks(:,2));

is2FP = TwoLambdaRatio > ((TR_GFP .* coincPksSort(:,2) + 5) ./ (coincPksSort(:,2) - 5));
TwoFP_GFPDom(1) = filterPeaks(coincidentPeaks(1), isGFPDom & is2FP);
TwoFP_GFPDom(2) = filterPeaks(coincidentPeaks(2), isGFPDom & is2FP);
OneFP_GFPDom(1) = filterPeaks(coincidentPeaks(1), isGFPDom & ~is2FP);
OneFP_GFPDom(2) = filterPeaks(coincidentPeaks(2), isGFPDom & ~is2FP);

is2FP = TwoLambdaRatio > ((TR_TDT .* coincPksSort(:,2) + 5) ./ (coincPksSort(:,2) - 5));
TwoFP_TDTDom(1) = filterPeaks(coincidentPeaks(1), ~isGFPDom & is2FP);
TwoFP_TDTDom(2) = filterPeaks(coincidentPeaks(2), ~isGFPDom & is2FP);
OneFP_TDTDom(1) = filterPeaks(coincidentPeaks(1), ~isGFPDom & ~is2FP);
OneFP_TDTDom(2) = filterPeaks(coincidentPeaks(2), ~isGFPDom & ~is2FP);

TwoFP(1) = combinePeaks(TwoFP_GFPDom(1), TwoFP_TDTDom(1));
TwoFP(2) = combinePeaks(TwoFP_GFPDom(2), TwoFP_TDTDom(2));
OneFP(1) = combinePeaks(OneFP_GFPDom(1), OneFP_TDTDom(1));
OneFP(2) = combinePeaks(OneFP_GFPDom(2), OneFP_TDTDom(2));

end

function [] = plot_4channels_clusters(data, time, all_peaks, coinc_pks_cluster, coinc_pks_other, fwd_pks, rev_pks, removed_pks, data_color, thresh_curve)
hold on
plot(time(1:4:end), data(1:4:end), 'color', data_color, 'DisplayName', 'DiFC Data')
plot(time(all_peaks.locs), all_peaks.pks, 'o', 'MarkerSize', 2, 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0], 'DisplayName', 'All Peaks');
plot(time(coinc_pks_cluster.locs), coinc_pks_cluster.pks, 'o', 'MarkerSize', 6, 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0], 'DisplayName', '2-color Coinc Peaks')
plot(time(coinc_pks_other.locs), coinc_pks_other.pks, 'o', 'MarkerSize', 7, 'Color', [0 0 0], 'DisplayName', 'Low Ratio Coinc Peaks');
plot(time(fwd_pks.locs), fwd_pks.pks, '>', 'MarkerSize', 6, 'Color', [1 0 1], 'DisplayName', 'Forward Matches');
plot(time(rev_pks.locs), rev_pks.pks, '<', 'MarkerSize', 6, 'Color', [0 0 1], 'DisplayName', 'Reverse Matches');
plot(time(removed_pks.locs), removed_pks.pks, 'x', 'MarkerSize', 6, 'Color', [0 0 1], 'DisplayName', '4-channel Motion')
plot(time, thresh_curve, 'Color', [0 1 1], 'DisplayName', 'Peak Threshold')
hold off
end
