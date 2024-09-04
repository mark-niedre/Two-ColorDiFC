%% Identify and plot 2FP and 1FP detections from 2lambda peaks
processedDir = 'Two-colorDiFC_Data\';

%% MM Phantom Data (Single Cells)
MMPhantom_GFP = 'Phantom_MM_GFP';
MMPhantom_tdTomato = 'Phantom_MM_TDT';
[MMPhantomData_GFP] = load2ColorData({MMPhantom_GFP}, processedDir);
[MMPhantomData_tdT] = load2ColorData({MMPhantom_tdTomato}, processedDir);
MM_Phantom_1stQ_gfp_ampl = quantile([MMPhantomData_GFP.fwd_pk_ampl_gfp(:,1); MMPhantomData_GFP.fwd_pk_ampl_gfp(:,2)],.25);
MM_Phantom_1stQ_tdt_ampl = quantile([MMPhantomData_tdT.fwd_pk_ampl_tdt(:,1); MMPhantomData_tdT.fwd_pk_ampl_tdt(:,2)],.25);

% Plot 2lambda Peaks - GFP and TDT Single Cells
Plot2LambdaPeaks(MMPhantomData_GFP, [MM_Phantom_1stQ_gfp_ampl MM_Phantom_1stQ_tdt_ampl],'GFP MM Phantom Data');
axis([0 600 0 600]);
a = get(gcf, 'Position'); set(gcf, 'Position', [a(1:2) 389.2 288]);
Plot2LambdaPeaks(MMPhantomData_tdT, [MM_Phantom_1stQ_gfp_ampl MM_Phantom_1stQ_tdt_ampl],'tdTomato MM Phantom Data');
axis([0 600 0 600]);
a = get(gcf, 'Position'); set(gcf, 'Position', [a(1:2) 389.2 288]);

%% 4T1 Phantom Data (Clusters)
BC4T1_Phantom_2cCluster = {'Phantom_4T1_2Color_Clusters'};
BC4T1_Phantom_both_cluster = {'Phantom_4T1_Both_GFP_TDT_1' 'Phantom_4T1_Both_GFP_TDT_2' 'Phantom_4T1_Both_GFP_TDT_3' 'Phantom_4T1_Both_GFP_TDT_4' };
BC4T1_Phantom_GFP = {'Phantom_4T1_GFP_1' 'Phantom_4T1_GFP_2'};
BC4T1_Phantom_TDT = {'Phantom_4T1_TDT_1' 'Phantom_4T1_TDT_2' 'Phantom_4T1_TDT_3'};

[BC4T1_GFP] = load2ColorData(BC4T1_Phantom_GFP, processedDir);
[BC4T1_tdT] = load2ColorData(BC4T1_Phantom_TDT, processedDir);
BC4T1_Phantom_1stQ_gfp_ampl = quantile([BC4T1_GFP.fwd_pk_ampl_gfp(:,1); BC4T1_GFP.fwd_pk_ampl_gfp(:,2)],.25);
BC4T1_Phantom_1stQ_tdt_ampl = quantile([BC4T1_tdT.fwd_pk_ampl_tdt(:,1); BC4T1_tdT.fwd_pk_ampl_tdt(:,2)],.25);

% Plot 2lambda Peaks - GFP-only or TDT-only 1FP Clusters
Plot2LambdaPeaks(BC4T1_GFP, [0 0],'GFP 4T1 Phantom Data');
axis([0 600 0 600]);
a = get(gcf, 'Position'); set(gcf, 'Position', [a(1:2)  389.2 288]);
Plot2LambdaPeaks(BC4T1_tdT, [0 0],'tdTomato 4T1 Phantom Data');
axis([0 600 0 600]);
a = get(gcf, 'Position'); set(gcf, 'Position', [a(1:2)  389.2 288]);

% Plot 2lambda Peaks - Mixed 1FP Clusters
[BC4T1_both] = load2ColorData(BC4T1_Phantom_both_cluster, processedDir);
% [~, ~, G_tmp, T_tmp, G_1FP_3cells, T_1FP_3cells] = Plot2LambdaPeaks(BC4T1_both, [0 0],'4T1 Phantom Data (Both)');
Plot2LambdaPeaks(BC4T1_both, [0 0],'4T1 Phantom Data (Both)');
axis([0 600 0 600]);
a = get(gcf, 'Position'); set(gcf, 'Position', [a(1:2)  389.2 288]);

% Plot 2lambda Peaks - 2FP Clusters
[BC4T1_2c] = load2ColorData(BC4T1_Phantom_2cCluster, processedDir);
% [~, ~, G_tmp, T_tmp, G_1FP_3cells, T_1FP_3cells] = Plot2LambdaPeaks(BC4T1_2c, [0 0],'4T1 Phantom Data (2 Color)');
Plot2LambdaPeaks(BC4T1_2c, [0 0],'4T1 Phantom Data (2 Color)');
axis([0 600 0 600]);
a = get(gcf, 'Position'); set(gcf, 'Position', [a(1:2)  389.2 288]);

%% MM Mice
MM_files = {'MM_Mouse1_1' 'MM_Mouse1_2' 'MM_Mouse1_3' 'MM_Mouse1_4' 'MM_Mouse1_5' 'MM_Mouse1_6' 'MM_Mouse1_7' 'MM_Mouse1_8'...
    'MM_Mouse2_1' 'MM_Mouse2_2' 'MM_Mouse2_3' 'MM_Mouse2_4' 'MM_Mouse2_5'...
    'MM_Mouse3_1' 'MM_Mouse3_2' 'MM_Mouse3_3' 'MM_Mouse3_4'...
    'MM_Mouse4_1' 'MM_Mouse4_2' 'MM_Mouse4_3' 'MM_Mouse4_4'};

[allMMMouseData] = load2ColorData(MM_files, processedDir);
MM_Mouse_1stQ_gfp_ampl = quantile([allMMMouseData.fwd_pk_ampl_gfp(:,1); allMMMouseData.fwd_pk_ampl_gfp(:,2)],.25);
MM_Mouse_1stQ_tdt_ampl = quantile([allMMMouseData.fwd_pk_ampl_tdt(:,1); allMMMouseData.fwd_pk_ampl_tdt(:,2)],.25);

% Manually identified false 2FP detections
high_ratio_keep = [1 2 4 5 8 9 16 17 18 19 21 22 24 25 26 28 29 30 32 33 34 36 39 44 45 46 47 48 50 51 52 54 57 60 63];
allMMMouseData.TwoFP_ampl = allMMMouseData.TwoFP_ampl(high_ratio_keep,:);
allMMMouseData.TwoFP_color = allMMMouseData.TwoFP_color(high_ratio_keep);

% Plot 2lambda Peaks - MM Mice
[gfp_2FP, tdt_2FP, ~, ~, G_1FP_3cells, T_1FP_3cells] = Plot2LambdaPeaks(allMMMouseData, [MM_Mouse_1stQ_gfp_ampl MM_Mouse_1stQ_tdt_ampl],'MM Mouse Data');
axis([0 600 0 600])
Plot2LambdaPeaks(allMMMouseData, [MM_Mouse_1stQ_gfp_ampl MM_Mouse_1stQ_tdt_ampl],'MM Mouse Data');
axis([0 100 0 100]);

% Plot histogram of estimated 2FP sizes
figure('DefaultAxesFontSize', 14); 
histogram([gfp_2FP; tdt_2FP], 'BinWidth', 3, 'FaceColor', .3*[1 1 1]);
xlabel('Estimated Cluster Size (# cells)'); ylabel('Number of Detections');
set(gca, 'LineWidth',1); 
a = get(gcf, 'Position'); set(gcf, 'Position', [a(1:2)  384.0000  344.4000])

% Plot histogram of estimated 1FP sizes
gfp_color = .5*[0 .5 .22];
tdt_color = [.98 .48 .67];
figure('DefaultAxesFontSize', 14); hold on
histogram(G_1FP_3cells, 'BinWidth', 3, 'DisplayName', 'GFP', 'FaceColor', gfp_color);
histogram(T_1FP_3cells, 'BinWidth', 3, 'DisplayName', 'tdTomato', 'FaceColor', tdt_color);
xlabel('Estimated Cluster Size (# cells)'); ylabel('Number of Detections');
set(gca, 'LineWidth',1); legend; 
a = get(gcf, 'Position'); set(gcf, 'Position', [a(1:2)  384.0000  344.4000]); box on

% Mean and median 2FP and 1FP sizes
fprintf('2FP\n\tMean 2FP size: %f\n', mean([gfp_2FP; tdt_2FP]));
fprintf('\tMedian 2FP size: %f\n', median([gfp_2FP; tdt_2FP]));
fprintf('1FP\n\tMean 1FP size: %f\n', mean([G_1FP_3cells; T_1FP_3cells]));
fprintf('\tMedian 1FP size: %f\n', median([G_1FP_3cells; T_1FP_3cells]));

%% Control Mice
ctrlDataFiles = {'NSG_Control_1' 'NSG_Control_2' 'NSG_Control_3' 'NSG_Control_4'};
[Control_data_tail] = load2ColorData(ctrlDataFiles, processedDir);

% Plot 2lambda Peaks - Control Mice
Plot2LambdaPeaks(Control_data_tail, [MM_Mouse_1stQ_gfp_ampl MM_Mouse_1stQ_tdt_ampl],'Control Mice');

%%
% Plot GFP and tdTomato Coincident Peak Ratios vs Amplitude
function [GFP_2FP_sizes, tdTomato_2FP_sizes, GFP_2FP_sizes_g_t,...
    tdTomato_2FP_sizes_g_t, GFP_1FP_3cells_sizes, tdTomato_1FP_3cells_sizes] = ...
    Plot2LambdaPeaks(allData, singleCellAmpl, fig_title)
figure('Name',fig_title);
fprintf('\n%s\n', fig_title);
gfp_color = [0 .69 .31];
tdt_color = [.8 .38 .54];

% For the sake of accuracy, reidentify 2FP and 1FP peaks
all_2lambda = [allData.TwoFP_ampl; allData.OneFP_ampl]; % All 2lambda peaks
sorted_2lambda = sort(all_2lambda, 2); % Larger ampl peaks in column 2

disp1 = size(all_2lambda,1);
fprintf('All 2lambda: %d\n', disp1);

% Limiting to 1500 mV amplitude peaks due PMT saturation above ~3000 mV
all_2lambda = all_2lambda(sorted_2lambda(:,2) < 1500,:);
sorted_2lambda = sorted_2lambda(sorted_2lambda(:,2) < 1500,:);

disp2 = size(all_2lambda,1);
fprintf('Larger than 1500 mV: %d\n', disp1-disp2);
fprintf('2lambda (-1500mV): %d\n', disp2);

gfp_2lambda_flag = (all_2lambda(:,1) > all_2lambda(:,2)); % Identifies GFP-dominant peaks
ratio_2lambda = sorted_2lambda(:,1) ./ sorted_2lambda(:,2); %small peak ampl / large peak ampl

gfp_2lambda_ampl_large = all_2lambda(gfp_2lambda_flag, 1);
gfp_2lambda_ampl_small = all_2lambda(gfp_2lambda_flag, 2);
gfp_2lambda_ratio = ratio_2lambda(gfp_2lambda_flag);

tdt_2lambda_ampl_large = all_2lambda(~gfp_2lambda_flag, 2);
tdt_2lambda_ampl_small = all_2lambda(~gfp_2lambda_flag, 1);
tdt_2lambda_ratio = ratio_2lambda(~gfp_2lambda_flag);

% Threshold for 2FP detections
x = 11:1:1580;
TR_GFP = 0.057;
TR_TDT = 0.089;
y_GFP = (TR_GFP .* x + 5) ./ (x - 5);
y_TDT = (TR_TDT .* x + 5) ./ (x - 5);

% Find coincident detection that are consistent with bleed
% 1FP
gfp_2lambda_bleed_flag = gfp_2lambda_ratio <= ((TR_GFP .* gfp_2lambda_ampl_large + 5) ./ (gfp_2lambda_ampl_large - 5));
tdt_2lambda_bleed_flag = tdt_2lambda_ratio <= ((TR_TDT .* tdt_2lambda_ampl_large + 5) ./ (tdt_2lambda_ampl_large - 5));

at_least_3_gfp = (all_2lambda(:,1) >= 3*singleCellAmpl(1));
at_least_3_tdt = (all_2lambda(:,2) >= 3*singleCellAmpl(2));

% 2FP
gfp_2FP_flag = ~gfp_2lambda_bleed_flag;
tdt_2FP_flag = ~tdt_2lambda_bleed_flag;

% 1FP, >=3 cells
gfp_1FP_3cells = gfp_2lambda_bleed_flag & at_least_3_gfp(gfp_2lambda_flag);
tdt_1FP_3cells = tdt_2lambda_bleed_flag & at_least_3_tdt(~gfp_2lambda_flag);

% 1FP, <3 cells
gfp_1FP_small = gfp_2lambda_bleed_flag & ~at_least_3_gfp(gfp_2lambda_flag);
tdt_1FP_small = tdt_2lambda_bleed_flag & ~at_least_3_tdt(~gfp_2lambda_flag);

% GFP
X_2FP_gfp = [gfp_2lambda_ampl_large(gfp_2FP_flag) gfp_2lambda_ampl_small(gfp_2FP_flag)];
X_1FP_3cells_gfp = [gfp_2lambda_ampl_large(gfp_1FP_3cells) gfp_2lambda_ampl_small(gfp_1FP_3cells)];
X_1FP_small_gfp = [gfp_2lambda_ampl_large(gfp_1FP_small) gfp_2lambda_ampl_small(gfp_1FP_small)];

% tdTomato
X_2FP_tdt = [tdt_2lambda_ampl_small(tdt_2FP_flag) tdt_2lambda_ampl_large(tdt_2FP_flag)];
X_1FP_3cells_tdt = [tdt_2lambda_ampl_small(tdt_1FP_3cells) tdt_2lambda_ampl_large(tdt_1FP_3cells)];
X_1FP_small_tdt = [tdt_2lambda_ampl_small(tdt_1FP_small) tdt_2lambda_ampl_large(tdt_1FP_small)];

fprintf('\t2FP: %d\n', size(X_2FP_gfp,1)+size(X_2FP_tdt,1));
fprintf('\t1FP 3Cells: %d\n', size(X_1FP_3cells_gfp,1)+size(X_1FP_3cells_tdt,1));
fprintf('\t1FP Small: %d\n', size(X_1FP_small_gfp,1)+size(X_1FP_small_tdt,1));

nexttile();
plot2LambdaScatter(X_2FP_gfp, X_1FP_3cells_gfp, X_1FP_small_gfp, gfp_color)
plot2LambdaScatter(X_2FP_tdt, X_1FP_3cells_tdt, X_1FP_small_tdt, tdt_color)
% Plot threshold for 2FP peaks
plot(x, y_GFP.*x, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off')
plot(y_TDT.*x, x, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off')
plot(3.*[singleCellAmpl(1) singleCellAmpl(1)], [0 3.*singleCellAmpl(1).*((TR_GFP .* 3.* singleCellAmpl(1) + 5) ./ (3.* singleCellAmpl(1) - 5))], '-', 'Color', .2*[1 1 1], 'LineWidth',2, 'HandleVisibility', 'off')
plot([0 3.*singleCellAmpl(2).*((TR_TDT .* 3.* singleCellAmpl(2) + 5) ./ (3.* singleCellAmpl(2) - 5))], 3.*[singleCellAmpl(2) singleCellAmpl(2)], '-', 'Color', .2*[1 1 1], 'LineWidth',2, 'HandleVisibility', 'off')
xlabel('Green Peak Amplitude (mV)', 'FontSize', 16); ylabel('Orange Peak Amplitude (mV)', 'FontSize', 16);

% % 2FP Cluster sizes
% % GFP
GFP_2FP_sizes_g_t = [ floor(X_2FP_gfp(:,1) ./ singleCellAmpl(1)) floor(X_2FP_gfp(:,2) ./ singleCellAmpl(2)) ];
if isempty(GFP_2FP_sizes_g_t)
    GFP_2FP_sizes_g_t = double.empty(0,2);
end
GFP_2FP_sizes = GFP_2FP_sizes_g_t(:,1) + GFP_2FP_sizes_g_t(:,2);
% tdTomato
tdTomato_2FP_sizes_g_t = [ floor(X_2FP_tdt(:,2) ./ singleCellAmpl(2)) floor(X_2FP_tdt(:,1) ./ singleCellAmpl(1)) ];
if isempty(tdTomato_2FP_sizes_g_t)
    tdTomato_2FP_sizes_g_t = double.empty(0,2);
end
tdTomato_2FP_sizes = tdTomato_2FP_sizes_g_t(:,1) + tdTomato_2FP_sizes_g_t(:,2);

% % 1FP Cluster sizes
% % GFP
if isempty(X_1FP_3cells_gfp)
    GFP_1FP_3cells_sizes = double.empty(0,2);
else
    GFP_1FP_3cells_sizes = floor(X_1FP_3cells_gfp(:,1) ./ singleCellAmpl(1));
end
% tdTomato
if isempty(X_1FP_3cells_tdt)
    tdTomato_1FP_3cells_sizes = double.empty(0,2);
else
    tdTomato_1FP_3cells_sizes = floor(X_1FP_3cells_tdt(:,2) ./ singleCellAmpl(2));
end

end

function [] = plot2LambdaScatter(X_2FP, X_1FP_3cells, X_1FP_small, marker_color)
%Plot scatter points for 2FP and 1FP detections
if size(X_2FP, 2) == 0
    X_2FP = double.empty(0,2);
end
if size(X_1FP_3cells, 2) == 0
    X_1FP_3cells = double.empty(0,2);
end
if size(X_1FP_small, 2) == 0
    X_1FP_small = double.empty(0,2);
end

hold on
scatter(X_2FP(:,1), X_2FP(:,2), 80, 'Marker', 'o', 'MarkerEdgeColor', [.28 .61 .9], 'MarkerFaceColor', [.28 .61 .9], 'MarkerFaceAlpha', .3);
scatter(X_1FP_3cells(:,1), X_1FP_3cells(:,2), 20, 'Marker', 'o', 'MarkerEdgeColor', marker_color, 'MarkerFaceColor', marker_color, 'MarkerFaceAlpha', .3);
scatter(X_1FP_small(:,1), X_1FP_small(:,2), 20, 'Marker', '^', 'MarkerEdgeColor',marker_color, 'MarkerFaceColor', marker_color, 'MarkerFaceAlpha', .3);

box on
set(gca, 'FontSize', 16, 'LineWidth', 1.5)
plot([0 1600], [0 1600], ':',  'Color', .2*[1 1 1], 'LineWidth', 1.5, 'HandleVisibility', 'off')
axis equal
end
