%% Plot MM CTCs Over Time
processedDir = 'Two-colorDiFC_Data\';

%% Load Data
MM_Mouse1_files = {'MM_Mouse1_8' 'MM_Mouse1_7' 'MM_Mouse1_5' 'MM_Mouse1_4' 'MM_Mouse1_3'};
difc_days_MM = [43 39 37 35 31];
difc_day_labels_MM = {'43' '39' '37'  '35' '31'};

detection_times_gfp = cell(1,length(difc_days_MM));
detection_times_tdt = cell(1,length(difc_days_MM));
scan_lengths = zeros(1,length(difc_days_MM));
avgs_gfp = scan_lengths;
avgs_tdt = scan_lengths;
for i = 1:length(MM_Mouse1_files)
    load([processedDir MM_Mouse1_files{i}], 'out_dat');
    avgs_gfp(i) = 20 * out_dat.fwd_peaks_color1(1).count / (out_dat.scan_length_minutes);
    avgs_tdt(i) = 20 * out_dat.fwd_peaks_color2(1).count / (out_dat.scan_length_minutes);
    detection_times_gfp{i} = out_dat.fwd_peaks_color1(1).locs./2000;
    detection_times_tdt{i} = out_dat.fwd_peaks_color2(1).locs./2000;
    scan_lengths(i) = out_dat.scan_length_minutes;
end
trial_length = max(scan_lengths);

%% Plot Rasters
gfp_color = [0 .69 .31];
tdt_color = [.8 .38 .54];
fig = figure('DefaultAxesFontSize',14);

% Two rasters
% GFP
t = tiledlayout(1,3);
t.TileSpacing = 'compact';
nexttile(1);
rasterplot(detection_times_gfp, trial_length, difc_day_labels_MM, '', 'Days after MM Injection', gfp_color);
set(gca, 'linewidth',1, 'XTick', [0 45])
xlim([0 45.1]);
title('GFP CTCs');
% tdTomato
nexttile(2);
rasterplot(detection_times_tdt, trial_length, difc_day_labels_MM, 'Time (min)', '', tdt_color);
set(gca, 'linewidth',1, 'XTick', [0 45], 'YTickLabel', {})
xlim([0 45.1]);
title('tdTomato CTCs');

% One raster
nexttile();
rasterplot_2colors(detection_times_gfp, detection_times_tdt, gfp_color, tdt_color, trial_length, difc_day_labels_MM, '', '');
set(gca, 'linewidth',1, 'XTick', [0 45], 'YTickLabel', {})
xlim([0 45.1]);
title('All CTCs');

%% Plot CTC Detections
figure('DefaultAxesFontSize',14);
plot(fliplr(difc_days_MM), fliplr(avgs_gfp),'o-', 'Color', gfp_color, 'DisplayName', 'GFP', 'MarkerFaceColor', gfp_color, 'MarkerSize', 12, 'LineWidth', 3);
xlabel('Days after MM Injection', 'FontSize', 18);
ylabel('CTCs/mL PB', 'FontSize', 18);
hold on
plot(fliplr(difc_days_MM), fliplr(avgs_tdt),'x:', 'Color', tdt_color, 'DisplayName', 'tdTomato', 'MarkerEdgeColor', tdt_color, 'MarkerSize', 13, 'LineWidth', 3);
ylim([0 70]);
xlim([27 45]);
set(gca, 'linewidth',1)
legend('Location', 'northwest');

