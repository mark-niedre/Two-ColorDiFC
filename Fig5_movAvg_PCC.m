%% Plot Moving Avg
processedDir = 'Two-colorDiFC_Data\';
MM_files = {'MM_Mouse1_1' 'MM_Mouse1_2' 'MM_Mouse1_3' 'MM_Mouse1_4' 'MM_Mouse1_5' 'MM_Mouse1_6' 'MM_Mouse1_7' 'MM_Mouse1_8'...
    'MM_Mouse2_1' 'MM_Mouse2_2' 'MM_Mouse2_3' 'MM_Mouse2_4' 'MM_Mouse2_5'...
    'MM_Mouse3_1' 'MM_Mouse3_2' 'MM_Mouse3_3' 'MM_Mouse3_4'...
    'MM_Mouse4_1' 'MM_Mouse4_2' 'MM_Mouse4_3' 'MM_Mouse4_4'};

% Identifying the files with the minimum CTC detection rate of 0.5 CTCs/min
inds_min_rate = [8 10 11 12 13 16 17 19 21];
MM_files = MM_files(inds_min_rate);

% Percents of the blood volume to use
percent_bloodVols = [1 5 10 20];
% Convert to time (of DiFC scan) assuming 50 uL per minute and 2000uL total blood volume
interval_lengths = percent_bloodVols .* 120 ./ 5;
interval_length = interval_lengths(2);
fs = 2000;

%% Calculate Pearson Correlation Coefficients of Moving Averages (PCCs)
PCCs = NaN(length(MM_files),4);
p_vals = NaN(length(MM_files),4);
num_gfp = NaN(length(MM_files),1);
num_tdt = NaN(length(MM_files),1);
scan_lengths = NaN(length(MM_files),1);

for i = 1:length(MM_files)
    load([processedDir MM_files{i}], 'out_dat');
    gfp_times = sort(out_dat.fwd_peaks_color1(1).locs) ./fs;
    tdt_times = sort(out_dat.fwd_peaks_color2(1).locs) ./fs;
    
    num_gfp(i) = length(gfp_times);
    num_tdt(i) = length(tdt_times);
    scan_lengths(i) = out_dat.scan_length_minutes;
    
    for j = 1:4
        gfp_per_interval = Count_CTCs_per_interval(gfp_times, interval_lengths(j), out_dat.scan_length_samples./fs, fs);
        tdt_per_interval = Count_CTCs_per_interval(tdt_times, interval_lengths(j), out_dat.scan_length_samples./fs, fs);
        
        % Down sampling
        gfp_per_interval = gfp_per_interval(1:2000:end);
        tdt_per_interval = tdt_per_interval(1:2000:end);
        gfp_per_interval = gfp_per_interval(1:60:end);
        tdt_per_interval = tdt_per_interval(1:60:end);
        
        [A_tmp, p_tmp] = corrcoef(gfp_per_interval,tdt_per_interval);
        PCCs(i,j) = A_tmp(2);
        p_vals(i,j) = p_tmp(2);
    end
end

%% PCC vs CTC Detection Rate
figure('DefaultAxesFontSize', 14);
set(gcf, 'Position', 1000 * [2.8802   -0.0238    0.3832    0.3744])
% set(gcf, 'Position', 1000 * [2.4394   -0.0182    0.4656    0.3744])
hold on;
ctc_rate = (num_gfp + num_tdt) ./ scan_lengths;
inds_p_big = p_vals(:,2) >= 0.05;
inds_p_005 = (p_vals(:,2) >= 0.01) & (p_vals(:,2) < 0.05);
inds_p_001 = (p_vals(:,2) >= 0.001) & (p_vals(:,2) < 0.01);
inds_p_0001 = (p_vals(:,2) < 0.001);
scatter(ctc_rate, PCCs(:,2), 100,'MarkerFaceColor', [0 0 0], 'MarkerFaceAlpha', .4,'MarkerEdgeColor', [0 0 0], 'Marker', 'o')
plot([0 8], [0 0], 'k-', 'HandleVisibility', 'off')
ylim([-.5 1])
xlim([0 8])
box on
xlabel('CTC Detection Rate (cells/min)')
ylabel('Pearson Correlation Coefficient')

%% Plot Moving Averages for File 1
load([processedDir MM_files{1}], 'out_dat');
gfp_times = sort(out_dat.fwd_peaks_color1(1).locs) ./fs;
tdt_times = sort(out_dat.fwd_peaks_color2(1).locs) ./fs;
gfp_per_interval = Count_CTCs_per_interval(gfp_times, interval_lengths(2), out_dat.scan_length_samples./fs, fs);
tdt_per_interval = Count_CTCs_per_interval(tdt_times, interval_lengths(2), out_dat.scan_length_samples./fs, fs);
gfp_per_interval = gfp_per_interval(1:2000:end);
tdt_per_interval = tdt_per_interval(1:2000:end);
gfp_per_interval = gfp_per_interval(1:60:end);
tdt_per_interval = tdt_per_interval(1:60:end);

avg_gfp = length(gfp_times)./ out_dat.scan_length_minutes;
avg_tdt = length(tdt_times)./ out_dat.scan_length_minutes;

time = linspace(0,out_dat.scan_length_samples./fs, out_dat.scan_length_samples-1)';
time_movavg = time(interval_length*fs/2:end-interval_length*fs/2+1)/60;
time_movavg = time_movavg(1:2000:end);
time_movavg = time_movavg(1:60:end);

plot2colorMovAvg(gfp_times, tdt_times, gfp_per_interval, tdt_per_interval, avg_gfp, avg_tdt, interval_length,  time_movavg, out_dat.scan_length_minutes);
set(gcf, 'Position', 1000 * [2.4386    0.3610    0.8760    0.4176])

plot2colorMovAvgScatter(gfp_per_interval, tdt_per_interval);
set(gcf, 'Position', 1000 * [2.4394   -0.0182    0.4656    0.3744])

%%
function [] = plot2colorMovAvgScatter(gfp_per_interval, tdt_per_interval)
figure('DefaultAxesFontSize',14);
hold on;
scatter(gfp_per_interval, tdt_per_interval, 'MarkerFaceColor', .3*[1 1 1], 'MarkerEdgeColor', [0 0 0])
a = 15;
% a = 1.1*max(max(gfp_per_interval), max(tdt_per_interval));
% plot([0 a], [0 a], 'k--');
axis equal
axis([0 a 0 a ]);
box on
ylabel('tdTomato CTCs per 2 Min')
xlabel('GFP CTCs per 2 Min')
end

function [] = plot2colorMovAvg(gfp_times, tdt_times, gfp_per_interval, tdt_per_interval, avg_gfp, avg_tdt, interval_length, time_movavg, scan_length)
gfp_color = [0 .69 .31];
tdt_color = [.8 .38 .54];

x_axis_span = [0.01 scan_length + 0.01];


avg_gfp_movavg = avg_gfp * interval_length / 60;
num_bins_gfp = ceil(max(gfp_per_interval))+1;
xl2_g = 1.1*(num_bins_gfp-1);

avgs_tdt_movavg = avg_tdt * interval_length / 60;
num_bins_tdt = ceil(max(tdt_per_interval))+1;
xl2_t = 1.1*(num_bins_tdt-1);

figure('DefaultAxesFontSize',14);
tiledlayout(3,1);

nexttile(1);
hold on;
if ~isempty(gfp_times)
    plot([gfp_times/60 gfp_times/60]',[0.55 1]','-', 'HandleVisibility', 'off', 'Color', gfp_color);
end

if ~isempty(tdt_times)
    plot([tdt_times/60 tdt_times/60]',[0 0.45]','-', 'HandleVisibility', 'off', 'Color', tdt_color);
end

ylim([-0.12 1.12]);
xlim([0 scan_length + 0.01]);
set(gca,'ytick',[],'yticklabel',[],'xtick',[0 floor(scan_length)], 'xticklabel', [])
ylabel(sprintf('CTC\nDetections'));
box on
nexttile(2);
hold on
plot(x_axis_span, [avg_gfp_movavg avg_gfp_movavg], 'k-', 'LineWidth', 1.5);
plot(time_movavg, gfp_per_interval','-','LineWidth', 1.5, 'Color', gfp_color);
xlim([0 scan_length + 0.01]);
set(gca,'xtick',[0 floor(scan_length)], 'xticklabel', [])
box on
ylim([0 xl2_g])
ylabel(sprintf('GFP CTCs\nper 2 Min'))

nexttile(3);
hold on
plot(x_axis_span, [avgs_tdt_movavg avgs_tdt_movavg], 'k-', 'LineWidth', 1.5);
plot(time_movavg, tdt_per_interval','-','LineWidth', 1.5, 'Color', tdt_color);
xlim([0 scan_length + 0.01]);
set(gca,'xtick',[0 floor(scan_length)])
box on
ylim([0 xl2_t])
ylabel(sprintf('TDT CTCs\nper 2 Min'))
xlabel('Time (min)')
end


