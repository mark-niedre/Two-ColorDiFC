%% Phantom MM and 4T1 SNR
load('Phantom_MM_GFP', 'out_dat'); out_dat_gfp_MM = out_dat;
load('Phantom_MM_TDT', 'out_dat'); out_dat_tdt_MM = out_dat;
gfp_snr_MM = [ 20*log10(out_dat_gfp_MM.fwd_peaks_color1(1).pks./ out_dat_gfp_MM.noise_color1(1)); 20*log10(out_dat_gfp_MM.fwd_peaks_color1(2).pks./ out_dat_gfp_MM.noise_color1(2))];
tdt_snr_MM = [ 20*log10(out_dat_tdt_MM.fwd_peaks_color2(1).pks./ out_dat_tdt_MM.noise_color2(1)); 20*log10(out_dat_tdt_MM.fwd_peaks_color2(2).pks./ out_dat_tdt_MM.noise_color2(2))];

load('Phantom_4T1_GFP_3', 'out_dat'); out_dat_gfp_4T1 = out_dat;
load('Phantom_4T1_TDT_4', 'out_dat'); out_dat_tdt_4T1 = out_dat;
gfp_snr_4T1 = [ 20*log10(out_dat_gfp_4T1.fwd_peaks_color1(1).pks ./ out_dat_gfp_4T1.noise_color1(1)); 20*log10(out_dat_gfp_4T1.fwd_peaks_color1(2).pks ./ out_dat_gfp_4T1.noise_color1(2)) ];
tdt_snr_4T1 = [ 20*log10(out_dat_tdt_4T1.fwd_peaks_color2(1).pks ./ out_dat_tdt_4T1.noise_color2(1)); 20*log10(out_dat_tdt_4T1.fwd_peaks_color2(2).pks ./ out_dat_tdt_4T1.noise_color2(2)) ];

%% 1 plot - SNR based on avg of both channels
gfp_color = .5*[0 .5 .22];
tdt_color = [.98 .48 .67];

f1 = figure;
histogram(gfp_snr_MM, 'Normalization', 'probability', 'BinWidth', 2.5, 'DisplayName', 'GFP', 'FaceColor', gfp_color)
xlabel('SNR (dB)');
hold on
histogram(tdt_snr_MM, 'Normalization', 'probability', 'BinWidth', 2.5, 'DisplayName', 'tdTomato', 'FaceColor', tdt_color)
set(gca, 'FontSize', 16, 'LineWidth', 1)
ylim([0 .32]);
xlim([0 60])
legend('Location', 'northeast')
ylabel('Probability');

a = get(f1, 'Position');
set(f1, 'Position', [a(1:2) 397.6 297.6])

%% SNR of each channel
f2 = figure('DefaultAxesFontSize',20);
histogram(gfp_snr_4T1, 'Normalization', 'probability', 'BinWidth', 2.5, 'DisplayName', 'GFP', 'FaceColor', gfp_color)
xlabel('SNR (dB)');
hold on
histogram(tdt_snr_4T1, 'Normalization', 'probability', 'BinWidth', 2.5, 'DisplayName', 'tdTomato', 'FaceColor', tdt_color)
set(gca, 'FontSize', 16, 'LineWidth', 1)
ylim([0 .2]);
xlim([0 60])
ylabel('Probability');

a = get(f2, 'Position');
set(f2, 'Position', [a(1:2) 397.6 297.6])