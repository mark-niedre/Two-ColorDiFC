function [not_coinc_peaks1, not_coinc_peaks2, not_coinc_peaks3, not_coinc_peaks4, rm_coinc_peaks1, rm_coinc_peaks2, rm_coinc_peaks3, rm_coinc_peaks4, coinc_pk_count] = remove4ChannelCoincPeaks(peaks1, peaks2, peaks3, peaks4, time, nameValArgs)
%REMOVECOINCPEAKS Returns peaks structs with peaks coincident to eachother
% removed. Also cleans up each peak struct individually
% NOTE: When multiple peaks in channel 1 are found to be "coincident"
% with a peak in another channel, the earlier of the peaks is considered to
% be the coincident one. For the rest of the channels, I loop through all
% possible combinations of peaks in channel 2-4 for each channel 1 peak. I
% take the combination with the lowest combined distance from channel 1.
%
% Currently v2


arguments
    peaks1 struct
    peaks2 struct
    peaks3 struct
    peaks4 struct
    time(:,1) double
    
    nameValArgs.CoincidenceWindow double = 0.03;
end
coinc_window = nameValArgs.CoincidenceWindow;

fs = 1 ./ (time(2)-time(1)); % Sampling frequency
cws = coinc_window * fs;

% find any peaks that match each other between the two fibers within the
% user-specified coincidence window
matches1 = false(peaks1.count, 1);
matches2 = false(peaks2.count, 1);
matches3 = false(peaks3.count, 1);
matches4 = false(peaks4.count, 1);
for ii = 1:peaks1.count
    pk1_loc = peaks1.locs(ii);
    peaks2_tmp = peaks2.locs(~matches2);
    peaks3_tmp = peaks3.locs(~matches3);
    peaks4_tmp = peaks4.locs(~matches4);
    coinc2_loc = peaks2_tmp(peaks2_tmp>pk1_loc-cws & peaks2_tmp<pk1_loc+cws);
    coinc3_loc = peaks3_tmp(peaks3_tmp>pk1_loc-cws & peaks3_tmp<pk1_loc+cws);
    coinc4_loc = peaks4_tmp(peaks4_tmp>pk1_loc-cws & peaks4_tmp<pk1_loc+cws);
    if(~isempty(coinc2_loc) && ~isempty(coinc3_loc) && ~isempty(coinc4_loc))
        coinc_fit = Inf;
        best_fit_peaks = [];
        for jj = 1:length(coinc2_loc)
            fit_value = abs(pk1_loc - coinc2_loc(jj));
            coinc3_loc2 = coinc3_loc(coinc3_loc>coinc2_loc(jj)-cws & coinc3_loc<coinc2_loc(jj)+cws);
            coinc4_loc2 = coinc4_loc(coinc4_loc>coinc2_loc(jj)-cws & coinc4_loc<coinc2_loc(jj)+cws);
            if(~isempty(coinc3_loc2) && ~isempty(coinc4_loc2))
                for kk = 1:length(coinc3_loc2)
                    fit_value = fit_value + abs(pk1_loc - coinc3_loc2(kk));
                    coinc4_loc3 = coinc4_loc2(coinc4_loc2>coinc3_loc2(kk)-cws & coinc4_loc2<coinc3_loc2(kk)+cws);
                    if(~isempty(coinc4_loc3))
                        fit_value = fit_value + abs(pk1_loc - coinc4_loc3);
                        [fit_value, ind] = min(fit_value);
                        if fit_value < coinc_fit
                            coinc_fit = fit_value;
                            best_fit_peaks = [coinc2_loc(jj) coinc3_loc2(kk) coinc4_loc3(ind)];
                        end
                    end
                end
            end
        end
        if coinc_fit < Inf
            matches1(ii) = true;
            matches2(peaks2.locs == best_fit_peaks(1)) = true;
            matches3(peaks3.locs == best_fit_peaks(2)) = true;
            matches4(peaks4.locs == best_fit_peaks(3)) = true;
        end
    end
end

% remove all matches from peaks1
not_coinc_peaks1 = filterPeaks(peaks1, ~matches1);
not_coinc_peaks2 = filterPeaks(peaks2, ~matches2);
not_coinc_peaks3 = filterPeaks(peaks3, ~matches3);
not_coinc_peaks4 = filterPeaks(peaks4, ~matches4);
% Also return the 4-channel coincident matches (aka motion artifacts)
rm_coinc_peaks1 = filterPeaks(peaks1, matches1);
rm_coinc_peaks2 = filterPeaks(peaks2, matches2);
rm_coinc_peaks3 = filterPeaks(peaks3, matches3);
rm_coinc_peaks4 = filterPeaks(peaks4, matches4);

% Count the number of peaks removed
coinc_pk_count = sum(matches1) + sum(matches2) + sum(matches3) + sum(matches4);

end

