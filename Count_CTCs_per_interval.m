
function [detections_per_interval] = Count_CTCs_per_interval(detection_times, interval_length, scan_length, fs)
% detection_times in seconds
% interval_length in seconds
% scan_length in seconds
% fs in Hz
if ~isempty(detection_times) && detection_times(end) > scan_length
    error('1 or more detection times exceed the input scan length.')
end
detections = zeros(round(scan_length*fs),1);
detections(round(detection_times * fs)) = 1;
detections_per_interval = movsum(detections, interval_length*fs);
% Take only full-sized intervals (remove intervals that only count portions
% of the beginning and end of scan)
detections_per_interval = detections_per_interval(interval_length*fs/2:end-interval_length*fs/2);
end