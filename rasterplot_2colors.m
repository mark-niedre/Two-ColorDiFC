% Original code for rasterplot.m by:
% Rajiv Narayan
% askrajiv@gmail.com
% Boston University, Boston, MA
%
% Edited by Amber Williams Jan11,2019
% Last edited Oct28, 2019

function rasterplot_2colors(event_times1, event_times2, color1, color2, trial_length, day_labels, x_axis_label, y_axis_label)
% INPUTS:
% event_times: cell of arrays of events (in samples) (# trials x #events
%       per trial) eg {[1 40 550] [2] [60 400 440 500]};
% trial_length: length of the longest trial (time in minutes)
% day_labels: cell of strings that will label each day of raster plots. 
%       eg: {'1' '4' '10'}

trial_length = trial_length * 60;
numtrials = size(event_times1,2);

% Plot variables %
plotwidth=1;     % spike thickness
trialgap=1.25;    % distance between trials

[xx1, yy1, xl] = raster_func(event_times1, numtrials, trial_length);
[xx2, yy2,~] = raster_func(event_times2, numtrials, trial_length);

yyaxis left;
axes(gca);
hold on;
plot(xx1, yy1, '-', 'Color', color1, 'linewidth',plotwidth);
plot(xx2, yy2, '-', 'Color', color2, 'linewidth',plotwidth);


axis ([xlim,0,(numtrials)*trialgap]);
set(gca, 'ytick',0.5:trialgap:(trialgap*numtrials+0.5),'tickdir','out');
yticklabels(day_labels);
yl = ylim;
ylim([(1-trialgap) yl(2)]);
ylabel(y_axis_label, 'FontSize', 24);
set(gca,'YColor','k')
yyaxis right;
set(gca, 'ytick',0.5:trialgap:(trialgap*numtrials+0.5),'tickdir','in');
ylim([(1-trialgap) yl(2)]);
yticklabels([]);
set(gca,'YColor','k')

xlim(xl);
xlabel(x_axis_label, 'FontSize', 24);
box on;

end

function [xx, yy, xl] = raster_func(event_times, numtrials, trial_length)


times = [];
for i = 1:numtrials
    x = event_times{i}';
    if size(x,2) > 0
        times = [times (x + trial_length.*(i-1))];
    end
end
trialgap=1.25;    % distance between trials

trials=ceil(times/trial_length);
reltimes=mod(times,trial_length);
reltimes(~reltimes)=trial_length;
numspikes=length(times);
xx=ones(3*numspikes,1)*nan;
yy=ones(3*numspikes,1)*nan;

yy(1:3:3*numspikes)=(trials-1)*trialgap;
yy(2:3:3*numspikes)=yy(1:3:3*numspikes)+1;

% Scale the time axis to minutes
xx(1:3:3*numspikes)=reltimes/60;
xx(2:3:3*numspikes)=reltimes/60;
xl=[0,trial_length/60];
end
