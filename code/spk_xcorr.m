
% Fast calculation of auto- and cross-correlations for spike trains
%
% Inputs...
%   events1 : time stamps for 1st signal
%   events2 : time stamps for 2nd signal
%   interval_min and interval_max: limits on the interval between events
%   num_bins : number of bins
%
% Outputs...
%   correlogram : histogram of intervals with edges = linspace(interval_min,interval_max,num_bins)
%   event_pairs : matrix of individual events with columns [interval n m]
%   where n and m are the indices of events for events1 and events2,
%   respectively

function [correlogram,event_pairs] = spk_xcorr(events1,events2,interval_min,interval_max,num_bins)

n1 = length(events1);
n2 = length(events2);

if n1 == 0 || n2 == 0
    correlogram = [];
    event_pairs = [];
    return;
end

% preallocate using a rough estimate of # of time difference required (assuming independence)
maxTTT = (interval_max-interval_min);
eN = ceil((max(n1, n2))^2 * maxTTT * 2 / max(events1(end), events1(end)));
event_pairs = zeros(10 * eN, 3);

% compute all the time differences
lastStartIdx = 1;
k = 1;
for n = 1:n1
    incIdx = 0;
    for m = lastStartIdx:n2
        timeDiff = events2(m) - events1(n);
        if timeDiff >= interval_min
            if incIdx==0
                incIdx = m;
            end
            if timeDiff <= interval_max
                event_pairs(k,:) = [timeDiff n m];
                k = k + 1;
            else % this is the ending point
                break;
            end
        end
    end
    if incIdx>0
        lastStartIdx = incIdx;
    end
end
event_pairs = event_pairs(1:(k-1),:);
edges = linspace(interval_min,interval_max,num_bins+1);
correlogram = histcounts(event_pairs(:,1),edges);
correlogram = correlogram';