

%%


TT = readtable("spike_data_with_mni_reduced");
fs = 256;
field= 'all_spike_times';

n = height(TT);
spike_times_regions = cell(n,1);

for i = 1:n
    all_t = [];
        s = TT.(field){i};
        if ~isempty(s)
            t = str2num(s); 
            all_t = [all_t, t];
        end
   spike_times_channels{i} = all_t;

end


% Binary spike trains
nonpadded_binary_trains_regions = cell(n,1);
for i = 1:n
    t = spike_times_channels{i};
    if isempty(t)
        nonpadded_binary_trains_regions{i} = [];
    else
        idx = unique(round(t * fs));
        v = zeros(idx(end),1);
        v(idx) = 1;
        nonpadded_binary_trains_regions{i} = v;
    end
end


%% choose region


% choose region from 1 to number n
region = 1;  
maxLen = length(nonpadded_binary_trains_regions{region});

binary_trains_regions = false(n, maxLen);

% padding operation so that all other regions have the same size binary
% spike train to train a specific region
for i = 1:n
    binary_train_region = nonpadded_binary_trains_regions{i};
    L = length(binary_train_region);

    lengt_current = min(L, maxLen);
    binary_trains_regions(i, 1:lengt_current) = binary_train_region(1:lengt_current);
end

%% turn into time-lagged matrix 


% first bin the binary spike train

binSize = 5/fs; % this is approx. 20 ms if we want to set a specific amount of ms 
% (desired_ms*256)/1000 here desired_ms = 20


% this is just to get correct dims the increments for the first 2
% dimensions of next variable containing the increments pe region
[increments_region, ~] =spikeIncrements(binary_trains_regions(region,:),binSize,fs);


% increments for all regions based on the binsize - 1 spike occured in bin , else 0
increments_all_regions = zeros(length(increments_region),n);
for i=1:n
    [increments_region, ~] =spikeIncrements(binary_trains_regions(i,:),binSize,fs);
    increments_all_regions(:,i) = increments_region;

end


 


%% Creating input matrices for a region model 

% setting self and cross coupling past history (now 0.2 and 0.1 secs
% respecitvely)
pastTimeSelf = 0.2;  
pastTimeCross = 0.1;  

numDelays_self  = floor(pastTimeSelf/binSize);
numDelays_cross = floor(pastTimeCross/binSize);



y_region = increments_all_regions(:, region);      
self_history = historicalVec(numDelays_self, y_region);  



cross_history = [];

for i = 1:n
    if i == region
        continue;                                  
    end
    increment_curr = increments_all_regions(:, i);             
    curr_region_history = historicalVec(numDelays_cross, increment_curr);      
    cross_history = [cross_history, curr_region_history];           % concatenating regions
end


X = [self_history, cross_history];   % getting self and cross history
y = y_region;                       




