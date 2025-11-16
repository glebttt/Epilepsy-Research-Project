

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




%% %% __________ FIT GLM (Poisson) _________________________
fprintf('Fitting Poisson GLM...\n');


[b, dev, stats] = glmfit(X, y, 'poisson','link','log');


% Extract intercept and weigths beta(0) and beta(1)
const = b(1);
weights = b(2:end);




numCoupledNeurons = size(X,2) - numDelays_self;  
numCoupledNeurons = numCoupledNeurons / numDelays_cross; 

neuronNames = cell(1,numCoupledNeurons);
for k = 1:numCoupledNeurons
    neuronNames{k} = ['Neuron ' num2str(k)];
end



selfFilt = weights(1 : numDelays_self);


couplingFilts = [];
idx = numDelays_self + 1;

for n = 1:numCoupledNeurons
    couplingFilts(:, n) = weights(idx : idx + numDelays_cross - 1);
    idx = idx + numDelays_cross;
end

%  PREVISION FIRING RATE 
lambda = exp(const + X * weights);

fprintf("lambda : %d",lambda);

figure; clf;

% --- 1. Plot spike-history of the targeted neuron
subplot(2,2,1);
stem( - (1:numDelays_self), selfFilt, 'LineWidth', 2 );
title('Spike-history filter'); xlabel('lag (bins)'); ylabel('weight');
grid on;

% --- 2. Plot coupling filters 
subplot(2,2,2);
hold on;
for n = 1:numCoupledNeurons
    plot( - (1:numDelays_cross), couplingFilts(:,n), 'LineWidth', 2 );
end
hold off;
title('Coupling filters (incoming)'); 
xlabel('lag (bins)'); ylabel('weight');
legend(neuronNames, 'Location', 'best');  % se hai i nomi
grid on;

% --- 3. Plot spike vs firing rate
subplot(2,2,3);
binsToPlot = 1:200;   % cambia come vuoi
plot(binsToPlot, y(binsToPlot), 'k', 'LineWidth', 1.5); hold on;
plot(binsToPlot, lambda(binsToPlot), 'r', 'LineWidth', 2);
title('Spikes vs Predicted Firing Rate');
xlabel('time bin'); ylabel('spikes / rate');
legend('Observed', 'Predicted');
grid on;

% --- 4. Deviance / LL / AIC
subplot(2,2,4);
text(0.1, 0.8, sprintf('Deviance: %.2f', dev), 'FontSize', 12);
logL = sum( y .* log(lambda + eps) - lambda );
AIC = 2*length(b) - 2*logL;
text(0.1, 0.6, sprintf('Log-Likelihood: %.2f', logL), 'FontSize', 12);
text(0.1, 0.4, sprintf('AIC: %.2f', AIC), 'FontSize', 12);
axis off;
title('Model Statistics');

sgtitle(sprintf('GLM neuron %d', region));



