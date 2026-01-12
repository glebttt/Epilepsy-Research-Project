%% ----------------------------------------------------------
% Simulation Study (Granger Causality: Linear and Nonlinear)
% -----------------------------------------------------------



%% ----------------------------
% 1. Synthetic Data Generation
% -----------------------------

% clear the workspace 
% clear the command window 
% fix a random seed for reproducible simulation 
% (ensures the same spike trains every run)

clear; clc; rng(1);

%% -------------------------------
% Simulation parameters
% -------------------------------
N  = 7;                 % number of brain regions
dt = 0.004;             % time bin (4 ms)
Tsec = 600;             % total duration (seconds)
T  = round(Tsec / dt);  % total number of time bins

% this is a discrete-time approximation of a continuous point process
% bin size is small enough that at most one spike per bin is reasonable 

Kh = 25;                % self-history length (100 ms) (25 bins)
Kc = 25;                % coupling history length (100 ms) (25 bins)

% when deciding whether a neuron spikes now, the model looks back 100 ms
% into the past
% it does not care about spikes older than 100ms 

% self-history: if a neuron spiked recently
% that affects how likely it is to fire again shortly after 

% coupling history: if another neuron fired recently, it might influence me
% now
% this models communication between neurons or brain regions 

% these self and coupling histories capture:

% 1. refractory periods: 
% After a neuron fires, for a few milliseconds it becomes temporarily less
% excitable i.e. less likely to fire 
% "I just spoke, give me a moment before I speak again"

% 2. short-latency synaptic effects: 
% synapse: the connection where one neuron influences another 
% latency: means delay 
% short-latency synaptic effects means that when one neuron fires, it
% briefly and quickly increases the chance that another neuron will fire,
% with a small delay of a few milliseconds 


fprintf('Simulating %d regions for %.1f seconds (%d bins)\n', ...
        N, Tsec, T);

%% -------------------------------
% Baseline firing rates (Hz)
% -------------------------------

baseline_rate = 10 * ones(N,1);      % moderate firing

% N = number of regions 
% ones(N,1) --> column vector of ones, size Nx1
% multiply by 10 --> every neuron has 10 Hz baseline firing rate
% so each neuron wants to fire around 10 times per second if nothing else
% is happening 

beta0 = log(baseline_rate * dt);     % baseline log-rate

% dt are the 4 ms time bins 
% in each bin, a neuron can fire 0 or 1 spikes 
% baseline_rate * dt because we want to know how many spikes do I expect in
% each 4 ms bin 
% on average, in a single 4 ms bin, the neuron will fire 0.04 spikes 
% since 0.04 < 1, most bins are empty (no spike)
% occasionally, a spike occurs, roughly 10 spikes per second, matching the
% target 10 Hz
% why take the log? 
% the neurons firing follows a Poisson process in discrete time
% this is an assumption of the Poisson GLM
% baseline term: when there is no influence from self or other neurons, the
% neuron should fire at its baseline rate 
% a little bit of maths to show that beta0 is the log :)



%% -------------------------------
% Network adjacency matrix (GROUND TRUTH)
% A(j,i) = 1 means j -> i
% -------------------------------
A = zeros(N,N);

A(1,2) = 1;
A(1,3) = 1;
A(1,4) = 1;
A(2,5) = 1;
A(3,5) = 1;
A(4,6) = 1;
A(5,7) = 1;
A(6,7) = 1;


% each 1 indicates that neuron j influences neuron i via the coupling
% filter C(j, i, :)
% each 0 means no direct influence 
% the matrix encodes which neurons are connected and in which direction 
% I am generating the data, so I already know the true network connections



%% -------------------------------
% Self-history filters (refractory)
% -------------------------------
Ah = 2.0;        % refractory strength
tau_h = 5;       % decay (bins)

h = zeros(N,Kh);
for i = 1:N
    h(i,:) = -Ah * exp(-(1:Kh)/tau_h);
end

% 1. FILTER INTUITION
% a filter is a set of numbers that says how past spikes influence current
% firing 
% recent spikes might suppress or increase future spikes 

% modelling refractoriness: after a neuron spikes, it is less likely to
% fire for a short while 

% so the filter will be negative, strongest immediately after a spike, and
% then decays towards zero


% 2. CODE EXPLANATION
% Ah is the strength of the refractory effect
% bigger number --> neuron is more strongly suppressed right after a spike

% tau_h: this is how fast the effect decays over time (in number of bins)
% larger tau --> slower decay (suppression lasts longer)
% smaller tau --> shorter refractory period 

% h is a matrix: N regions x Kh time bins 
% each row will store the self-history filter for that region 
% Kh = 25 (each region remembers 25 past bins ~ 100 ms)

% the loop:
% for each row (1:N):
% 1:Kh --> row vector [1, 2, 3, ..., Kh]
% exp(...) computes the decaying exponential for each lag 
% multiplying by -Ah gives the self-history filter values for each lag 
% result: h(i,:) is a row vector of length Kh representing the refractory
% effect at each past time step 
% i.e. each element in h(i, :) tells you how strongly a spike k steps ago
% suppresses the current spike 
% right after a spike: strong suppression
% later: weak suppression 
% after ~100 ms, almost no effect

% 3. HOW IT IS USED IN SIMULATION
% for computing eta, multiplying each past spike by the filter value h(i,k) 
% (for region i  and time step k)
% adds a negative contribution --> reduces current firing probability
% this implements the refractory period!!! :)

%% -------------------------------
% Coupling filters (propagation)
% -------------------------------
B = 1.5;         % coupling strength
d = 5;           % delay (bins = 20 ms)
sigma = 2;       % smoothness

% coupling filter: how neuron j affects neuron i 
% B: how strong the effect is 
% d: how long it takes for the effect to happen
% sigma: how spread out the effect is over time 

C = zeros(N,N,Kc);
% C is a 3D matrix 
% C(j, i, k) = effect of neuron j on neuron i at lag k
% C will store how each presynaptic neuron affects each postsynaptic neuron
% over the past 100 ms
% presynaptic: the neuron that sends a signal (spike)
% postsynaptic: the neuron that receives the signal 

kvec = 1:Kc;
% vector [1, 2, ... , Kc] (time lags), where Kc is the coupling history
% length 



for j = 1:N
    for i = 1:N
        if A(j,i) == 1
            C(j,i,:) = B * exp(-(kvec - d).^2 / (2*sigma^2));
        end
    end
end

% loop explanation:
% loop over all neuron pairs (j, i)
% only create a filter if A(j:i) == 1
% A is the adjacency matrix (network connectivity)
% So we only model connections that exist 
% exp(...) is the Gaussian curve over lags 
% therefore, spike from neuron j increases probability that neuron i spikes
% after a short delay
% the effect is strongest at lag d, decays around it like a bell curve
% makes sense to use Gaussian since influence rises as you approach lag = d
% peaks at lag = d
% slowly decays after lag = d 


% HOW IT IS USED IN SIMULATION
% for each presynaptic region j, each past spike Y(j, t-k):
% multiply by the filter weight C(j, i, k)
% adds to the postsynaptic regions linear predictor eta
% this models spike propagation over time



%% -------------------------------
% Initialise spike trains
% -------------------------------
Y = zeros(N,T);


% preallocaes spike matrix 
% rows = regions 
% columns = time bins 

%% -------------------------------
% Recursive simulation
% -------------------------------


t_start = max(Kh, Kc) + 1;
% model requires up to Kh past bins for self-history
% model required up to Kc past bins for coupling 
% e.g. at t=1 there is no past, t=10 u have past up til t-1 = 9

% time
for t = t_start:T
    % only start simulating spikes once all required history terms exist
    % region
    for i = 1:N
        % simulating spike at time t for each region 
        
        % Linear predictor
        eta = beta0(i);
        
        % Self-history contribution
        % history 
        for k = 1:Kh
            % filter times the spike at lag k
            % Y(i, t-k) means did region i spike exactly k time bins ago
            eta = eta + h(i,k) * Y(i, t-k);
        end
        
        % Coupling contribution
        % region
        for j = 1:N
            % if j not equal to i 
            if j ~= i && A(j,i) == 1
                % history
                for k = 1:Kc
                    eta = eta + C(j,i,k) * Y(j, t-k);
                end
            end
        end
        
        % Conditional intensity
        lambda = exp(eta);
        
        % Poisson spike generation
        spike = poissrnd(lambda);
        
        % Enforce at most one spike per bin
        Y(i,t) = min(spike, 1);
    end
end

% NOTE:
% time bins: t = 1, ..., 25 are never simulated
% t_start = 26 
% such a small part of the duration that its negligible 

fprintf('Simulation complete.\n');

%% -------------------------------
% Basic sanity checks
% -------------------------------
firing_rates = sum(Y,2) / Tsec;

fprintf('\nMean firing rates (Hz):\n');
for i = 1:N
    fprintf('Region %d: %.2f Hz\n', i, firing_rates(i));
end


% expected rate is 10 Hz so they should be around that 
% depending on which regions affect which
% e.g. in region 1 no incoming excitation --> only self-refractory effect
% rate slightly below baseline 

% e.g. in region 7 excitation from region 5 and region 6
% rate above baseline since 5 and 6 already have elevated rates 
% cascade amplification effect: which means spikes propogate through a
% directed network, and each stage increases the firing rate of the next,
% even if individual connections are not very strong 
% meaning mean firing rates increase downstream due to accumulated excitation

% purpose of basic sanity checks:
% ensures baseline rates are correct
% model is stable (no runaway excitation)


%% --------------------------------
% plotting the binary spike trains
% ---------------------------------

fs = 1/dt;
time = (0:T-1) / fs;

% using zoomed windows for visualisation
tmin = 0;   % seconds
tmax = 600;

idx = time >= tmin & time <= tmax;

for i = 1:5
    figure;
    plot(time(idx), Y(i,idx), '-o', ...
         'Color', 'b', ...
         'MarkerSize', 4, ...
         'LineWidth', 1);
    xlabel('Time (s)');
    ylabel('Spike response');
    title(sprintf('Region %d Spike Train', i));
    ylim([-0.1 1.1]);
    grid on;
end

% plots binary spike trains, one figure per region 


%% ----------
% to do
% -----------

% to incorporate clinical data: 
% define beta0 using mean firing rates from clinical data
% filter lengths (self and cross-coupling)
% Ah and B values 
% A (network design)












%% ------------------------------------
% 2. Smoothing the Binary Spike Trains
% -------------------------------------

% Binary spike trains are 0/1 sequences
% this is not ideal for classical Granger causality (GC) methods
% because GC assumes continuous-valued, roughly Gaussian, stationary
% signals

% Smoothing transforms the spike train into an instantaneous firing rate
% estimate 



% kernel width in milliseconds
sigma_ms = 20;                  
% 20 ms is a typical choice because:
% it is long enough to smooth out single-spike variability
% it is short enough to preserve temporal resolution for interactions
% between neurons

% convert ms to time bins
sigma = sigma_ms / 1000 / dt;   
% dividing by 1000 converts milliseconds to seconds 
% dividing by dt = 0.004 converts seconds to number of time bins in your
% discrete-time spike train 
% this tells you how many bins your Gaussian kernel has


L = ceil(5*sigma);
% Gaussian theoretically extends from - to + infinity
% for computation, truncate it to where its effectively zero 
% L: 5*sigma is a common rule of thumb 
% ceil just rounds up to the nearest integer
x = -L:L;
% the kernel spans from -L to L 
g = exp(-x.^2/(2*sigma^2));

g = g / sum(g);
% normalise so that total area = 1 (since we truncated)
% ensures the smoothed firing rate has the same mean scale as original
% spikes

Y_smooth = zeros(size(Y));
% empty matrix to store the smoothed values 
% Y_smooth will be in Hz

for i = 1:N
    Y_smooth(i,:) = conv(Y(i,:), g, 'same') / dt;
end

% loop explanation
% Y(i,:) is a binary array = 0,1,0,0,...
% each element is 1 if the neuron spiked in that bin, 0 otherwise
% bin size is dt = 0.004s (4 ms)
% dividing by dt converts the smoothed spike counts per bin to firing rate
% in Hz 
% now its meaningful in spikes per second 

% conv: convolution spreads each spike across nearby bins, making the
% signal smooth instead of spikes at exact points

% the 'same' argument keeps the output the same length as the original
% spike train 



%% ----------------------------------
% plotting the smoothed spike train
% -----------------------------------

region = 1;   % choose region to plot

figure;
plot(time(idx), Y_smooth(region,idx), 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Firing rate (Hz)');
title(sprintf('Region %d Smoothed Firing Rate', region));
grid on;











%% ---------------------------
% 3. Linear Granger Causality
% ----------------------------



%% --------------
% MVGC toolbox
% ---------------

addpath(genpath('/Users/lisakelly/Desktop/masters/y1/project_sem_1/synthetic/linear GC/github_repo'));

maxLag = 1;

X = Y_smooth';
% transposing since compute_granger_causality_Lisa expects n_timepoints x
% n_channels 

[GC_matrix, F_matrix, p_matrix] = compute_granger_causality_Lisa(X, maxLag);





%% ------------------
% downsampling
% -------------------

ds = 100;  % downsample factor (4 ms → 40 ms)
X_ds = Y_smooth(:, 1:ds:end)';
X_ds = double(X_ds);


maxLag = 1;
[GC_matrix, F_matrix, p_matrix] = compute_granger_causality(X_ds, maxLag);









%% -------------------------------
% 4. Nonlinear Granger Causality
% --------------------------------

