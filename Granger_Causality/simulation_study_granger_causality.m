%% ----------------------------------------------------------
% Simulation Study (Granger Causality: Linear and Nonlinear)
% -----------------------------------------------------------
clear; close all; clc; rng(1);

% Change the path so ALL the folders are included
addpath(genpath(pwd));

%% ----------------------------
% 1. Synthetic Data Generation
% -----------------------------
plot_flag = 0; 
N  = 7;                 % number of brain regions
dt = 0.004;             % time bin (4 ms)
Tsec = 600;             % total duration (seconds)


[Y, A] = generate_synthetic_data(N, dt, Tsec, plot_flag);

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
sigma_ms = 5;                  
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
% for computation, truncate it to where it's effectively zero 
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
% now it's meaningful in spikes per second 

% conv: convolution spreads each spike across nearby bins, making the
% signal smooth instead of spikes at exact points

% the 'same' argument keeps the output the same length as the original
% spike train 

%% ---------------------------
% 3. Linear Granger Causality
% ----------------------------

% TODO: INCREASE PRECISION (0.44) FOR SYNTHETIC DATA + TEST IT W/ CLINICAL
% DATA

X = Y_smooth'; % input GC matrix

maxLag = 5;  % number of past time points to include
% Too small maxLag: you may miss slow causal effects; the model may underfit.
% Too large maxLag: increases dimensionality → more parameters → overfitting, more variance, less reliable statistics.

[GC_matrix, F_matrix, p_matrix] = compute_linear_granger_causality(X, maxLag);


% METRICS 
TP = sum((GC_matrix==1) & (A==1), 'all');
FP = sum((GC_matrix==1) & (A==0), 'all');
FN = sum((GC_matrix==0) & (A==1), 'all');
TN = sum((GC_matrix==0) & (A==0), 'all');

fprintf('TP=%d, FP=%d, FN=%d, TN=%d\n', TP, FP, FN, TN);
precision = TP / (TP + FP);
recall = TP / (TP + FN);

fprintf('Precision=%.2f, Recall=%.2f\n', precision, recall);

%% -------------------------------
% 4. Nonlinear Granger Causality
% --------------------------------

% TODO: IMPLEMENT AUXILIARY FUNCTIONS + TEST IT BOTH SYNTHETIC & CLINICAL
% DATA

maxLag = 5;  % number of past time points to include
Q = 8;  % quantization levels (typical: 4–10)
nSurr  = 40;    % number of surrogates
alpha  = 0.05;  % significance level

results_nl = compute_nonLinear_granger_causality(X, maxLag, Q, nSurr, alpha);

GC_nonlinear = results_nl.sig;   % binary significance matrix
C_strength   = results_nl.C;     % coupling strength matrix

disp('Nonlinear GC matrix (significant):');
disp(GC_nonlinear);
disp('Nonlinear GC strength matrix:');
disp(C_strength);

% METRICS
TP_nl = sum((GC_nonlinear==1) & (A==1), 'all');
FP_nl = sum((GC_nonlinear==1) & (A==0), 'all');
FN_nl = sum((GC_nonlinear==0) & (A==1), 'all');
TN_nl = sum((GC_nonlinear==0) & (A==0), 'all');
fprintf('NONLINEAR GC: TP=%d, FP=%d, FN=%d, TN=%d\n', TP_nl, FP_nl, FN_nl, TN_nl);

precision_nl = TP_nl / (TP_nl + FP_nl + eps);
recall_nl    = TP_nl / (TP_nl + FN_nl + eps);
fprintf('NONLINEAR GC: Precision=%.2f, Recall=%.2f\n', precision_nl, recall_nl);

