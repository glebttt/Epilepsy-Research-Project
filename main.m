%% Main script

clear; close all; clc;
% Change the path so ALL the folders are included
addpath(genpath('C:\Users\Usuario\Desktop\GITHUB\Epilepsy-Research-Project'));

%% Preprocessing data
file_path = "C:\Users\Usuario\Desktop\GITHUB\Epilepsy-Research-Project\Connectivity_analysis_source_data\P1\P1_S6\spike_data_with_mni_reduced.csv";
%file_path = readtable("softpriors_counts_and_times_per_channel.xlsx",'Sheet','filtered_datafeature_based_o_5');

fs = 256; % 256 samples every 1 second
% bin_size = 0.02; % 20 ms
bin_size = 1; % increase bin size for now to reduce computational complexity
plot_flag = 0; % 0: skip visualizations 1: visualize

[binned_binary_trains, padded_binary_trains, total_num_spike] = preprocess_spike_data(file_path, fs, bin_size, plot_flag);

%% Granger Causality
X = binned_binary_trains.'; % input GC matrix

% Provisional input matrix (too avoid computational complexity)
% X_provisional = binned_binary_trains(1:12, :); % example to practice with
% X = X_provisional.';

maxLag = 2;  % number of past time points to include
% Too small maxLag: you may miss slow causal effects; the model may underfit.
% Too large maxLag: increases dimensionality → more parameters → overfitting, more variance, less reliable statistics.

[GC_matrix, F_matrix, p_matrix] = compute_granger_causality(X, maxLag);

%% Network from GC
% No self-edges
GC_matrix_noSelf = GC_matrix;
N = size(GC_matrix,1);
for i = 1:N
    GC_matrix_noSelf(i,i) = 0;
end

G = digraph(GC_matrix_noSelf);
nodeLabels = { ...
    'LA1-LA6', 'LA6-LA9', 'LA9-LA12', ...
    'LH1-LH6', 'LH6-LH9', 'LH9-LH12', ...
    'LP1-LP6', 'LP6-LP9', 'LP9-LP12', ...
    'LS1-LS6', 'LS6-LS9', ...
    'LT1-LT3', 'LT3-LT5', ...
    'RA1-RA6', 'RA6-RA10', 'RA10-RA14', ...
    'RH1-RH6', 'RH6-RH10', 'RH10-RH14', ...
    'RP1-RP6', 'RP6-RP9', 'RP9-RP12', ...
    'RS1-RS6', 'RS6-RS9', ...
    'RT1-RT3', 'RT3-RT6'};
G.Nodes.Name = nodeLabels(:);

figure;
h = plot(G,'Layout','subspace3');
h.NodeLabel = G.Nodes.Name;
title('Granger Causality Network - Epilepsy Channels');

%% Poisson Point Process GLM


%% 3P network
