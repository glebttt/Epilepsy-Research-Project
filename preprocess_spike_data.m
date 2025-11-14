%% Preprocessing interictal spikes data function

function [binned_binary_trains, padded_binary_trains, total_num_spike] = preprocess_spike_data(file_path, fs, bin_size, plot_flag)
% Reads spike times from a CSV and converts them into binned & padded binary spike trains.
%
% Inputs:
%   file_path : path to .csv file
%   fs        : sampling frequency (Hz)
%   bin_size  : bin size in seconds (e.g., 0.02 for 20 ms)
%   plot_flag : 1 = plot, 0 = NO plots
%
% Outputs:
%   binned_binary_trains : [n_channels × n_bins] matrix
%   padded_binary_trains : [n_channels × maxlen] logical matrix
%   total_num_spike      
    %% Load data
    TT = readtable(file_path);
    n = height(TT);
    spike_times_channels = cell(n,1);

    %% Extract spike times from each channel
    for i = 1:n
        s = TT.all_spike_times{i};
        if ~isempty(s)
            t = str2num(s); 
            spike_times_channels{i} = t;
        end
    end

    %% Convert to binary spike trains
    binary_trains = cell(n,1);
    for i = 1:n
        t = spike_times_channels{i};
        if isempty(t)
            binary_trains{i} = [];
        else
            idx = unique(round(t * fs)); % returns the output sorted in ascending order by default
            v = zeros(idx(end),1);
            v(idx) = 1;
            binary_trains{i} = v;
        end
    end

    if plot_flag
        for i = 1:5 % first 5 channels
            figure;
            time = (0:length(binary_trains{i})-1)/fs;
            stem(time, binary_trains{i});
            xlabel('Time (s)');
            ylabel('Spike response');
            title(sprintf('Channel %d', i));
        end
    end

    %% Padding: all binary trains to same length
    % Ensuring that all binary spike trains are same size as the biggest
    % channel
   
    % maxlen = 0;
    % for i = 1:n
    %     L = length(binary_trains{i});
    %     if L > maxlen
    %         maxlen = L;
    %     end
    % end
    maxlen = max(cellfun(@length,binary_trains));
    padded_binary_trains = false(n, maxlen);
    for i = 1:n
        L = length(binary_trains{i});
        if L > 0
            padded_binary_trains(i,1:L) = binary_trains{i};
        end
    end

    %% Binning
    bin_samples = round(bin_size * fs);
    [n_channels, n_timepoints] = size(padded_binary_trains);
    n_bins = floor(n_timepoints / bin_samples); % ensures we only take complete bins
    binned_binary_trains = zeros(n_channels, n_bins); % preallocate the output matrix

    for ch = 1:n_channels
        for b = 1:n_bins
            idx_start = (b-1)*bin_samples + 1;
            idx_end = b*bin_samples;
            % check if there is at least one spike in this bin
            binned_binary_trains(ch,b) = any(padded_binary_trains(ch, idx_start:idx_end));
        end
    end
    %% Obtain the total nº of spikes 
    total_num_spike = 0;
    for i=1:length(binary_trains)
        total_num_spike = total_num_spike + sum(binary_trains{i});
    end
end


% % Getting the combined spike train 
% combined_train = false(maxlen,1);
% % combined_train = zeros(maxlen,1);  
% 
% for i = 1:n
%     binary_train_channel = padded_binary_trains(i,:);
%     binary_train_channel = binary_train_channel.';
% 
%     L = length(binary_train_channel);
%     combined_train(1:L) = combined_train(1:L) | binary_train_channel;  
% end
% 
% % Plot combined spike train
% if plot_flag
%     time = (0:length(combined_train)-1)/fs;
%     figure;
%     stem(time, combined_train);
%     xlabel('Time (s)');
%     ylabel('Combined spike response');
%     title('Combined multi-channel spike train');
% end
