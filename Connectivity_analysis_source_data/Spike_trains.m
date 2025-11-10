

TT = readtable("softpriors_counts_and_times_per_channel.xlsx",'Sheet','filtered_datafeature_based_o_5');
fs = 256;
fields = {'p_spike_1_0_times','p_poly_1_0_times','p_sw_1_0_times'};

n = height(TT);
spike_times_channels = cell(n,1);

for i = 1:n
    all_t = [];
    for f = 1:numel(fields)
        s = TT.(fields{f}){i};
        if ~isempty(s)
            t = str2num(s); 
            all_t = [all_t, t];
        end
    end
    spike_times_channels{i} = all_t;

end

% Binary spike trains
binary_trains = cell(n,1);
for i = 1:n
    t = spike_times_channels{i};
    if isempty(t)
        binary_trains{i} = [];
    else
        idx = unique(round(t * fs));
        v = zeros(idx(end),1);
        v(idx) = 1;
        binary_trains{i} = v;
    end
end

% for i = 1:5
%     figure;
%     time = (0:length(binary_trains{i})-1)/fs;
%     stem(time, binary_trains{i});
%     xlabel('Time (s)');
%     ylabel('Spike response');
%     title(sprintf('Channel %d', i));
% end




%% This code is for combining all channels into a single spike train (to be changed to selected channels)
% finidng the longest of the binary spike trains
maxlen = 0;
for i = 1:n
    L = length(binary_trains{i});
    if L > maxlen
        maxlen = L;
    end
end


% ensuring that all binary spike trains are same size as the biggest
% channel
padded_binary_trains = false(n, maxlen);
for i = 1:n
    binary_train_channel = binary_trains{i};
    L = length(binary_train_channel);
    padded_binary_trains(i, 1:L) = binary_train_channel;   
    
end



% getting the combined spike train 
combined_train = false(maxlen,1);

% combined_train = zeros(maxlen,1);  

% 
for i = 1:n
    binary_train_channel = padded_binary_trains(i,:);
    binary_train_channel = binary_train_channel.';

    L = length(binary_train_channel);
    combined_train(1:L) = combined_train(1:L) | binary_train_channel;  


end


% Plot combined spike train
% time = (0:length(combined_train)-1)/fs;
% figure;
% stem(time, combined_train);
% xlabel('Time (s)');
% ylabel('Combined spike response');
% title('Combined multi-channel spike train');


%%

total_sum = 0;
for i=1:length(binary_trains)
    total_sum = total_sum + sum(binary_trains{i});
end
