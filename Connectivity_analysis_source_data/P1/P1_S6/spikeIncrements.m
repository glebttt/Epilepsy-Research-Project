% code slighlty modified for binning turns (zeros or ones)

function [incVec,t] = spikeIncrements(spikeVec,binSize,Fs)
% function to take in a spike vector and turn it into an increments vector
% based on the given binsize. 
% We pad with zeros at the end as necessary to fill out the final bin.
% spikeVec is samples x 1
% binSize is in seconds
% Fs is in Hz
% 2/24/2022
% Ani Wodeyar 

if ~iscolumn(spikeVec)
    spikeVec = spikeVec';
end
% how many incs present
numIncs = ceil(length(spikeVec)/(binSize*Fs));
t = binSize/2:binSize:numIncs*binSize;
numZerosNeeded = floor((numIncs * binSize *Fs) - length(spikeVec));
% mean over spikes present in each bin:
incVec = any(reshape([spikeVec; zeros(numZerosNeeded,1)], binSize*Fs, numIncs), 1)';

