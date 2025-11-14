% this code was provided by Professor Anirudh Wodeyar
% it creates a time-lagged Hankel matrix with the specified amount of 
% past bins to use


function vec_Hist = historicalVec(numDelays,vec)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

vec_Hist= zeros(length(vec),numDelays);
for i = 1:numDelays
        vec_Hist(i+1:end,i) = vec(1:end-i);
end

vec_Hist = vec_Hist(:,end:-1:1, :); %reverse so that when in design mat, the oldest in time is farthest back
vec_Hist = reshape(vec_Hist, size(vec_Hist,1), size(vec_Hist,2)); %reshape so that this is a design matrix again



% % this code was provided by Professor Anirudh Wodeyar
% % it creates a time-lagged Hankel matrix with the specified amount of 
% % past bins to use
% 
% 
% function vec_Hist = historicalVec(numDelays,regions,vec)
% %UNTITLED Summary of this function goes here
% %   Detailed explanation goes here
% 
% vec_Hist= zeros(length(vec),numDelays , regions);
% 
% for r = 1:regions% go over all oscillators
%     for i = 1:numDelays
%         vec_Hist(i+1:end,i, r) = vec(1:end-i,r);
%     end
% end
% vec_Hist = vec_Hist(:,end:-1:1, :); %reverse so that when in design mat, the oldest in time is farthest back
% vec_Hist = reshape(vec_Hist, size(vec_Hist,1), size(vec_Hist,2) * size(vec_Hist,3)); %reshape so that this is a design matrix again
