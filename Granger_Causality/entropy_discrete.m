function H = entropy_discrete(X)
% ENTROPY_DISCRETE  Entropy of discrete multivariate variable
%
% Input:
%   X : N x D matrix of integers (symbols)
%
% Output:
%   H : entropy in nats

    if isempty(X)
        H = 0;
        return;
    end

    % convert rows to unique symbols
    [~,~,idx] = unique(X, 'rows');
    counts = accumarray(idx,1);

    p = counts / sum(counts);
    p(p==0) = [];

    H = -sum(p .* log(p));
end
