function Hc = conditional_entropy_corrected(xi, Xk)
% CORRECTED conditional entropy Hc(xi | Xk)
%
% Inputs:
%   xi : N x 1 (discrete)
%   Xk : N x K (discrete embedding)
%
% Output:
%   Hc : corrected conditional entropy

    % H(xi | Xk)
    H_joint = entropy_discrete([xi Xk]);
    H_cond  = entropy_discrete(Xk);

    H = H_joint - H_cond;

    % bias correction
    if isempty(Xk)
        n_unique = 1;
    else
        [~,~,idx] = unique(Xk, 'rows');
        counts = accumarray(idx,1);
        n_unique = sum(counts == 1) / length(xi);
    end

    H_xi = entropy_discrete(xi);

    Hc = H + n_unique * H_xi;
end
