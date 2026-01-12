function Xq = quantize_uniform(X, Q)
% QUANTIZE_UNIFORM  Uniform quantization of time series
%
% Inputs:
%   X : T x M matrix (time x variables)
%   Q : number of quantization levels
%
% Output:
%   Xq : T x M matrix with values in {0,1,...,Q-1}

    % normalize each channel: mean 0, var 1
    Xn = zscore(X);

    [T, M] = size(X);
    Xq = zeros(T, M);

    for m = 1:M
        xmin = min(Xn(:,m));
        xmax = max(Xn(:,m));

        if xmax == xmin
            Xq(:,m) = 0;
            continue;
        end

        r = (xmax - xmin) / Q;

        % quantization
        q = floor((Xn(:,m) - xmin) / r);

        % safety (values in [0, Q-1])
        q(q < 0)   = 0;
        q(q >= Q)  = Q-1;

        Xq(:,m) = q;
    end
end
