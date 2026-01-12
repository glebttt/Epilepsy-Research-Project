%% Linear granger causality function

function [GC_matrix, F_matrix, p_matrix] = compute_linear_granger_causality(X, maxLag)
% COMPUTE_GRANGER Compute pairwise Granger causality
% 
% Inputs:
%   X       - n_timepoints x n_channels matrix (rows = time, cols = channels)
%   maxLag  - number of past time points (VAR order) 
%
% Outputs:
%   GC_matrix - binary matrix: 1 if significant Granger causality, 0 otherwise
%   F_matrix  - F-statistics for all pairs
%   p_matrix  - p-values for all pairs
%
    [T, N] = size(X);
    if T <= 2
        error('Too few timepoints (T=%d). Need more data.', T);
    end

    %% Check stationarity
    isStat = false(1,N);
    for ch = 1:N
        try
            [tstat, ~] = adf(X(:,ch), 0, 0);
            isStat(ch) = (tstat < -2.86); % approx 5% threshold
        catch
            isStat(ch) = false;
        end
    end

    if any(~isStat)
        warning('Some channels appear non-stationary. VAR fitting may be invalid.');
    end

    %% Choose VAR order p: AIC metric
    Xv = X'; % N x T
    % reshape for trials: variables x time x trials
    Xv = reshape(Xv, N, T, 1);
    
    regmode = 'OLS';
    verb = false;
    
    [aic, ~, moaic, ~] = tsdata_to_infocrit(Xv, maxLag, regmode, verb);
    p_opt = moaic; % optimal p value

    % Penalizes too large models (many parameters)
    % Penalizes too small models (poor fit)
    % The optimal p minimizes AIC.

    %% Estimate VAR model
    try
        [A, SIG] = tsdata_to_var(Xv, p_opt, 'OLS');
        % A = VAR coefficient matrices, size N×N×p
        % SIG = residual covariance matrix, size N×N
    catch ME
        error('VAR model estimation failed: %s', ME.message);
    end
    %% Check stability 
    try % ensures all eigenvalues of the VAR
        % stab = var_check(A, SIG); % DEPRECATED
        [n,~,p] = size(A);
        A1 = [reshape(A,n,n*p); eye(n*(p-1)), zeros(n*(p-1),n)];
        
        % Spectral radius
        rho = max(abs(eig(A1)));
        
        if rho >= 1
            warning('VAR model is unstable: spectral radius = %f', rho);
        else
            disp('VAR model is stable');
        end
    catch ME
        error('Stability check failed: %s', ME.message);
    end
    
    %% Granger F-statistic matrix and p-values
    try
        F_matrix = zeros(N,N);
        p_matrix = zeros(N,N);
        
        for i = 1:N
            for j = 1:N
                if i == j % NO make sense self-coupling
                    F_matrix(i,j) = 0;
                    p_matrix(i,j) = 0;
                    continue;
                end
                x = i;    % source
                y = j;    % target
                [F_matrix(i,j), p_matrix(i,j)] = var_to_mvgc(A, SIG, x, y, Xv, 'OLS', 'F');
            end
        end
    catch ME
        error('Error computing Granger causality F-statistics: %s', ME.message);
    end
    
    % 1 = significant Granger causality
    % 0 = not significant
    alpha = 0.05;
    GC_matrix = (p_matrix < alpha) & ~eye(N); % Avoid to count self-coupling as significant
    % Assuming a significance level of 0.05, if the p-value is lesser than 0.05, 
    % then we do NOT reject the null hypothesis that X does NOT granger cause Y
end