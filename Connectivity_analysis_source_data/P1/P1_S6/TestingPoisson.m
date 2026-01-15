
%%


% change this manually to number of nodes in network and factor_n 
n = 4;
factor_n = 20;

importFileName = sprintf('test_%d_nodes_factor_%d.mat',n,factor_n);
load(importFileName) 
% this is for one factor and one network design 



numDelays_self  = 5;
numDelays_cross = 5; 


for region=1:n
    % region = 2; % loop over this

    y_region = Y(region,:).';

    self_history = historicalVec(numDelays_self, y_region);



    cross_history = [];


    n = size(Y,1);
    for i = 1:n
        if i == region
            continue;
        end
        increment_curr = Y(i, :).';
        curr_region_history = historicalVec(numDelays_cross, increment_curr);
        cross_history = [cross_history, curr_region_history];           % concatenating regions
    end


    X = [self_history, cross_history];   % getting self and cross history
    y = y_region;


    % mdl = fitglm(X,y,'linear','Distribution','poisson','Link','log');


    % [B, FitInfo] = lassoglm(X, y, 'poisson', 'Standardize', false); % get 100 lambdas
    L = max(numDelays_self, numDelays_cross);
    X = X(L+1:end, :);
    y = y_region(L+1:end);


    otherRegions = setdiff(1:n, region, 'stable');

    group = region *ones(1, numDelays_self);   % group 1 = self

    g = region;
    for i = 1:length(otherRegions)
        group = [group, otherRegions(i) * ones(1, numDelays_cross)];

    end % debug for other regions


    % Sanity check
    assert(length(group) == size(X,2), 'Group vector length must equal ncol(X).');
    %
    % save('glm_input_4_nodes_region2_factor20.mat', 'X', 'y', 'group', ...
    %      'region', 'numDelays_self', 'numDelays_cross');


    fileName = sprintf('glm_input_4_nodes_region%d_factor20.mat',region);
    save(fileName, 'X', 'y', 'group', ...
        'region', 'numDelays_self', 'numDelays_cross');


end


% Save for R
% save('glm_input_2_node_2.mat', 'X', 'y', 'region', 'numDelays_self', 'numDelays_cross')
%% gather statistics

n = 4;
factor_n = 20;

tol = 1e-8;

% Outputs (one entry per region)
bestAIC     = nan(1,n);
bestLambda  = nan(1,n);
bestIdx     = nan(1,n);
bestDF      = nan(1,n);

AIC_all = cell(1,n);
LL_all  = cell(1,n);
df_all  = cell(1,n);
Lambda_all = cell(1,n);

% Optional: store best-fit coefficients per region
Intercept_best = cell(1,n);
B_best         = cell(1,n);


B_all       = cell(1,n);   % <-- NEW: full beta paths
Intercept_all  = cell(1,n);   % <-- NEW

for region = 1:n

    RFileName = sprintf("grpreg_%d_nodes_region%d_factor%d.mat", n, region, factor_n);
    S = load(RFileName);  % safer than load() into workspace

    % If you did NOT save X,y in the R output, you must load them from the original glm_input file:
    % (uncomment next 3 lines if needed)
    % inFile = sprintf("glm_input_%d_nodes_region%d_factor%d.mat", n, region, factor_n);
    % In = load(inFile);
    % S.X = In.X; S.y = In.y;

    X = double(S.X);
    y = double(S.y(:));
    Lambda = S.Lambda(:).';   % 1 x L
    Beta   = double(S.Beta);  % (p+1) x L

    Intercept = Beta(1, :);   % 1 x L
    B         = Beta(2:end,:);% p x L

    L = numel(Lambda);

    LL  = zeros(1, L);
    df  = zeros(1, L);
    AIC = zeros(1, L);

    for k = 1:L
        eta = Intercept(k) + X * B(:, k);   % T x 1
        mu  = exp(eta);                     % Poisson rate

        % Degrees of freedom: intercept + nonzero coefficients
        df(k) = 1 + sum(abs(B(:, k)) > tol);

        % Poisson log-likelihood (ignoring constant log(y!))
        LL(k) = y' * log(mu) - sum(mu);

        % AIC
        AIC(k) = -2 * LL(k) + 2 * df(k);
    end

    % Best by AIC
    [bestAIC(region), idx] = min(AIC);
    bestIdx(region)    = idx;
    bestLambda(region) = Lambda(idx);
    bestDF(region)     = df(idx);

    % Store curves
    AIC_all{region}    = AIC;
    LL_all{region}     = LL;
    df_all{region}     = df;
    Lambda_all{region} = Lambda;

    % Store best coefficients (optional)
    Intercept_best{region} = Intercept(idx);
    B_best{region}         = B(:, idx);


    B_all{region}      = B;        % <-- full coefficient path
    Intercept_all{region} = Intercept;   % <-- redundant but convenient


    fprintf('Region %d: best AIC=%.6f, best lambda=%.6g, df=%d\n', ...
            region, bestAIC(region), bestLambda(region), bestDF(region));
end

% % Save everything in ONE file
% outFile = sprintf("grpreg_eval_%d_nodes_factor%d_allRegions.mat", n, factor_n);
% save(outFile, ...
%      'bestAIC','bestLambda','bestIdx','bestDF', ...
%      'AIC_all','LL_all','df_all','Lambda_all', ...
%      'Intercept_best','B_best', ...
%      'n','factor_n','tol');
% 
% 




% fprintf('Number of ')



% X()

% figure;
% subplot(2,1,1)
% stem(y_region), title('Raw y')
% 
% subplot(2,1,2)
% stem(X(:,1)), title('Self lag 1 (previous point)')


% n = 4;
% factor_n = 20;
% 
% 
% 
% 
% 
% for region=1:n
% 
%     RFileName = sprintf("grpreg_%d_nodes_region%d_factor%d.mat",n, region, factor_n);
% 
% 
%     load(RFileName)
% 
% 
%     % Ensure correct shapes
%     X = double(X);
%     y = double(y(:));
%     Lambda = Lambda(:).';      % make it 1 x L
% 
%     Intercept = Beta(1, :);    % 1 x L
%     B = Beta(2:end, :);        % p x L
% 
%     L = length(Lambda);
%     tol = 1e-8;
% 
%     LL  = zeros(1, L);
%     df  = zeros(1, L);
%     AIC = zeros(1, L);
% 
%     for k = 1:L
%         eta = Intercept(k) + X * B(:, k);   % T x 1
%         mu  = exp(eta);                     % Poisson rate
% 
% 
%         % Degrees of freedom: intercept + nonzero coefficients
%         df(k) = 1 + sum(abs(B(:, k)) > tol);
% 
% 
%         ratepred_pGLM = exp(Intercept(k) + X*B(:,k)); % rate under exp nonlinearity
%         LL(k) = y'*log(ratepred_pGLM) - sum(ratepred_pGLM);
% 
% 
%         % AIC
%         AIC(k) = -2 * LL(k)  + 2 * df(k);
% 
%         % mu = exp(eta);
%         %
%         % % Bernoulli spike probability
%         % p  = 1 - exp(-mu);
%         % p  = min(max(p,1e-12), 1-1e-12);
%         %
%         % % Log-likelihood (Bernoulli)
%         % LL(k) = sum( y .* log(p) + (1-y).*log(1-p) );
%         %
%         % % Degrees of freedom (intercept + active betas)
%         % df(k) = 1 + sum(abs(B(:,k)) > tol);
%         %
%         % % AIC
%         % AIC(k) = -2 * LL(k) + 2 * df(k);
%     end
% 
%     % Find minimum AIC
%     [bestAIC, idx] = min(AIC);
%     bestLambda = Lambda(idx);
% 
%     fprintf('Minimum AIC = %.6f\n', bestAIC);
%     fprintf('Best lambda = %.6g\n', bestLambda);
%     fprintf('Degrees of freedom = %d\n', df(idx));
% 
%     % if (n == 2) B2 = B;    
%     % if (n == 3) B3 = B;    
%     % if (n == 1) B1 = B;
% 
% end


%%
figure;
plot(AIC,'-o');


% 






