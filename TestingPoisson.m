load('test.mat')
%%


patientNo = 1;

numDelays_self  = 5;
numDelays_cross = 5;



region = 1;
y_region = Y(region,:);

self_history = historicalVec(numDelays_self, y_region.');  


cross_history = [];


n = size(Y,1);
for i = 1:n
    if i == region
        continue;                                  
    end
    increment_curr = Y(i, :);             
    curr_region_history = historicalVec(numDelays_cross, increment_curr.');      
    cross_history = [cross_history, curr_region_history];           % concatenating regions
end


X = [self_history, cross_history];   % getting self and cross history
y = y_region;   

whos X y y_region


%% ------- Dataset sanity ---------
fprintf('\n=== Data sanity ===\n');
fprintf('Region: %d\n', region);
fprintf('T bins: %d\n', length(y));
fprintf('X size: %d x %d (samples x features)\n', size(X,1), size(X,2));
fprintf('y spike count: %d\n', sum(y));
fprintf('y mean per bin: %.4g\n', mean(y));
fprintf('y rate (Hz): %.4f\n', mean(y) / dt);   % requires dt in workspace
fprintf('X nonzero fraction: %.4f\n', nnz(X) / numel(X));


% mdl = fitglm(X,y,'linear','Distribution','poisson','Link','log');


[B, FitInfo] = lassoglm(X, y, 'poisson', 'Standardize', false); % get 100 lambdas

fprintf('\n=== Lambda path ===\n');

% Show what FitInfo actually contains (useful once, then you can comment it out)
disp('FitInfo fields:');
disp(fieldnames(FitInfo));

nLambda = length(FitInfo.Lambda);
fprintf('Num lambdas: %d\n', nLambda);

% Pick indices robustly
idx_mid = round(nLambda/2);

% Try to get min-deviance index, otherwise compute it
idx_min = [];
if isfield(FitInfo, 'IndexMinDeviance')
    idx_min = FitInfo.IndexMinDeviance;
elseif isfield(FitInfo, 'Deviance')
    [~, idx_min] = min(FitInfo.Deviance);
else
    idx_min = idx_mid; % fallback if deviance not available
end

% Try to get 1SE index if available
idx_1se = [];
if isfield(FitInfo, 'Index1SE')
    idx_1se = FitInfo.Index1SE;
end

fprintf('idx_min: %d, lambda(min): %.4g\n', idx_min, FitInfo.Lambda(idx_min));
fprintf('idx_mid: %d, lambda(mid): %.4g\n', idx_mid, FitInfo.Lambda(idx_mid));
if ~isempty(idx_1se)
    fprintf('idx_1se: %d, lambda(1SE): %.4g\n', idx_1se, FitInfo.Lambda(idx_1se));
end

% Build a list of indices you'll report consistently elsewhere
report_idx = unique([idx_min, idx_mid, idx_1se]);
report_idx = report_idx(~isnan(report_idx) & report_idx>=1 & report_idx<=nLambda);



%% gather statistics


fprintf('Number of ')


%% ------------- Deviance ----------------
fprintf('\n=== Deviance diagnostics (manual computation) ===\n');

% Make sure y is a column vector once, globally for this region
y_col = y(:);
y_safe = y_col;
y_safe(y_safe == 0) = eps;

% Null model mean
mu0 = mean(y_col) * ones(size(y_col));
mu0 = max(mu0, eps);

% Null deviance
dev0 = 2 * sum( y_col .* log(y_safe ./ mu0) - (y_col - mu0) );

fprintf('Sanity sizes: y=%dx%d, mu0=%dx%d, X=%dx%d\n', ...
    size(y_col,1), size(y_col,2), size(mu0,1), size(mu0,2), size(X,1), size(X,2));

for ii = report_idx
    bhat = B(:,ii);

    mu = exp(FitInfo.Intercept(ii) + X * bhat);
    mu = mu(:);              % force column
    mu = max(mu, eps);

    dev = 2 * sum( y_col .* log(y_safe ./ mu) - (y_col - mu) );
    fracExplained = 1 - dev/dev0;

    fprintf('idx=%d  lambda=%.3g  dev=%.4g  nullDev=%.4g  devExplained=%.2f%%\n', ...
        ii, FitInfo.Lambda(ii), dev, dev0, 100*fracExplained);
end


%% ------------- Predicted rate vs Observed Rate --------------
fprintf('\n=== Predictive mean check ===\n');
for ii = report_idx
    bhat = B(:,ii);
    mu = exp(FitInfo.Intercept(ii) + X * bhat);  % expected spikes per bin
    
    fprintf('idx=%d  mean(y)=%.4g  mean(mu)=%.4g  mean(mu)/dt=%.4g Hz\n', ...
        ii, mean(y), mean(mu), mean(mu)/dt);
end


ii = idx_min;  % pick one
bhat = B(:,ii);

figure;
stem(1:numDelays_self, bhat(1:numDelays_self), 'filled');
xlabel('Self-history lag');
ylabel('Coefficient');
title(sprintf('Region %d self-history coeffs (lambda=%.3g)', region, FitInfo.Lambda(ii)));
grid on;

figure;
stem(bhat(numDelays_self+1:end), 'filled');
xlabel('Cross-history feature index');
ylabel('Coefficient');
title(sprintf('Region %d cross-history coeffs (lambda=%.3g)', region, FitInfo.Lambda(ii)));
grid on;



