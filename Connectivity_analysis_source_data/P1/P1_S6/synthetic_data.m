%% ============================================================
%  Multivariate Point-Process GLM Spike Train Simulation
%  Synthetic Epilepsy Propagation Network
% =============================================================

clear; clc; rng(1);

%% -------------------------------
% Simulation parameters
% -------------------------------
N  = 2;                 % number of brain regions
dt = 0.004;             % time bin (4 ms)
Tsec = 3600;             % total duration (seconds)
T  = round(Tsec / dt);  % number of time bins

Kh = 5;                % self-history length (100 ms)
Kc = Kh;                % coupling history length (100 ms)

fprintf('Simulating %d regions for %.1f seconds (%d bins)\n', ...
        N, Tsec, T);

%% -------------------------------
% Baseline firing rates (Hz)
% -------------------------------


% come back
% to be determiened based on impirical data later 
rate_hz = zeros(N,1);
rate_hz(1)   = 0.3;    % hub: 15/min
rate_hz(2:4) = 0.15;   % others: 7.5/min 

beta0 = log(rate_hz * dt);   
   % baseline log-rate

%% -------------------------------
% Network adjacency matrix
% A(j,i) = 1 means j -> i
% -------------------------------
A = zeros(N,N);
A(1,2) = 1; % uncomment to revert
% A(1,3) = 1;
% A(1,4) = 1;

%% -------------------------------
% Self-history filters (refractory)
% -------------------------------
Ah = 3.5;        % refractory strength
tau_h = 1.2;       % decay (bins) => 2 bins = 16 ms

h = zeros(N,Kh);
for i = 1:N
    h(i,:) = -Ah * exp(-(1:Kh)/tau_h);
end

%% -------------------------------
% Coupling filters (propagation)
% -------------------------------

d = 2;           % delay (bins) => 3 bins = 12 ms
sigma = 0.8;     % smoothness (bins) 

C = zeros(N,N,Kc);
kvec = 1:Kc;

B = zeros(N,N);
B(1,2) = 1.1;
B(1,3) = 0.9;
B(1,4) = 1.0;

for j = 1:N
    for i = 1:N
        if A(j,i) ~= 0
            C(j,i,:) = B(j,i) * exp(-(kvec - d).^2 / (2*sigma^2));
        
        else
            disp('zero');
        end
    end
end
%% -------------------------------
% Initialise spike trains
% -------------------------------
Y = zeros(N,T);

%% -------------------------------
% Recursive simulation
% -------------------------------
t_start = max(Kh, Kc) + 1;

for t = t_start:T
    for i = 1:N
        
        % Linear predictor
        eta = beta0(i);
        
        % Self-history contribution
        for k = 1:Kh
            eta = eta + h(i,k) * Y(i, t-k);
        end
        
        % Coupling contribution
        for j = 1:N
            if j ~= i && A(j,i) == 1
                for k = 1:Kc
                    eta = eta + C(j,i,k) * Y(j, t-k);
                end
            end
        end
        
        % Conditional intensity
        lambda = exp(eta);
        
        % Poisson spike generation
        spike = poissrnd(lambda);
        
        % Enforce at most one spike per bin
        Y(i,t) = min(spike, 1);
    end
end

fprintf('Simulation complete.\n');

%% -------------------------------
% Basic sanity checks
% -------------------------------
firing_rates = sum(Y,2) / Tsec;

fprintf('\nMean firing rates (Hz):\n');
for i = 1:N
    fprintf('Region %d: %.2f Hz\n', i, firing_rates(i));
end



%% --------------------------------
% plotting the binary spike trains
% ---------------------------------

fs = 1/dt;
time = (0:T-1) / fs;

% using zoomed windows for visualisation
tmin = 0;   % seconds
tmax = Tsec;

idx = time >= tmin & time <= tmax;

for i = 1:N
    figure;
    plot(time(idx), Y(i,idx), '-o', ...
         'Color', 'b', ...
         'MarkerSize', 4, ...
         'LineWidth', 1);
    xlabel('Time (s)');
    ylabel('Spike response');
    title(sprintf('Region %d Spike Train', i));
    ylim([-0.1 1.1]);
    grid on;
end


%% ----------
% to do
% -----------

% to incorporate clinical data: 
% define beta0 using mean firing rates from clinical data
% filter lengths (self and cross-coupling)
% Ah and B values 
% A (network design)



save('test_2_node.mat','Y');