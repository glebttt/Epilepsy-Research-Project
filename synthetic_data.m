%% ============================================================
%  Multivariate Point-Process GLM Spike Train Simulation
%  Synthetic Epilepsy Propagation Network
% =============================================================

clear; clc; rng(1);

%% -------------------------------
% Simulation parameters
% -------------------------------
N  = 7;                 % number of brain regions
dt = 0.004;             % time bin (4 ms)
Tsec = 3600;             % total duration (seconds)
T  = round(Tsec / dt);  % number of time bins

Kh = 25;                % self-history length (100 ms)
Kc = 25;                % coupling history length (100 ms)

fprintf('Simulating %d regions for %.1f seconds (%d bins)\n', ...
        N, Tsec, T);

%% -------------------------------
% Baseline firing rates (Hz)
% -------------------------------

% to be determiened based on impirical data later 
rate_hz = zeros(N,1);
rate_hz(1)   = 0.25;   % 15/min
rate_hz(2:4) = 0.12;   % 7.2/min
rate_hz(5:7) = 0.05;   % 3/min

rate_hz = rate_hz .* (1 + 0.30*randn(N,1));
rate_hz = max(rate_hz, 1e-3);

beta0 = log(rate_hz * dt);   
   % baseline log-rate

%% -------------------------------
% Network adjacency matrix
% A(j,i) = 1 means j -> i
% -------------------------------
A = zeros(N,N);

A(1,2) = 1;
A(1,3) = 1;
A(1,4) = 1;
A(2,5) = 1;
A(3,5) = 1;
A(4,6) = 1;
A(5,7) = 1;
A(6,7) = 1;

%% -------------------------------
% Self-history filters (refractory)
% -------------------------------
Ah = 3.0;        % refractory strength
tau_h = 5;       % decay (bins)

h = zeros(N,Kh);
for i = 1:N
    h(i,:) = -Ah * exp(-(1:Kh)/tau_h);
end

%% -------------------------------
% Coupling filters (propagation)
% -------------------------------
d = 5;           % delay (bins = 20 ms)
sigma = 2;       % smoothness

C = zeros(N,N,Kc);
kvec = 1:Kc;

B = zeros(N,N);

% strong/medium/weak links
B(1,2) = 1;
B(1,3) = 1;
B(1,4) = 1;
B(2,5) = 0.5;
B(3,5) = 0.5;

B(4,6) = 0.25;
B(5,7) = 0.25;
B(6,7) = 0.25;






for j = 1:N
    for i = 1:N
        if A(j,i) == 1
            C(j,i,:) = B(j,i) * exp(-(kvec - d).^2 / (2*sigma^2));
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



save('initial_synthetic_network.mat','Y');