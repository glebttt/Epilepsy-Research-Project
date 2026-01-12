%% Synthetic spike train generator (network GLM simulation)

function [Y, A] = generate_synthetic_data(N, dt, Tsec, plot_flag)
% GENERATE_SYNTHETIC_DATA Simulate multivariate spike trains from a known network
%
% This function generates synthetic binary spike trains using a discrete-time
% Poisson GLM with self-history (refractory effects) and coupling (network
% interactions) filters. The underlying network connectivity is known and
% returned as the adjacency matrix A (ground truth), allowing evaluation of
% connectivity inference methods such as Granger causality.
%
% Inputs:
%   N         - number of neurons / brain regions
%   dt        - time bin size (seconds)
%   Tsec      - total simulation duration (seconds)
%   plot_flag - if 1, plot example spike trains; if 0, no plotting
%
% Outputs:
%   Y - N x T binary spike matrix (rows = regions, columns = time bins)
%   A - N x N adjacency matrix (ground truth connectivity)
%       A(j,i) = 1 means region j causally influences region i
    %% -------------------------------
    % Simulation parameters
    % -------------------------------
    T  = round(Tsec / dt);  % total number of time bins
    
    % this is a discrete-time approximation of a continuous point process
    % bin size is small enough that at most one spike per bin is reasonable 
    
    Kh = 25;                % self-history length (100 ms) (25 bins)
    Kc = 25;                % coupling history length (100 ms) (25 bins)
    
    % when deciding whether a neuron spikes now, the model looks back 100 ms
    % into the past
    % it does not care about spikes older than 100ms 
    
    % self-history: if a neuron spiked recently
    % that affects how likely it is to fire again shortly after 
    
    % coupling history: if another neuron fired recently, it might influence me
    % now
    % this models communication between neurons or brain regions 
    
    % these self and coupling histories capture:
    
    % 1. refractory periods: 
    % After a neuron fires, for a few milliseconds it becomes temporarily less
    % excitable i.e. less likely to fire 
    % "I just spoke, give me a moment before I speak again"
    
    % 2. short-latency synaptic effects: 
    % synapse: the connection where one neuron influences another 
    % latency: means delay 
    % short-latency synaptic effects means that when one neuron fires, it
    % briefly and quickly increases the chance that another neuron will fire,
    % with a small delay of a few milliseconds 

    %% -------------------------------
    % Baseline firing rates (Hz)
    % -------------------------------
    
    baseline_rate = 10 * ones(N,1);      % moderate firing
    
    % N = number of regions 
    % ones(N,1) --> column vector of ones, size Nx1
    % multiply by 10 --> every neuron has 10 Hz baseline firing rate
    % so each neuron wants to fire around 10 times per second if nothing else
    % is happening 
    
    beta0 = log(baseline_rate * dt);     % baseline log-rate
    
    % dt are the 4 ms time bins 
    % in each bin, a neuron can fire 0 or 1 spikes 
    % baseline_rate * dt because we want to know how many spikes do I expect in
    % each 4 ms bin 
    % on average, in a single 4 ms bin, the neuron will fire 0.04 spikes 
    % since 0.04 < 1, most bins are empty (no spike)
    % occasionally, a spike occurs, roughly 10 spikes per second, matching the
    % target 10 Hz
    % why take the log? 
    % the neuron's firing follows a Poisson process in discrete time
    % this is an assumption of the Poisson GLM
    % baseline term: when there is no influence from self or other neurons, the
    % neuron should fire at its baseline rate 
    % a little bit of maths to show that beta0 is the log :)

    %% -------------------------------
    % Network adjacency matrix (GROUND TRUTH)
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
    
    % each 1 indicates that neuron j influences neuron i via the coupling
    % filter C(j, i, :)
    % each 0 means no direct influence 
    % the matrix encodes which neurons are connected and in which direction 
    % I am generating the data, so I already know the true network connections

    %% -------------------------------
    % Self-history filters (refractory)
    % -------------------------------
    Ah = 2.0;        % refractory strength
    tau_h = 5;       % decay (bins)
    
    h = zeros(N,Kh);
    for i = 1:N
        h(i,:) = -Ah * exp(-(1:Kh)/tau_h);
    end
    
    % 1. FILTER INTUITION
    % a filter is a set of numbers that says how past spikes influence current
    % firing 
    % recent spikes might suppress or increase future spikes 
    
    % modelling refractoriness: after a neuron spikes, it is less likely to
    % fire for a short while 
    
    % so the filter will be negative, strongest immediately after a spike, and
    % then decays towards zero
    
    
    % 2. CODE EXPLANATION
    % Ah is the strength of the refractory effect
    % bigger number --> neuron is more strongly suppressed right after a spike
    
    % tau_h: this is how fast the effect decays over time (in number of bins)
    % larger tau --> slower decay (suppression lasts longer)
    % smaller tau --> shorter refractory period 
    
    % h is a matrix: N regions x Kh time bins 
    % each row will store the self-history filter for that region 
    % Kh = 25 (each region remembers 25 past bins ~ 100 ms)
    
    % the loop:
    % for each row (1:N):
    % 1:Kh --> row vector [1, 2, 3, ..., Kh]
    % exp(...) computes the decaying exponential for each lag 
    % multiplying by -Ah gives the self-history filter values for each lag 
    % result: h(i,:) is a row vector of length Kh representing the refractory
    % effect at each past time step 
    % i.e. each element in h(i, :) tells you how strongly a spike k steps ago
    % suppresses the current spike 
    % right after a spike: strong suppression
    % later: weak suppression 
    % after ~100 ms, almost no effect
    
    % 3. HOW IT IS USED IN SIMULATION
    % for computing eta, multiplying each past spike by the filter value h(i,k) 
    % (for region i  and time step k)
    % adds a negative contribution --> reduces current firing probability
    % this implements the refractory period!!! :)

    %% -------------------------------
    % Coupling filters (propagation)
    % -------------------------------
    B = 1.5;         % coupling strength
    d = 5;           % delay (bins = 20 ms)
    sigma = 2;       % smoothness
    
    % coupling filter: how neuron j affects neuron i 
    % B: how strong the effect is 
    % d: how long it takes for the effect to happen
    % sigma: how spread out the effect is over time 
    
    C = zeros(N,N,Kc);
    % C is a 3D matrix 
    % C(j, i, k) = effect of neuron j on neuron i at lag k
    % C will store how each presynaptic neuron affects each postsynaptic neuron
    % over the past 100 ms
    % presynaptic: the neuron that sends a signal (spike)
    % postsynaptic: the neuron that receives the signal 
    
    kvec = 1:Kc;
    % vector [1, 2, ... , Kc] (time lags), where Kc is the coupling history
    % length 
    
    for j = 1:N
        for i = 1:N
            if A(j,i) == 1
                C(j,i,:) = B * exp(-(kvec - d).^2 / (2*sigma^2));
            end
        end
    end
    
    % loop explanation:
    % loop over all neuron pairs (j, i)
    % only create a filter if A(j:i) == 1
    % A is the adjacency matrix (network connectivity)
    % So we only model connections that exist 
    % exp(...) is the Gaussian curve over lags 
    % therefore, spike from neuron j increases probability that neuron i spikes
    % after a short delay
    % the effect is strongest at lag d, decays around it like a bell curve
    % makes sense to use Gaussian since influence rises as you approach lag = d
    % peaks at lag = d
    % slowly decays after lag = d 
    
    % HOW IT IS USED IN SIMULATION
    % for each presynaptic region j, each past spike Y(j, t-k):
    % multiply by the filter weight C(j, i, k)
    % adds to the postsynaptic region's linear predictor eta
    % this models spike propagation over time

    Y = zeros(N,T);
    % preallocaes spike matrix 
    % rows = regions 
    % columns = time bins 

    %% -------------------------------
    % Recursive simulation
    % -------------------------------    
    t_start = max(Kh, Kc) + 1;
    % model requires up to Kh past bins for self-history
    % model required up to Kc past bins for coupling 
    % e.g. at t=1 there is no past, t=10 u have past up til t-1 = 9
    
    % time
    for t = t_start:T
        % only start simulating spikes once all required history terms exist
        % region
        for i = 1:N
            % simulating spike at time t for each region 
            
            % Linear predictor
            eta = beta0(i);
            
            % Self-history contribution
            % history 
            for k = 1:Kh
                % filter times the spike at lag k
                % Y(i, t-k) means did region i spike exactly k time bins ago
                eta = eta + h(i,k) * Y(i, t-k);
            end
            
            % Coupling contribution
            % region
            for j = 1:N
                % if j not equal to i 
                if j ~= i && A(j,i) == 1
                    % history
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
    
    % NOTE:
    % time bins: t = 1, ..., 25 are never simulated
    % t_start = 26 
    % such a small part of the duration that it's negligible 

    %% -------------------------------
    % Basic sanity checks
    % -------------------------------
    firing_rates = sum(Y,2) / Tsec;
    
    fprintf('\nMean firing rates (Hz):\n');
    for i = 1:N
        fprintf('Region %d: %.2f Hz\n', i, firing_rates(i));
    end
    
    % expected rate is 10 Hz so they should be around that 
    % depending on which regions affect which
    % e.g. in region 1 no incoming excitation --> only self-refractory effect
    % rate slightly below baseline 
    
    % e.g. in region 7 excitation from region 5 and region 6
    % rate above baseline since 5 and 6 already have elevated rates 
    % cascade amplification effect: which means spikes propogate through a
    % directed network, and each stage increases the firing rate of the next,
    % even if individual connections are not very strong 
    % meaning mean firing rates increase downstream due to accumulated excitation
    
    % purpose of basic sanity checks:
    % ensures baseline rates are correct
    % model is stable (no runaway excitation)

    %% -------------------------------------------------------------
    % Plotting the binary spike trains (only if the flag, plot is 1)
    % --------------------------------------------------------------
    if plot_flag 
        fs = 1/dt;
        time = (0:T-1) / fs;
        
        % using zoomed windows for visualisation
        tmin = 0;   % seconds
        tmax = 600;
        
        idx = time >= tmin & time <= tmax;
        
        for i = 1:5
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
    end
end