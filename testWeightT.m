function [ scenario ] = testWeightT( sigmoidX_T )

% Initialization code
% Initialize time variables
tau = 0.01;
dt = .0001;
t_final = 500.*tau;  % Creates a range for the time
time_global = 0: dt: t_final;

% Initialize voltage storages
Vreset = -0.065;
Vth = -0.05;
V_A = Vreset .* ones(1, length(time_global));
V_B = Vreset .* ones(1, length(time_global));
V_T = Vreset .* ones(1, length(time_global));
V_S = Vreset .* ones(1, length(time_global));

% V_A(1) = Vreset;
% V_B(1) = Vreset;
% V_T(1) = Vreset;
% V_S(1) = Vreset;

% Initialize spike time storages
binary_spike_A = zeros( 1, length( time_global ) );
binary_spike_B = zeros( 1, length( time_global ) );
binary_spike_T = zeros( 1, length( time_global ) );
binary_spike_S = zeros( 1, length( time_global ) );

A_spikeTime = 0;
B_spikeTime = 0;
S_spikeTime = 0;

% Initialize current for S
Ie_S = zeros( 1, length( time_global ) );
Ie_S(1) = 1E-9;

% Initialize weighting functions
C_plus = 0.5;
C_minus = -1.05.*C_plus;
lambda = 0.01;			% temporary value, change according to literature

nonlin_exp = 2;
W_max = 10;			% temporary value, change according to calculation

W_A = zeros( 1, length( time_global ) );
W_B = zeros( 1, length( time_global ) );
W_T = zeros( 1, length( time_global ) );		% for atrophy of T as S learns

sigX_A = -1.5;		% CHOOSE WISELY
sigX_B = -1.5;		% CHOOSE WISELY
sigX_T = sigmoidX_T;		% CHOOSE WISELY

alpha = 1;
W_A(1) = W_max ./ (1 + exp(-alpha.*sigX_A));
W_B(1) = W_max ./ (1 + exp(-alpha.*sigX_B));
W_T(1) = W_max ./ (1 + exp(-alpha.*sigX_T));

% W_A(1) = W_max .* (sigX_A.^nonlin_exp) ./ (1 + sigX_A.^nonlin_exp);
% W_B(1) = W_max .* (sigX_B.^nonlin_exp) ./ (1 + sigX_B.^nonlin_exp);
% W_T(1) = W_max .* (sigX_T.^nonlin_exp) ./ (1 + sigX_T.^nonlin_exp);

mean_time_A = 100;
mean_time_B = 100;
Spike_countdown_A = zeros(1, length( time_global ) );
Spike_countdown_B = zeros(1, length( time_global ) );
Spike_countdown_A(1) = poissrnd(mean_time_A);
Spike_countdown_B(1) = poissrnd(mean_time_B);

% Creating the spiking patterns of A, B, and T:

for currentTime = 2:(t_final/dt + 1)
    
    if Spike_countdown_A(currentTime) == 0;
        binary_spike_A(currentTime) = 1;
        V_A(currentTime) = Vth;
        Spike_countdown_A(currentTime) = poissrnd(mean_time_A) + 1;
    end
    
    if Spike_countdown_B(currentTime) == 0;
        binary_spike_B(currentTime) = 1;
        V_B(currentTime) = Vth;
        Spike_countdown_B(currentTime) = poissrnd(mean_time_B) + 1;
    end
    
    Spike_countdown_A(currentTime+1) = Spike_countdown_A(currentTime) - 1;
    Spike_countdown_B(currentTime+1) = Spike_countdown_B(currentTime) - 1;
    
    % At every time step, we let A, B, T, and S integrate their respective input currents
%     [ V_A(currentTime+1), binary_spike_A(currentTime+1) ] = neuronA(V_A(currentTime));
    
%     [ V_B(currentTime+1), binary_spike_B(currentTime+1) ] = neuronB(V_B(currentTime));
    
    [ V_T(currentTime), binary_spike_T(currentTime) ] = neuronT(V_T(currentTime-1), binary_spike_A(currentTime-1));
    
    [ V_S(currentTime), Ie_S(currentTime), binary_spike_S(currentTime) ] = neuronS( V_S(currentTime-1), Ie_S(currentTime-1), binary_spike_A(currentTime-1), binary_spike_B(currentTime-1), binary_spike_T(currentTime-1), W_A(currentTime-1), W_B(currentTime-1), W_T(currentTime-1) );
    
    % Now we need to check for spikes and store any spike times in memory
    if binary_spike_A(currentTime) == 1;
        A_spikeTime = currentTime;
    end
    if binary_spike_B(currentTime) == 1;
        B_spikeTime = currentTime;
    end
    if binary_spike_S(currentTime) == 1;
        S_spikeTime = currentTime;
    end
    
    W_A(currentTime) = W_A(currentTime-1);
    W_B(currentTime) = W_B(currentTime-1);
    W_T(currentTime) = W_T(currentTime-1);
    
    if binary_spike_S(currentTime) == 1 && A_spikeTime ~= 0;
        deltaT_A = A_spikeTime - currentTime;
        STDP_A = C_plus.*exp( lambda .* deltaT_A );
        sigX_A = sigX_A + STDP_A;
        W_A(currentTime) = W_max ./ (1 + exp(-alpha.*sigX_A));
        A_spikeTime = 0;
    elseif binary_spike_A(currentTime) == 1 && S_spikeTime ~= 0;
        deltaT_A = -(S_spikeTime - currentTime);
        STDP_A = C_minus.*exp( -lambda .* deltaT_A);
        sigX_A = sigX_A + STDP_A;
        W_A(currentTime) = W_max ./ (1 + exp(-alpha.*sigX_A));
        % S_spikeTime = 0;
    end
    
    if binary_spike_S(currentTime) == 1 && B_spikeTime ~= 0;
        deltaT_B = B_spikeTime - currentTime;
        STDP_B = C_plus.*exp( lambda .* deltaT_B );
        sigX_B = sigX_B + STDP_B;
        W_B(currentTime) = W_max ./ (1 + exp(-alpha.*sigX_B));
        B_spikeTime = 0;
    elseif binary_spike_B(currentTime) == 1 && S_spikeTime ~= 0;
        deltaT_B = -(S_spikeTime - currentTime);
        STDP_B = C_minus.*exp( -lambda .* deltaT_B);
        sigX_B = sigX_B + STDP_B;
        W_B(currentTime) = W_max ./ (1 + exp(-alpha.*sigX_B));
        % S_spikeTime = 0;
    end
    
    
end

scenario = 0;

% 1 = good learning (follow A, ignore B)
% 2 = hyperactive learning (follow both)
% 3 = opposite learning (ignore A, follow B)
% 4 = inhibitory learning (ignores both)
% 0 = ambiguous learning (somewhere in between)

if W_A(end) >= 9 && W_B(end) <= 1
    scenario = 1;
elseif W_A(end) >= 9 && W_B(end) >= 9
    scenario = 2;
elseif W_A(end) <= 1 && W_B(end) >= 9
    scenario = 3;
elseif W_A(end) <= 1 && W_B(end) <= 1
    scenario = 4;
end

figure(2)
plot(time_global, W_A, time_global, W_B, time_global, W_T)
legend('W_A', 'W_B', 'W_T')
xlabel('Time (s)'), ylabel('Synaptic Weights (a.u.)')
title('Synaptic Weights between Four Neurons over Time')

end