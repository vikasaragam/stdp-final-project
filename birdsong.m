function [ scenario, W_mat, W_T ] = birdsong( sigmoidX_T, C_neg )

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
V_C = Vreset .* ones(1, length(time_global));
V_T = Vreset .* ones(1, length(time_global));
V_1 = Vreset .* ones(1, length(time_global));
V_2 = Vreset .* ones(1, length(time_global));
V_3 = Vreset .* ones(1, length(time_global));

% Initialize spike time storages
binary_spike_A = zeros( 1, length( time_global ) );
binary_spike_B = zeros( 1, length( time_global ) );
binary_spike_C = zeros( 1, length( time_global ) );
binary_spike_1 = zeros( 1, length( time_global ) );
binary_spike_2 = zeros( 1, length( time_global ) );
binary_spike_3 = zeros( 1, length( time_global ) );

% T will output a 1, 2, or 3 if it was pushed over Vth by A, B, or C
% respectively; if no spike, T outputs a 0
quad_spike_T = zeros( 1, length( time_global ) );

A_spikeTime1 = 0;
A_spikeTime2 = 0;
A_spikeTime3 = 0;
B_spikeTime1 = 0;
B_spikeTime2 = 0;
B_spikeTime3 = 0;
C_spikeTime1 = 0;
C_spikeTime2 = 0;
C_spikeTime3 = 0;
S1_spikeTime = 0;
S2_spikeTime = 0;
S3_spikeTime = 0;

% Initialize current for S
Ie_1 = zeros( 1, length( time_global ) );
Ie_1(1) = 1E-9;
Ie_2 = zeros( 1, length( time_global ) );
Ie_2(1) = 1E-9;
Ie_3 = zeros( 1, length( time_global ) );
Ie_3(1) = 1E-9;

% Initialize weighting functions
C_plus = 0.5;
C_minus = -C_neg.*C_plus; % originally 1.15
lambda = 0.01;			% temporary value, change according to literature

W_max = 10;			% temporary value, change according to calculation

W_mat = zeros( 3, 3, length( time_global ) );
W_T = zeros( 1, length( time_global ) );

sigX_mat = -1.5 .* ones(3, 3);		% CHOOSE WISELY
sigX_T = sigmoidX_T;		% CHOOSE WISELY

alpha = 4;
W_mat(:, :, 1) = W_max ./ (1 + exp(-alpha.*sigX_mat));
W_T(1) = W_max ./ (1 + exp(-alpha.*sigX_T));

spike_interval = 150;
mean_time_A = 50;
mean_time_B = 100;
mean_time_C = 150;
Spike_countdown_A = zeros(1, length( time_global ) );
Spike_countdown_B = zeros(1, length( time_global ) );
Spike_countdown_C = zeros(1, length( time_global ) );
Spike_countdown_A(2) = mean_time_A;
Spike_countdown_B(2) = mean_time_B;
Spike_countdown_C(2) = mean_time_C;

% Creating the spiking patterns of A, B, and T:

for currentTime = 2:(t_final/dt + 1)
    
    if Spike_countdown_A(currentTime) == 0;
        binary_spike_A(currentTime) = 1;
        V_A(currentTime) = Vth;
        Spike_countdown_A(currentTime) = spike_interval;
    end
    
    if Spike_countdown_B(currentTime) == 0;
        binary_spike_B(currentTime) = 1;
        V_B(currentTime) = Vth;
        Spike_countdown_B(currentTime) = spike_interval;
    end

    if Spike_countdown_C(currentTime) == 0;
        binary_spike_C(currentTime) = 1;
        V_C(currentTime) = Vth;
        Spike_countdown_C(currentTime) = spike_interval;
    end
    
    Spike_countdown_A(currentTime+1) = Spike_countdown_A(currentTime) - 1;
    Spike_countdown_B(currentTime+1) = Spike_countdown_B(currentTime) - 1;
    Spike_countdown_C(currentTime+1) = Spike_countdown_C(currentTime) - 1;
    
    [ V_T(currentTime), quad_spike_T(currentTime) ] = ...
        neuronTbird(V_T(currentTime-1), binary_spike_A(currentTime-1), binary_spike_B(currentTime-1), binary_spike_C(currentTime-1));
    
    [ V_1(currentTime), Ie_1(currentTime), binary_spike_1(currentTime) ] = ...
        neuronSbird( 1, V_1(currentTime-1), Ie_1(currentTime-1), binary_spike_A(currentTime-1), binary_spike_B(currentTime-1), binary_spike_C(currentTime-1), quad_spike_T(currentTime-1), W_mat( :, 1, currentTime-1), W_T(currentTime-1) );

    [ V_2(currentTime), Ie_2(currentTime), binary_spike_2(currentTime) ] = ...
        neuronSbird( 2, V_2(currentTime-1), Ie_2(currentTime-1), binary_spike_A(currentTime-1), binary_spike_B(currentTime-1), binary_spike_C(currentTime-1), quad_spike_T(currentTime-1), W_mat( :, 2, currentTime-1), W_T(currentTime-1) );

    [ V_3(currentTime), Ie_3(currentTime), binary_spike_3(currentTime) ] = ...
        neuronSbird( 3, V_3(currentTime-1), Ie_3(currentTime-1), binary_spike_A(currentTime-1), binary_spike_B(currentTime-1), binary_spike_C(currentTime-1), quad_spike_T(currentTime-1), W_mat( :, 3, currentTime-1), W_T(currentTime-1) );
    
    % Now we need to check for spikes and store any spike times in memory
    if binary_spike_A(currentTime) == 1;
        A_spikeTime1 = currentTime;
        A_spikeTime2 = currentTime;
        A_spikeTime3 = currentTime;
    end
    if binary_spike_B(currentTime) == 1;
        B_spikeTime1 = currentTime;
        B_spikeTime2 = currentTime;
        B_spikeTime3 = currentTime;
    end
    if binary_spike_C(currentTime) == 1;
        C_spikeTime1 = currentTime;
        C_spikeTime2 = currentTime;
        C_spikeTime3 = currentTime;
    end
    if binary_spike_1(currentTime) == 1;
        S1_spikeTime = currentTime;
    end
    if binary_spike_2(currentTime) == 1;
        S2_spikeTime = currentTime;
    end
    if binary_spike_3(currentTime) == 1;
        S3_spikeTime = currentTime;
    end
    
    W_mat(:, :, currentTime) = W_mat(:, :, currentTime-1);
    W_T(currentTime) = W_T(currentTime-1);
    
    % Comparisons between A and students
    
    if binary_spike_1(currentTime) == 1 && A_spikeTime1 ~= 0;
        deltaT_A = A_spikeTime1 - currentTime;
        STDP_A = C_plus.*exp( lambda .* deltaT_A );
        sigX_mat(1, 1) = sigX_mat(1, 1) + STDP_A;
        W_mat(1, 1, currentTime) = W_max ./ (1 + exp(-alpha.*sigX_mat(1, 1)));
        A_spikeTime1 = 0;
    elseif binary_spike_A(currentTime) == 1 && S1_spikeTime ~= 0;
        deltaT_A = -(S1_spikeTime - currentTime);
        STDP_A = C_minus.*exp( -lambda .* deltaT_A);
        sigX_mat(1, 1) = sigX_mat(1, 1) + STDP_A;
        W_mat(1, 1, currentTime) = W_max ./ (1 + exp(-alpha.*sigX_mat(1, 1)));
    end
    if binary_spike_2(currentTime) == 1 && A_spikeTime2 ~= 0;
        deltaT_A = A_spikeTime2 - currentTime;
        STDP_A = C_plus.*exp( lambda .* deltaT_A );
        sigX_mat(1, 2) = sigX_mat(1, 2) + STDP_A;
        W_mat(1, 2, currentTime) = W_max ./ (1 + exp(-alpha.*sigX_mat(1, 2)));
        A_spikeTime2 = 0;
    elseif binary_spike_A(currentTime) == 1 && S2_spikeTime ~= 0;
        deltaT_A = -(S2_spikeTime - currentTime);
        STDP_A = C_minus.*exp( -lambda .* deltaT_A);
        sigX_mat(1, 2) = sigX_mat(1, 2) + STDP_A;
        W_mat(1, 2, currentTime) = W_max ./ (1 + exp(-alpha.*sigX_mat(1, 2)));
    end
    if binary_spike_3(currentTime) == 1 && A_spikeTime3 ~= 0;
        deltaT_A = A_spikeTime3 - currentTime;
        STDP_A = C_plus.*exp( lambda .* deltaT_A );
        sigX_mat(1, 3) = sigX_mat(1, 3) + STDP_A;
        W_mat(1, 3, currentTime) = W_max ./ (1 + exp(-alpha.*sigX_mat(1, 3)));
        A_spikeTime3 = 0;
    elseif binary_spike_A(currentTime) == 1 && S3_spikeTime ~= 0;
        deltaT_A = -(S3_spikeTime - currentTime);
        STDP_A = C_minus.*exp( -lambda .* deltaT_A);
        sigX_mat(1, 3) = sigX_mat(1, 3) + STDP_A;
        W_mat(1, 3, currentTime) = W_max ./ (1 + exp(-alpha.*sigX_mat(1, 3)));
    end
    
    % Comparisons between B and students
    if binary_spike_1(currentTime) == 1 && B_spikeTime1 ~= 0;
        deltaT_B = B_spikeTime1 - currentTime;
        STDP_B = C_plus.*exp( lambda .* deltaT_B );
        sigX_mat(2, 1) = sigX_mat(2, 1) + STDP_B;
        W_mat(2, 1, currentTime) = W_max ./ (1 + exp(-alpha.*sigX_mat(2, 1)));
        B_spikeTime1 = 0;
    elseif binary_spike_B(currentTime) == 1 && S1_spikeTime ~= 0;
        deltaT_B = -(S1_spikeTime - currentTime);
        STDP_B = C_minus.*exp( -lambda .* deltaT_B);
        sigX_mat(2, 1) = sigX_mat(2, 1) + STDP_B;
        W_mat(2, 1, currentTime) = W_max ./ (1 + exp(-alpha.*sigX_mat(2, 1)));
    end
    if binary_spike_2(currentTime) == 1 && B_spikeTime2 ~= 0;
        deltaT_B = B_spikeTime2 - currentTime;
        STDP_B = C_plus.*exp( lambda .* deltaT_B );
        sigX_mat(2, 2) = sigX_mat(2, 2) + STDP_B;
        W_mat(2, 2, currentTime) = W_max ./ (1 + exp(-alpha.*sigX_mat(2, 2)));
        B_spikeTime2 = 0;
    elseif binary_spike_B(currentTime) == 1 && S2_spikeTime ~= 0;
        deltaT_B = -(S2_spikeTime - currentTime);
        STDP_B = C_minus.*exp( -lambda .* deltaT_B);
        sigX_mat(2, 2) = sigX_mat(2, 2) + STDP_B;
        W_mat(2, 2, currentTime) = W_max ./ (1 + exp(-alpha.*sigX_mat(2, 2)));
    end
    if binary_spike_3(currentTime) == 1 && B_spikeTime3 ~= 0;
        deltaT_B = B_spikeTime3 - currentTime;
        STDP_B = C_plus.*exp( lambda .* deltaT_B );
        sigX_mat(2, 3) = sigX_mat(2, 3) + STDP_B;
        W_mat(2, 3, currentTime) = W_max ./ (1 + exp(-alpha.*sigX_mat(2, 3)));
        B_spikeTime3 = 0;
    elseif binary_spike_B(currentTime) == 1 && S3_spikeTime ~= 0;
        deltaT_B = -(S3_spikeTime - currentTime);
        STDP_B = C_minus.*exp( -lambda .* deltaT_B);
        sigX_mat(2, 3) = sigX_mat(2, 3) + STDP_B;
        W_mat(2, 3, currentTime) = W_max ./ (1 + exp(-alpha.*sigX_mat(2, 3)));
    end
    
    % Comparisons between C and students
    if binary_spike_1(currentTime) == 1 && C_spikeTime1 ~= 0;
        deltaT_C = C_spikeTime1 - currentTime;
        STDP_C = C_plus.*exp( lambda .* deltaT_C );
        sigX_mat(3, 1) = sigX_mat(3, 1) + STDP_C;
        W_mat(3, 1, currentTime) = W_max ./ (1 + exp(-alpha.*sigX_mat(3, 1)));
        C_spikeTime1 = 0;
    elseif binary_spike_C(currentTime) == 1 && S1_spikeTime ~= 0;
        deltaT_C = -(S1_spikeTime - currentTime);
        STDP_C = C_minus.*exp( -lambda .* deltaT_C);
        sigX_mat(3, 1) = sigX_mat(3, 1) + STDP_C;
        W_mat(3, 1, currentTime) = W_max ./ (1 + exp(-alpha.*sigX_mat(3, 1)));
    end
    if binary_spike_2(currentTime) == 1 && C_spikeTime2 ~= 0;
        deltaT_C = C_spikeTime2 - currentTime;
        STDP_C = C_plus.*exp( lambda .* deltaT_C );
        sigX_mat(3, 2) = sigX_mat(3, 2) + STDP_C;
        W_mat(3, 2, currentTime) = W_max ./ (1 + exp(-alpha.*sigX_mat(3, 2)));
        C_spikeTime2 = 0;
    elseif binary_spike_C(currentTime) == 1 && S2_spikeTime ~= 0;
        deltaT_C = -(S2_spikeTime - currentTime);
        STDP_C = C_minus.*exp( -lambda .* deltaT_C);
        sigX_mat(3, 2) = sigX_mat(3, 2) + STDP_C;
        W_mat(3, 2, currentTime) = W_max ./ (1 + exp(-alpha.*sigX_mat(3, 2)));
    end
    if binary_spike_3(currentTime) == 1 && C_spikeTime3 ~= 0;
        deltaT_C = C_spikeTime3 - currentTime;
        STDP_C = C_plus.*exp( lambda .* deltaT_C );
        sigX_mat(3, 3) = sigX_mat(3, 3) + STDP_C;
        W_mat(3, 3, currentTime) = W_max ./ (1 + exp(-alpha.*sigX_mat(3, 3)));
        C_spikeTime3 = 0;
    elseif binary_spike_C(currentTime) == 1 && S3_spikeTime ~= 0;
        deltaT_C = -(S3_spikeTime - currentTime);
        STDP_C = C_minus.*exp( -lambda .* deltaT_C);
        sigX_mat(3, 3) = sigX_mat(3, 3) + STDP_C;
        W_mat(3, 3, currentTime) = W_max ./ (1 + exp(-alpha.*sigX_mat(3, 3)));
    end
    
end

scenario = 0;

% 1 = good learning
% 0 = incomplete learning

if W_mat(1, 1, end) >= 9 && W_mat(1, 2, end) >= 9 && W_mat(2, 2, end) >= 9 ...
        && W_mat(2, 3, end) >= 9 && W_mat(3, 3, end) >= 9 &&...
        W_mat(3, 1, end) >= 9 && W_mat(1, 3, end) <= 1 && W_mat(2, 1, end) <= 1 && W_mat(3, 2, end) <= 1
    scenario = 1;
end

% for i = 1: 3
%     for j = 1: 3
%         plot( time_global, squeeze( W_mat(i,j,:) ) )
%         legend('W_mat')
%         xlabel('Time (s)'), ylabel('Synaptic Weights (a.u.)')
%         title('Synaptic Weights between Four Neurons over Time')
%     end
% end


% figure(1)
% plot( time_global, V_1, time_global, V_2, time_global, V_3 )
% 
% figure(2)
% plot( time_global, V_A, time_global, V_B, time_global, V_C )
% legend('W_mat', 'W_T')
% xlabel('Time (s)'), ylabel('Synaptic Weights (a.u.)')
% title('Synaptic Weights between Four Neurons over Time')


end