sigmoidX_T = ( -3: 0.5: 3 );
C_minus = ( 0.75: 0.05: 1.5 );
numTrials = 100;
learnRate = zeros( length(C_minus), length(sigmoidX_T) );

for k = 1: length(C_minus)
    for i = 1:length(sigmoidX_T);
        for j = 1:numTrials;
            learnRate(k,i) = learnRate(k,i) + birdsong( sigmoidX_T(i), C_minus(k) );
        end
    end
end

learnRate = learnRate./numTrials;

%% Plots of learning rates for fixed initial weight of teacher

plot( C_minus, learnRate(:,7:13) )
xlabel('Magnitude of STDP Inhibition Exponential (/0.5)');
ylabel('Probability of Proper Learning');
title('Proper Learning Probabilities vs. Inhibition Magnitude, Various Initial Teacher Weights');
legend( 'W = 5', 'W = 8.8080', 'W = 9.8201', 'W = 9.9753', 'W = 9.9966', 'W = 9.9995', 'W = 9.9999' );

%% PLOT OF WEIGHTS FOR ONE TEST TRIAL OF FIXED SIGMOIDX_T, C_MINUS
tau = 0.01;
dt = .0001;
t_final = 500.*tau;  % Creates a range for the time
time_global = 0: dt: t_final;

figure(3)
hold on
for i = 1: 3
    for j = 1: 3
        plot( time_global, squeeze( W_mat(i,j,:) ) )
        legend('W_mat')
        xlabel('Time (s)'), ylabel('Synaptic Weights (a.u.)')
        title('Synaptic Weights between Four Neurons over Time')
    end
end
hold off