sigmoidX_T = -3: 0.5: 3;
numTrials = 100;
learnings = zeros(numTrials, length(sigmoidX_T));

for i = 1:length(sigmoidX_T);
	for j = 1:numTrials;
		learnings(j, i) = testWeightT( sigmoidX_T(i) );
	end
end

learnRate = zeros(1, length(sigmoidX_T));

for i = 1:length(sigmoidX_T);
	learnRate(i) = length(find(learnings(:, i) == 1)) ./ numTrials;
end
%%
W_max = 10;
alpha = 1;
plot(W_max ./ (1 + exp(-alpha .* sigmoidX_T)), learnRate, 'o')
xlabel('Synaptic Weight of T and S (/1)')
ylabel('Probability of Student Following A and Ignoring B')
title('Probability of Good Learning vs. Synaptic Weight of T and S');
%%
alphas = 0.5: 0.5: 4.5;
plot( alphas, alphaVarLearnRate );
xlabel('\alpha of STDP Sigmoid (/1)');
ylabel('Probability of Student Following A and Ignoring B');
title('Probability of Good Learning vs. \alpha of STDP Sigmoid');
legend('W = 0.0001', 'W = 0.0005', 'W = 0.0034', 'W = 0.0247', ...
    'W = 0.1799', 'W = 1.1920', 'W = 5', 'W = 8.8080', 'W = 9.8201', ...
    'W = 9.9753', 'W = 9.9966', 'W = 9.9995', 'W = 9.9999', ...
    'Location', 'eastoutside' );

%%
Cminus = [ 0.9: 0.05: 1.25, 1.35: 0.1: 1.65 ];
plot( Cminus, CminusVarLearnRate );
xlabel('Magnitude of STDP Inhibition Exponential (/0.5)');
ylabel('Probability of Student Following A and Ignoring B');
title('Probability of Good Learning vs. Inhibition Magnitude, Various Initial Teacher Weights');
legend('W = 0.0001', 'W = 0.0005', 'W = 0.0034', 'W = 0.0247', ...
    'W = 0.1799', 'W = 1.1920', 'W = 5', 'W = 8.8080', 'W = 9.8201', ...
    'W = 9.9753', 'W = 9.9966', 'W = 9.9995', 'W = 9.9999', ...
    'Location', 'eastoutside' );