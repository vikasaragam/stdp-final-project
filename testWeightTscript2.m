alphas_T = ( 0.5: 0.5: 4 );
sigmoidX_T = ( 0: 0.5: 3 );
numTrials = 100;
learnRate = zeros(length(alphas_T), length(sigmoidX_T));

for k = 1:length(alphas_T);
    
    tic
    
    learnings = zeros(numTrials, length(sigmoidX_T));
    
    for i = 1:length(sigmoidX_T);
        for j = 1:numTrials;
            learnings(j, i) = testWeightT2( sigmoidX_T(i), alphas_T(k) );
        end
    end
    
%     learnRate = zeros(1, length(sigmoidX_T));
    
    for i = 1:length(sigmoidX_T);
        learnRate(k, i) = length(find(learnings(:, i) == 1)) ./ numTrials;
    end
    
    toc
    
end
%%
W_max = 10;
hold on
for k = 1:length(sigmoidX_T);
plot( alphas_T, learnRate(:, k) )
% plot( W_max ./ (1 + exp( -alphas_T .* sigmoidX_T(k)) ), learnRate(:, k) )
end
xlabel('\alpha of Teacher Decay Sigmoid (/1)')
ylabel('Probability of Student Following A and Ignoring B')
title('Probability of Good Learning vs. \alpha of Teacher Decay Sigmoid');
legend( 'W = 5', 'W = 8.8080', 'W = 9.8201', 'W = 9.9753', 'W = 9.9966', ...
    'W = 9.9995', 'W = 9.9999', 'Location', 'eastoutside' );
