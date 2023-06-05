%% STDP
dt = 0.0001;
bound = 0.02;
C_plus = 0.5;
C_minus = -1.05.*C_plus;
L = 100;

time_neg = -bound:dt:0;
time_pos = 0:dt:bound;

temp_left = C_plus.*exp( L.* time_neg );
temp_right = C_minus.*exp( -L.* time_pos);

figure(5)
plot(time_neg, temp_left, time_pos, temp_right);
xlabel('Change in Spike Times of S and Input Neuron (s)');
ylabel('Change in Nonlinearity Variable (/1)');
title('STDP Intermediate Function');

%% Nonlinearity
W_max = 10;
sigX = ( -3: 0.1: 3 );
Weight = W_max ./ (1 + exp(-alpha.*sigX));

figure(6)
plot( sigX, Weight );
xlabel('Change in Nonlinearity Variable (/1)');
ylabel('Weight (/1)');
title('Sigmoidal Nonlinearity between STDP and Weight');