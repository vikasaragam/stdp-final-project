function [ V_out, Ie_out, binary_spike ] = neuronSbird( number, V_in, Ie_in, binary_spike_A, binary_spike_B, binary_spike_C, quad_spike_T, W_in_ABC, W_in_T )
% Unlike A and B, the students external current will vary based on the timing of the students
% spikes with respect to the spikes from A and B. The external current follows a differential
% equation that incorporates a kick corresponding a spike from A or B, and an exponential
% decay (i.e. the decreasing influence of a spike from A or B over time).

dt = .0001;
tau_V = .01;        % membrane voltage time constant
tau_I = .001;        % membrane current time constant

Vreset = -0.065;   % reset voltage
Vth = -0.05;       % threshold voltage
Rm = 1E7;         % resistance

W_in_A = W_in_ABC(1);
W_in_B = W_in_ABC(2);
W_in_C = W_in_ABC(3);

I_const = 3E-9;

if number == 1 && ( quad_spike_T == 1 || quad_spike_T == 3 )
    dIe = ( -Ie_in + (rand >= 0.25).*I_const + ...
        (5E-9).*( W_in_A.*(binary_spike_A == 1) +  W_in_B.*(binary_spike_B == 1) + W_in_C.*(binary_spike_C == 1) ) + (5e-9).*W_in_T )./tau_I.*dt; % NOTE INFLUENCE OF TEACHER
    Ie_out = dIe + Ie_in;
elseif number == 2 && ( quad_spike_T == 1 || quad_spike_T == 2 )
    dIe = ( -Ie_in + (rand >= 0.25).*I_const + ...
        (5E-9).*( W_in_A.*(binary_spike_A == 1) +  W_in_B.*(binary_spike_B == 1) + W_in_C.*(binary_spike_C == 1) ) + (5e-9).*W_in_T )./tau_I.*dt; % NOTE INFLUENCE OF TEACHER
    Ie_out = dIe + Ie_in;
elseif number == 3 && ( quad_spike_T == 2 || quad_spike_T == 3 )
    dIe = ( -Ie_in + (rand >= 0.25).*I_const + ...
        (5E-9).*( W_in_A.*(binary_spike_A == 1) +  W_in_B.*(binary_spike_B == 1) + W_in_C.*(binary_spike_C == 1) ) + (5e-9).*W_in_T )./tau_I.*dt; % NOTE INFLUENCE OF TEACHER
    Ie_out = dIe + Ie_in;
else
    dIe = ( -Ie_in + (rand >= 0.25).*I_const + ...
        (5E-9).*( W_in_A.*(binary_spike_A == 1) +  W_in_B.*(binary_spike_B == 1) + W_in_C.*(binary_spike_C == 1) ) )./tau_I.*dt; % NO TEACHER INFLUENCE
    Ie_out = dIe + Ie_in;
end

binary_spike = 0;	% unless a spike is triggered, the neuron has not yet spiked

dV = ( Vreset - V_in + Rm.*Ie_out )./tau_V.*dt;
V_out =  V_in + dV;

if V_out > Vth
    binary_spike = 1;
    V_out = Vreset;
end