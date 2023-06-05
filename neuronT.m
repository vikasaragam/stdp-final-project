function [ V_out, binary_spike ] = neuronT( V_in, binary_spike_a )
% Full integrate and fire neuron, just with a very strong synaptic connection to A
% potentially with a time delay.

dt = .0001;
tau = .01;        % membrane time constant

Vreset = -0.065;   % reset voltage
Vth = -0.05;       % threshold voltage
Rm = 1E7;         % resistance

Ie = 1E-6;       % constant external current

binary_spike = 0;	% unless a spike is triggered, the neuron has not yet spiked

dV = ( Vreset - V_in + Rm.*Ie.*(binary_spike_a == 1))./tau.*dt;
V_out =  V_in + dV;

if V_out > Vth
    binary_spike = 1;
    V_out = Vreset;
end