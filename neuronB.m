function [ V_out, binary_spike ] = neuronB( V_in )
% Function outputs the updated voltage, V_0 + dV.
% Also outputs either the value 0 corresponding to no spike or 1 if the neuron spikes.

dt = .0001;
tau = .01;        % membrane time constant

Vreset = -0.065;   % reset voltage
Vth = -0.05;       % threshold voltage
Rm = 1E7;         % resistance

Ie = 2E-9;       % constant external current

binary_spike = 0;	% unless a spike is triggered, the neuron has not yet spiked

dV = ( Vreset - V_in + Rm.*Ie.*(rand>=.25) )./tau.*dt;
V_out =  V_in + dV;

if V_out > Vth
    binary_spike = 1;
    V_out = Vreset;
end