function [ V_out, quad_spike ] = neuronTbird( V_in, binary_spike_A, binary_spike_B, binary_spike_C )
% Full integrate and fire neuron, just with a very strong synaptic connection to A
% potentially with a time delay.

dt = .0001;
tau = .01;        % membrane time constant

Vreset = -0.065;   % reset voltage
Vth = -0.05;       % threshold voltage
Rm = 1E7;         % resistance

Ie = 1E-6;       % constant external current

quad_spike = 0;	% unless a spike is triggered, the neuron has not yet spiked

dV = ( Vreset - V_in + Rm.*Ie.*(binary_spike_A == 1) + Rm.*Ie.*(binary_spike_B == 1) + Rm.*Ie.*(binary_spike_C == 1) )./tau.*dt;
% NOTE: the above line assumes A, B, and C spike at different times
V_out =  V_in + dV;

if V_out > Vth && binary_spike_A == 1
    quad_spike = 1;
    V_out = Vreset;
elseif V_out > Vth && binary_spike_B == 1
    quad_spike = 2;
    V_out = Vreset;
elseif V_out > Vth && binary_spike_C == 1
    quad_spike = 3;
    V_out = Vreset;
end

end