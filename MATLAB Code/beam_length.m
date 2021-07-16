function [ L_beam_length ] = beam_length( E, I_moment_of_inertia, beam_frequency, distributed_load )
%BEAM_LENGTH Summary of this function goes here
%   Detailed explanation goes here

L_beam_length = (((((beam_frequency/0.56)^2)*distributed_load)^-1)*E*I_moment_of_inertia)^(1/4);

end

