function [ max_stress_beam ] = max_beam_stress( beam_height, max_moment, I_moment_of_inertia )
%MAX_BEAM_STRESS Summary of this function goes here
%   Detailed explanation goes here

max_stress_beam = (beam_height * max_moment) / I_moment_of_inertia;

end

