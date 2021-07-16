function [ M_max ] = max_moment_beam( distributed_load, L_beam_length)
%MAX_MOMENT_BEAM Summary of this function goes here
%   Detailed explanation goes here

M_max = (distributed_load*L_beam_length^2)/3;

end

