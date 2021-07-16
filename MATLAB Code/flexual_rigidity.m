function [ D ] = flexual_rigidity( E, diaphragm_thickness, poisson_ratio )
%FLEXUAL_RIGIDITY Summary of this function goes here
%   Detailed explanation goes here

D = (E*diaphragm_thickness^3)/(12*(1-(poisson_ratio)^2));

end

