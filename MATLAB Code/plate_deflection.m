function [ w_plate_deflection ] = plate_deflection( a_plate_width, b_plate_height, D, q0_pressure, x, y)
%PLATE_DEFORMATION Summary of this function goes here
%   Detailed explanation goes here

w_plate_deflection = (q0_pressure/((pi^4)*D))*(((1/(a_plate_width^2))+(1/(b_plate_height^2)))^-2)*(sin(pi*(x/a_plate_width)))*(sin(pi*(y/b_plate_height)));

end

