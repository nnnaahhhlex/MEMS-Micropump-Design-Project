function [ x_stress, y_stress ] = diaphragm_stress( q0_pressure, D, a_plate_width, b_plate_height,diaphragm_thickness, poisson_ratio, x, y)
%DIQPHRAGM_STRESS Summary of this function goes here
%   Detailed explanation goes here

% Located at the center of the diaphragm
% plate_diff_dx = -(17592186044416*(q0_pressure*pi^2)*sin((pi*(x/2))/a_plate_width)*sin((pi*(y/2))/b_plate_height))/(1713638851887625*D*(a_plate_width^2)*(1/a_plate_width^2 + 1/b_plate_height^2)^2);
% plate_diff_dy = -(17592186044416*(q0_pressure*pi^2)*sin((pi*(x/2))/a_plate_width)*sin((pi*(y/2))/b_plate_height))/(1713638851887625*D*(b_plate_height^2)*(1/a_plate_width^2 + 1/b_plate_height^2)^2);
plate_diff_dx = -(17592186044416*(q0_pressure*pi^2)*sin((pi*(x/2))/a_plate_width)*sin((pi*(y/2))/b_plate_height))/(1713638851887625*D*(a_plate_width^2)*(1/a_plate_width^2 + 1/b_plate_height^2)^2);
plate_diff_dy = -(17592186044416*(q0_pressure*pi^2)*sin((pi*(x/2))/a_plate_width)*sin((pi*(y/2))/b_plate_height))/(1713638851887625*D*(b_plate_height^2)*(1/a_plate_width^2 + 1/b_plate_height^2)^2);
 

x_stress = (((-12*D*(diaphragm_thickness/2))/diaphragm_thickness^3))*((plate_diff_dx)+(poisson_ratio*plate_diff_dy));
y_stress = (((-12*D*(diaphragm_thickness/2))/diaphragm_thickness^3))*((plate_diff_dy)+(poisson_ratio*plate_diff_dx));

end

