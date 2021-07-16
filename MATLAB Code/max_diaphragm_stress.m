function [ max_edge_stress ] = max_diaphragm_stress( q0_pressure, diaphragm_thickness, a_plate_width, b_plate_height )
%MAX_DIAPHRAGM_STRESS Summary of this function goes here
%   Detailed explanation goes here

max_edge_stress = (q0_pressure*b_plate_height^2)/((2*diaphragm_thickness^2)*(0.623*((b_plate_height/a_plate_width)^6)+1));


end

