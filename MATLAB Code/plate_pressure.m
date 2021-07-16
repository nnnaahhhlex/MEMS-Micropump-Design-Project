function [ q0_pressure ] = plate_pressure( a_plate_width, b_plate_height, D, stroke_volume )
%PLATE_PRESSURE Summary of this function goes here
%   Detailed explanation goes here

q0_pressure = (stroke_volume*(pi^4)*D*((1/a_plate_width^2)+(1/b_plate_height^2))^2)/((a_plate_width/pi)-(a_plate_width/pi)*cos(pi*(b_plate_height/a_plate_width))*((b_plate_height/pi)-(b_plate_height/pi)*cos(pi*(a_plate_width/b_plate_height))));

end

