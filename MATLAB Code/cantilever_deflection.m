function [ y_def ] = cantilever_deflection( E,I_moment_of_inertia,w, L,x_length )
%CANTILEVER_DEFLECTION Summary of this function goes here
%   Detailed explanation goes here

y_def = (w/(E*I_moment_of_inertia))*((-((x_length^4)/24))+((x_length*L^3)/6)-((L^4)/8));

end
