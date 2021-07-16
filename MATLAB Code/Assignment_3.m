clc
clear all
close all

%******* Constants ******

%***** Micropump design specifications *****

Q = 1; % Flowrate: 0.095 mL/min --- 1mL = 0.000001 m3
op_f = 1000; % Operational frequency: 1000 Hz
delta_p = 10000; % Generated pressure: 10000 Pa
res_f = op_f*10; % Diaphragm resonance frequency: 10000 Hz
cant_f = op_f*10*10; % Cantilever resonance frequency: 100000 Hz
safety_factor = 2;

stroke_volume = Q*(0.000001/60)/op_f; %Calculating the stoke volume and converting it to m3
fluid_velocity = 0.1; % Unit: m/s
compressibility = 45.8e-11; %Unit: Pa-1
chamber_volume = stroke_volume/(delta_p*compressibility); % Unit: m3

%****** Silicon mechanical properties *******

density = 2330; % Silicon densiy = 2330 kg/m3
poisson_ratio = 0.28; % Poisson Ratio = 0.28
E = 112*10^9; % Modulus of Elasticity = 112 GPa = 112*10^9 Pa (Values from solidworks)
yield_strength = 120*10^6; %Yield strength = 120 MPa = 120*10^6 Pa (Values from Solidworks)
design_stress = yield_strength/safety_factor

%------------------ Inlet/Outlet Size Calculation ----------------

inlet_area = Q*(0.000001/60)/fluid_velocity; %Calculating the cross-sectional area of inlet
inlet_side_length = sqrt(inlet_area); %Assuming the inlet is a square, we can calculate the side length

%================== Diaphragm Analysis     =========================
%==================       Actuator         ========================= 

%------------ Flexural Rigidity Calculation -----------

diaphragm_thickness = (0.0001:0.00005:0.001); %100 micron to 2000 micron thick plate
n = length(diaphragm_thickness); %Determine the size of the array

for i=1:1:n
%dt(1 (one) row, i col) i = n which is the size of the dd array

D(i) = flexual_rigidity(E,diaphragm_thickness(i),poisson_ratio);

end

%------------ Plate Pressure Calculation -----------

a_plate_width = (0.0001:0.0001:0.002); %100 micron to 2000 micron width plate
b_plate_height = (0.0001:0.0001:0.002) ;%100 micron to 2000 micron height plate
nn = length(a_plate_width);

%Make a value larger or equal to b value, never greater
%Selecting desired dimensions for the diaphragm

 a = 16;
 b = 10;
 d = 1;
 dt = d;


for i=1:1:nn %Rows
    
    for ii=1:1:nn %Col
        % pa(rows,col) - For each row, D and a width remain constant *Change value*->
        q0_pressure(i,ii) = plate_pressure( a_plate_width(i), b_plate_height(ii), D(d), stroke_volume); %Unit: Pa 
        actuation_force(i,ii) = q0_pressure(i,ii)*a_plate_width(i)*b_plate_height(ii); % Units: N
    end


end


%------------ Plate Deflection Calculation -----------
%Max deflection

for i=1:1:n %Rows
    
    for ii=1:1:nn %Col
        % pa(rows,col)
    w_max_plate_deflection(i,ii) = plate_deflection(a_plate_width(i), b_plate_height(ii), D(d), q0_pressure(i,ii), (a_plate_width(i)/2), (b_plate_height(ii)/2));

    end
           
end

%---------- Diaphragm Deflection Plotting Example ---------

 x=(0:a_plate_width(a)/30:a_plate_width(a));
 n = length(x);
 y=(0:b_plate_height(b)/30:b_plate_height(b));
 nn = length(y);
 [xx,yy] = meshgrid(x,y);

figure;

for i=1:1:n %Rows
    
    for ii=1:1:nn %Col
        % pa(rows,col)                                                
    z(i,ii) = -1*plate_deflection(a_plate_width(a), b_plate_height(b), D(d), q0_pressure(a,b), x(i), y(ii));

    end

end
%3D plot 

 surf(xx,yy,z,'FaceAlpha',0.5);
 title('Diaphragm Deflection');
 xlabel('Diaphragm Side Length a [m]');
 ylabel('Diaphragm Side Length b [m]');
 zlabel('Diaphragm Deflection [m]');

 %---------------------- Diaphragm Stress calculations -----------------------

%  syms x y q0_pressure a_plate_width b_plate_height D
%   
%  w_plate_deflection = (q0_pressure/((pi^4)*D))*(((1/(a_plate_width^2))+(1/(b_plate_height^2)))^-2)*(sin(pi*(x/a_plate_width)))*(sin(pi*(y/b_plate_height)));
%  
%  plate_diff_dx = diff(w_plate_deflection,x,2)
%  plate_diff_dy = diff(w_plate_deflection,y,2)


%Stress corresponds to the center of the diaphragm
[x_stress, y_stress] = diaphragm_stress(q0_pressure(a,b), D(d), a_plate_width(a), b_plate_height(b),diaphragm_thickness(dt), poisson_ratio, a_plate_width(a)/2, b_plate_height(b)/2)

%Max stress occurs at the middle of the longest edge
max_edge_stress = max_diaphragm_stress(q0_pressure(a,b), diaphragm_thickness(d), a_plate_width(a), b_plate_height(b))


%================== Cantilever Beam Analysis     ==========================
%==================        Check Valve           ==========================


%------------------ Moment of Inertia ----------------------

% Beam width needs to be at least equal to or larger than the inlet
% measurements

% Side length = 4.082482904638630e-04 or 408 microns
% Beam width > 408 microns


b_beam_width_minimum = inlet_side_length + 0.00001; % 410 Microns - beam_width value + 10 microns

% Largest deflection is desirable meaning the moment of inertia must be low
% Since the moment of inertia is highly influenced by the hight of the beam
% The beam high must small and manufacturable

b_beam_width = (b_beam_width_minimum: 0.00005 : b_beam_width_minimum+0.0002);
h_beam_height = (0.00005:0.00005:0.00025); % 50 microns to 250 microns with 50 micron intervals
n = length(b_beam_width);
nn = length(h_beam_height);


figure;
for i=1:1:n
    for ii=1:1:nn
I_moment_of_inertia(i,ii) = (b_beam_width(i))*((h_beam_height(ii))^3)/12;
    end
    
     plot(b_beam_width,I_moment_of_inertia(i,:));
          hold on;
           title('Change In Beam Width VS Moment Of Inertia');
          xlabel('Beam Width [m]');
           ylabel('Moment Of Inertia [m^4]');
           grid on;
    %legend('h = 50 microns','h = 100 microns','h = 150 microns','h = 200 microns','h = 250 microns');
end
hold off;
%---------------- Max Beam Stress -------------------------

for i=1:1:n %Rows
    
    for ii=1:1:n %Col
        %Solving for distributed load
        distributed_mass(1,ii) = b_beam_width(1,ii)* h_beam_height(1,ii)*density; %Units: kg/m
        distributed_load(1,ii) = b_beam_width(1,ii)*delta_p; %Units: N/m
        %Distributed load units = Force/distance 
        %We can assume the generated pressure is applied to the inlet and
        %outlet. Since pressure units are Force/area, it is possible to
        %multiply the pressure value by the width dimension to achieve a
        %distributed load.
        
        % L(rows,col)
        L_beam_length(i,ii) = beam_length( E, I_moment_of_inertia(i,ii), cant_f, distributed_mass(1,ii) );
    end

         
end


for i=1:1:n
    
    for ii=1:1:n
   
    if(L_beam_length(i,ii) >= b_beam_width_minimum) % Checking to see if the beam length is greater than the length of the inlet. If the beam is smaller, it means no seal is created
        %                                                                 (b_w,b_h)
        M_max(i,ii) = max_moment_beam( distributed_load(1,ii), L_beam_length(i,ii));
        
        max_stress_beam(i,ii) =  max_beam_stress( h_beam_height(i)/2, M_max(i,ii), I_moment_of_inertia(i,ii) );
        
    else
        
        M_max(i,ii) = 0;
        
    end
    
    end
end


%----------------- Beam deflection ---------------------- 

n = length(L_beam_length);

row = 1;
for b_w = 1:1:n
    
    for b_h = 1:1:n
        
        
        if(L_beam_length(b_w, b_h) >= b_beam_width_minimum)
        
        x_length(1) = 0;
        for i = 2:1:10
        x_length(i) = x_length(i-1) + (L_beam_length(b_w,b_h)/9); %Creating an x axis with even distribution 
        end

        col = 1;
        %figure;
        for i = 10:-1:1

         y_def(row,col) = cantilever_deflection( E,I_moment_of_inertia(b_w,b_h),distributed_load(1,b_w), L_beam_length(b_w, b_h),x_length(i));

          col = col + 1;

        end
%         
%           plot(x_length,y_def(23,:));
%           hold on;
%           title('Beam Deflection');
%           xlabel('Beam Length [m]');
%           ylabel('Beam Deflection [m]');
%           grid on;
%           hold off;
             
          row = row+1;
          
        else
         
         y_def(row) = 0;
         %fprintf('under\n')
         row = row+1;      


        end

        end 
        
         
end
figure;
%select b_w = 5 and b_h = 3
plot(x_length,y_def(21,:));
          hold on;
          title('Beam Deflection');
          xlabel('Beam Length [m]');
          ylabel('Beam Deflection [m]');
          grid on;
          hold off;


























