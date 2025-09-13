close; clear all; clc;

T_wi = 15+273;
T = T_wi; % initial temperature in storage
T_env = 293; % Ambient temperature of air

cd_water = 4180; % Approximate specific heat of water (It changes slighlty with temperature)
irradiance = 1000;
density = 1000; % Approximately
S_B = 5.67*10e-9; % Stefan-Boltzmann
h_c_air = 25.5; % Natural convection, 25-20 C, 1 atm

k_pvc = 0.19;
k_i = 0.022;
k_floor = 0.30;
k_al = 237;
k_PU = 0.13;

em_pvc = 0.92;
em_i = 0.8;     % Approximate values based on color
em_floor = 0.1; % and intuition. Exact values
em_PU = 0.9;    % could not be found on the Internet.

d_in = 0.0085;
d_out = 0.012;
L_s_to_c = 0.5;
L_c_to_s = 1;
R_cond_tubes = log(d_out/d_in)/(2*pi*k_PU*L_s_to_c);
R_cond_tube = log(d_out/d_in)/(2*pi*k_PU*L_c_to_s);
R_conv_tubes = 1/(h_c_air*2*pi*(d_out/2)*L_s_to_c);
R_conv_tube = 1/(h_c_air*2*pi*(d_out/2)*L_c_to_s);
R_rad_tubes = 1/(S_B*2*pi*em_PU*(d_out/2)*L_s_to_c);
R_rad_tube = 1/(S_B*2*pi*em_PU*(d_out/2)*L_c_to_s);
volume_tubes = pi*L_s_to_c*(d_in/2)^2;
volume_tube = pi*L_c_to_s*(d_in/2)^2;

d_c = 0.001;     % Not mentioned
d_al = 0.002;    %
d_table = 0.006;
A_c = 1.44;
R_cond_collector = 1/(1/(d_c/(A_c*k_pvc))+1/((d_c/(A_c*k_pvc))+(d_al/(A_c*k_al))+(d_table/(A_c*k_floor))));
R_conv_collector = 1/(2*h_c_air*A_c);
R_rad_collector = 1/(S_B*A_c*(em_floor+em_pvc));
volume_collector = 0.005; % Max volume is about 10 liters. Assuming half full to not waste water.
sunlit_area_c = A_c;

D_t = 0.11;
d_wall = 0.002; % Not mentioned
d_i = 0.02;
d_lid = 0.005;  %
d_floor = 0.01; %
d_tube = d_out;
H = 0.2;
A_floor = pi*(D_t/2)^2;
A_lid = A_floor - 2*pi*(d_tube/2)^2;
R_cond_storage = 1/(1/(log((D_t/2+d_wall)/(D_t/2))/(2*pi*k_pvc*H)+log((D_t/2+d_wall+d_i)/(D_t/2+d_wall))/(2*pi*k_i*H))+1/((d_wall/(A_floor*k_pvc))+(d_i/(A_floor*k_i))+(d_floor/(A_floor*k_floor)))+1/(d_lid/(A_lid*k_pvc)));
R_conv_storage = 1/((h_c_air*2*pi*((D_t/2+d_wall+d_i)/2)*H)+(h_c_air*2*pi*A_floor)+(h_c_air*2*pi*A_lid));
R_rad_storage = 1/((S_B*2*pi*em_i*((D_t/2+d_wall+d_i)/2)*H)+(S_B*2*pi*em_floor*A_floor)+(S_B*2*pi*em_pvc*A_lid));
volume_storage = 0.002;

for section = 1:10
% tube (storage to pump)
T1 = HeatLossPerSection (T, T_env, R_cond_tubes, R_conv_tubes, R_rad_tubes, density, volume_tubes, cd_water, irradiance, 0);

% tube (from pump to hose)
T2 = HeatLossPerSection (T1, T_env, R_cond_tubes, R_conv_tubes, R_rad_tubes, density, volume_tubes, cd_water, irradiance, 0);

% collector
T3 = HeatLossPerSection (T2, T_env, R_cond_collector, R_conv_collector, R_rad_collector, density, volume_collector, cd_water, irradiance, sunlit_area_c);

% tube (hose to storage)
T4 = HeatLossPerSection (T3, T_env, R_cond_tube, R_conv_tube, R_rad_tube, density, volume_tube, cd_water, irradiance, 0);

% storage
T = HeatLossPerSection (T4, T_env, R_cond_storage, R_conv_storage, R_rad_storage, density, volume_storage, cd_water, irradiance, 0);
end

function [T] = HeatLossPerSection (T, T_env, R_cond, R_conv, R_rad, density, volume_segment, cd_water, irradiance, sunlit_area)

    % Documentaion
    %
    %
    %
    %
    %

    T_w = T;
    
    % energy_balance = irradiance*sunlit_area - (1/R_cond) * (T_outer - T_w) - ... % conduction from water to outer part of segment
    %             (1/R_conv) * (T_outer - T_env) - ... % convection to air
    %             (1/R_rad) * (T_outer^4 - T_env^4) == 0; % radiation to surroundings

    % (1/R_rad) * (T_outer^4) + (1/R_cond + 1/R_conv) * T_outer - (irradiance*sunlit_area + (T_env^4)/R_rad + T_w/R_cond + T_env/R_conv) = 0; 

    T_outer = ((irradiance*sunlit_area)+T_w/R_cond+T_env/R_conv)/(1/R_cond + 1/R_conv);

    Q_cond = (1/R_cond) * (T_w - T_outer); % solve for heat conducted from water to outer surface

    dT = Q_cond / (density * volume_segment * cd_water); % change in temperature over segment

    % new temperature for next segment
    T = T - dT;
end