T_wi = 15;
T = T_wi; % initial temperature in storage
T_env = 293; % Ambient temperature of air
cd_water = 4.18; % Approximate specific heat of water (It changes slighlty with temperature)
irradiance = 1000;
density = 1000; % Approximately
S_B = 5.67*10e-9; % Stefan-Boltzmann

% tube (storage to pump)

%for section = 1:n % loop for # of segments
%T_tube1_inlet = T_wi;

%        energy_balance = - (1/R_cond) * (T_tube_outer - T_w) - ... % conduction from water to outer part of tube
%                (1/R_conv) * (T_tube_outer - T_env) - ...  % convection to air
%               (1/R_rad) * (T_tube_outer^4 - T_env^4) == 0; % radiation to surroundings

%        T_tube_outer = solve(energy_balance, T_tube_outer); % solve for the temperature at the outer part of the tube

%        Q_cond = (1/R_cond) * (T_tube_outer - T_w); % solve for heat conducted from water to tube

%        dT = Q_cond / (massSegment * cWater); % change in temperature over segment

%        T_tube1_outlet = T_tube1_inlet + dT;  % new temperature for next segment
%end

 %pump segment?

 % tube (from pump to hose)
 %   for section = 1:n % loop for # of segments
 %       T_w = T;
 %       energy_balance = - (1/R_cond) * (T_tube_outer - T_w) - ... % conduction from water to outer part of tube
 %               (1/R_conv) * (T_tube_outer - T_env) - ... % convection to air
 %               (1/R_rad) * (T_tube_outer^4 - T_env^4) == 0; % radiation to surroundings

  %      T_tube_outer = solve(energy_balance, T_tube_outer); % solve for the temperature at the outer part of the tube

%        Q_cond = (1/R_cond) * (T_tube_outer - T_w); % solve for heat conducted from water to tube

 %       dT = Q_cond / (massSegment * cWater); % change in temperature over segment

       % new temperature for next segment
  %      T = T + dT; 
   % end

% hose (tube connection to collector)

    %for section = 1:n % loop for # of segments
    %       T_w = T;
    %      energy_balance_hose1 = irradiance * areaHose1 - (1/R_cond_hose) * (T_hose1_outer - T_w) - ... % conduction from water to outer part of tube
    %                (1/R_conv_hose) * (T_hose1_outer - T_env) - ... % convection to air
    %                (1/R_rad_hose) * (T_hose1_outer^4 - T_env^4) == 0; % radiation to surroundings
    
     %       T_hose_outer = solve(energy_balance_hose1, T_hose1_outer); % solve for the temperature at the outer part of the hose
    
      %      Q_cond_hose1 = (1/R_cond_hose) * (T_hose1_outer - T_w); % solve for heat conducted from water to hose
    
      %      dT = Q_cond_hose1 / (massSegmentHose1 * cWater); % change in temperature over segment
    
            % new temperature for next segment
       %     T = T + dT; 
    %end

% collector

    %for section = 1:n % loop for # of segments
     %       T_w = T;
     %       energy_balance_collectorTop = irradiance * areaCollector - (1/R_cond_coll) * (T_coll_outer - T_w) - ... % conduction from water to outer part of tube
     %               (1/R_conv_coll) * (T_coll_outer - T_env) - ... % convection to air
     %               (1/R_rad_coll) * (T_coll_outer^4 - T_env^4) == 0; % radiation to surroundings

      %      energy_balance_collectorBottom = - (1/R_cond_coll) * (T_coll_outer - T_w) - ... % conduction from water to bottom of the table
       %             (1/R_conv_coll) * (T_coll_outer - T_alum) - ... % convection to air
       %             (1/R_rad_coll) * (T_coll_outer^4 - T_env^4) == 0; % radiation to surroundings

       %     energy_balance_collector = energy_balance_collectorBottom + energy_balance_collectorTop;
    
        %    T_coll_outer = solve(energy_balance_collector, T_coll_outer); % solve for the temperature at the outer part of the collector
    
         %   Q_cond_coll = (1/R_cond_coll) * (T_coll_outer - T_w); % solve for heat conducted from water to hose
    
          %  dT = Q_cond_coll / (massSegmentColl1 * cWater); % change in temperature over segment
    
            % new temperature for next segment
           % T = T + dT; 
    %end

% hose (collector to tube connection)

% tube (hose to storage)

% storage

h_air = 25.5; % Natural convection, 25-20 C, 1 atm

k_pvc = 0.19;
k_i = 0.022;
k_floor = 0.30;
k_al = ;
k_PU = 0.13;

em_pvc = ;
em_i = ;
em_wood = ;
em_PU = ;

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

R_cond_hose = ;
R_conv_hose = ;
R_rad_hose = ;
volume_hose = ;

R_cond_collector = ;
R_conv_collector = ;
R_rad_collector = ;
volume_collector = ;

D_t = 0.11;
d_wall = 0.002;
d_i = 0.02;
d_lid = 0.005;
d_floor = 0.01;
H = 0.2;
A_floor = pi*(D_t/2)^2;
A_lid = A_floor - 2*pi*(d_tube/2)^2;
R_cond_storage = 1/(1/(log((D_t/2+d_wall)/(D_t/2))/(2*pi*k_pvc*H)+log((D_t/2+d_wall+d_i)/(D_t/2+d_wall))/(2*pi*k_i*H))+1/((d_wall/(A_floor*k_pvc))+(d_i/(A_floor*k_i))+(d_wood/(A_floor*k_floor)))+1/(d_lid/(A_lid*k*pvc)));
R_conv_storage = 1/(h_c_air*2*pi*((D_t/2+d_wall+d_i)/2)*H)+1/(h_c_air*2*pi*A_floor)+1/(h_c_air*2*pi*A_lid);
R_rad_storage = 1/(S_B*2*pi*em_i*((D_t/2+d_wall+d_i)/2)*H)+1/(S_B*2*pi*em_wood*A_floor)+1/(S_B*2*pi*em_pvc*A_lid);
volume_storage = 2*10^-3;

for section = 1:10
% tube (storage to pump)
T1 = HeatLossPerSection (T, T_env, R_cond_tube, R_conv_tube, R_rad_tube, mass_tubes, cd_water, irradiance, 0);
% pump segment?
% tube (from pump to hose)
T2 = HeatLossPerSection (T1, T_env, R_cond_tube, R_conv_tube, R_rad_tube, density, volume_tubes, cd_water, irradiance, 0);
% hose (tube connection to collector)
T3 = HeatLossPerSection (T2, T_env, R_cond_hose, R_conv_hose, R_rad_hose, density, volume_hose, cd_water, irradiance, sunlit_area_h);
% collector
T4 = HeatLossPerSection (T3, T_env, R_cond_collector, R_conv_collector, R_rad_collector, density, volume_collector, cd_water, irradiance, sunlit_area_c);
% hose (collector to tube connection)
T5 = HeatLossPerSection (T4, T_env, R_cond_hose, R_conv_hose, R_rad_hose, density, volume_hose, cd_water, irradiance, sunlit_area_h);
% tube (hose to storage)
T6 = HeatLossPerSection (T5, T_env, R_cond_tube, R_conv_tube, R_rad_tube, density, volume_tube, cd_water, irradiance, 0);
% storage
T = HeatLossPerSection (T6, T_env, R_cond_storage, R_conv_storage, R_rad_storage, density, volume_storage, cd_water, irradiance, 0);
end

function [T] = HeatLossPerSection (T, T_env, R_cond, R_conv, R_rad, density, volume_segment, cd_water, irradiance, sunlit_area)

    %
    %
    %
    %
    %
    %

    T_w = T;
    
    energy_balance = irradiance*sunlit_area - (1/R_cond) * (T_outer - T_w) - ... % conduction from water to outer part of segment
                (1/R_conv) * (T_outer - T_env) - ... % convection to air
                (1/R_rad) * (T_outer^4 - T_env^4) == 0; % radiation to surroundings

    T_outer = solve(energy_balance, T_outer); % solve for the temperature at the outer part of the tube

    Q_cond = (1/R_cond) * (T_w - T_outer); % solve for heat conducted from water to outer surface

    dT = Q_cond / (density * volume_segment * cd_water); % change in temperature over segment

    % new temperature for next segment
    T = T - dT;
end