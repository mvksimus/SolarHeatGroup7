
T = T_wi; % initial temperature at pump

 % tube (from pump to hose)
    for section = 1:n % loop for # of segments
        T_w = T;
        energy_balance = - (1/R_cond) * (T_tube_outer - T_w) - ... % conduction from water to outer part of tube
                (1/R_conv) * (T_tube_outer - T_env) - ... % convection to air
                (1/R_rad) * (T_tube_outer^4 - T_env^4) == 0; % radiation to surroundings

        T_tube_outer = solve(energy_balance, T_tube_outer); % solve for the temperature at the outer part of the tube

        Q_cond = (1/R_cond) * (T_tube_outer - T_w); % solve for heat conducted from water to tube

        dT = Q_cond / (massSegment * cWater); % change in temperature over segment

        % new temperature for next segment
        T = T + dT; 
    end

