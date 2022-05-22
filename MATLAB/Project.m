% Copyright (c) 2022 Vrund Patel, USA
%
%% ASE 366L: Applied Orbital Mechanics
% Course Project - SPRING 2022

% Vrund Patel

format long

%% Constants:

MU_EARTH = 398600.4415; % km^3/s^2
radiusEarth = 6378.1363; % km
J2 = 0.0010826267;
J3 = -0.0000025327;
Cd = 2.0;
Cr = 1.5; 
AreaMassRatio = 0.01 * 1e-6; % km^2/kg
MU_SUN = 1.327e11; % km^3/s^2
MU_MOON = 3903; % km^3/s^2
w_earth = 2 * pi / 86164; %rad/s

%% Task 1:

% tp = March 1, 2020, 12:00:00 UTC

% Assuming UT1 = UTC
JD_UT1 = cal_to_JD(2020, 3, 1, 12, 0, 0);

timespan = 0 : 10 * 60 : 5 * 24 * 3600; % timepsan in seconds - 10 min increments up to 5 days
MJD_TAI = JD_UT1 - 2400000.5;
ode_options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

outTimes = MJD_TAI + timespan / 86400; % days
outTimes_mins = outTimes * 24 * 60; % minutes

%% LEO Orbit:

Y0_LEO = COEstoRV(6763, 0.001, 50, 0, 0, 0, 0, 0, 0, MU_EARTH);

% Two-body only:
[T1_LEO, Y1_LEO] = ode45(@orbitPropagator, timespan, Y0_LEO, ode_options, MJD_TAI, MU_EARTH, 0, 0, 0, 0, 0, 0, 0, 0, 0);

% Two-body with J2:
[T2_LEO, Y2_LEO] = ode45(@orbitPropagator, timespan, Y0_LEO, ode_options, MJD_TAI, MU_EARTH, 0, 0, J2, 0, radiusEarth, 0, 0, 0, 0);

% Two-body with J3:
[T3_LEO, Y3_LEO] = ode45(@orbitPropagator, timespan, Y0_LEO, ode_options, MJD_TAI, MU_EARTH, 0, 0, 0, J3, radiusEarth, 0, 0, 0, 0);

% Two-body with Drag:
[T4_LEO, Y4_LEO] = ode45(@orbitPropagator, timespan, Y0_LEO, ode_options, MJD_TAI, MU_EARTH, 0, 0, 0, 0, radiusEarth, w_earth, Cd, 0, AreaMassRatio);

% Two-body with Sun third-body:
[T5_LEO, Y5_LEO] = ode45(@orbitPropagator, timespan, Y0_LEO, ode_options, MJD_TAI, MU_EARTH, MU_SUN, 0, 0, 0, 0, 0, 0, 0, 0);

% Two-body with Moon third-body:
[T6_LEO, Y6_LEO] = ode45(@orbitPropagator, timespan, Y0_LEO, ode_options, MJD_TAI, MU_EARTH, 0, MU_MOON, 0, 0, 0, 0, 0, 0, 0);

% Two-body with Solar Radiation Pressure:
[T7_LEO, Y7_LEO] = ode45(@orbitPropagator, timespan, Y0_LEO, ode_options, MJD_TAI, MU_EARTH, 0, 0, 0, 0, 0, 0, 0, Cr, AreaMassRatio);

twobody_J2_LEO = zeros(length(timespan), 1);
twobody_J3_LEO = zeros(length(timespan), 1);
twobody_Drag_LEO = zeros(length(timespan), 1);
twobody_Sun_LEO = zeros(length(timespan), 1);
twobody_Moon_LEO = zeros(length(timespan), 1);
twobody_SRP_LEO = zeros(length(timespan), 1);

for i = 1: length(timespan)
    
    twobody_J2_LEO(i) = norm( Y2_LEO(i, 1:3) - Y1_LEO(i, 1:3) );
    twobody_J3_LEO(i) = norm( Y3_LEO(i, 1:3) - Y1_LEO(i, 1:3) );
    twobody_Drag_LEO(i) = norm( Y4_LEO(i, 1:3) - Y1_LEO(i, 1:3) );
    twobody_Sun_LEO(i) = norm( Y5_LEO(i, 1:3) - Y1_LEO(i, 1:3) );
    twobody_Moon_LEO(i) = norm( Y6_LEO(i, 1:3) - Y1_LEO(i, 1:3) );
    twobody_SRP_LEO(i) = norm( Y7_LEO(i, 1:3) - Y1_LEO(i, 1:3) );

end

figure(1);
hold on
abscissa = outTimes_mins - outTimes_mins(1);
set(gca, 'YScale', 'log')
grid on

plot(abscissa, twobody_J2_LEO, 'r', 'DisplayName', 'Two-Body + J2')
plot(abscissa, twobody_J3_LEO, 'b', 'DisplayName', 'Two-Body + J3')
plot(abscissa, twobody_Drag_LEO, 'g', 'DisplayName', 'Two-Body + Drag')
plot(abscissa, twobody_Sun_LEO, 'c', 'DisplayName', 'Two-Body + Sun')
plot(abscissa, twobody_Moon_LEO, 'k', 'DisplayName', 'Two-Body + Moon')
plot(abscissa, twobody_SRP_LEO, 'm', 'DisplayName', 'Two-Body + SRP')

xlabel("Minutes Since Epoch")
ylabel("Error in Distance [km]")
%xticks([0 1440 2880 4320 5760 7200])
legend

%% MEO Orbit:

Y0_MEO = COEstoRV(26560, 0.001, 55, 0, 0, 0, 0, 0, 0, MU_EARTH);

% Two-body only:
[T1_MEO, Y1_MEO] = ode45(@orbitPropagator, timespan, Y0_MEO, ode_options, MJD_TAI, MU_EARTH, 0, 0, 0, 0, 0, 0, 0, 0, 0);

% Two-body with J2:
[T2_MEO, Y2_MEO] = ode45(@orbitPropagator, timespan, Y0_MEO, ode_options, MJD_TAI, MU_EARTH, 0, 0, J2, 0, radiusEarth, 0, 0, 0, 0);

% Two-body with J3:
[T3_MEO, Y3_MEO] = ode45(@orbitPropagator, timespan, Y0_MEO, ode_options, MJD_TAI, MU_EARTH, 0, 0, 0, J3, radiusEarth, 0, 0, 0, 0);

% Two-body with Drag:
[T4_MEO, Y4_MEO] = ode45(@orbitPropagator, timespan, Y0_MEO, ode_options, MJD_TAI, MU_EARTH, 0, 0, 0, 0, radiusEarth, w_earth, Cd, 0, AreaMassRatio);

% Two-body with Sun third-body:
[T5_MEO, Y5_MEO] = ode45(@orbitPropagator, timespan, Y0_MEO, ode_options, MJD_TAI, MU_EARTH, MU_SUN, 0, 0, 0, 0, 0, 0, 0, 0);

% Two-body with Moon third-body:
[T6_MEO, Y6_MEO] = ode45(@orbitPropagator, timespan, Y0_MEO, ode_options, MJD_TAI, MU_EARTH, 0, MU_MOON, 0, 0, 0, 0, 0, 0, 0);

% Two-body with Solar Radiation Pressure:
[T7_MEO, Y7_MEO] = ode45(@orbitPropagator, timespan, Y0_MEO, ode_options, MJD_TAI, MU_EARTH, 0, 0, 0, 0, 0, 0, 0, Cr, AreaMassRatio);

twobody_J2_MEO = zeros(length(timespan), 1);
twobody_J3_MEO = zeros(length(timespan), 1);
twobody_Drag_MEO = zeros(length(timespan), 1);
twobody_Sun_MEO = zeros(length(timespan), 1);
twobody_Moon_MEO = zeros(length(timespan), 1);
twobody_SRP_MEO = zeros(length(timespan), 1);

for i = 1: length(timespan)
    
    twobody_J2_MEO(i) = norm( Y2_MEO(i, 1:3) - Y1_MEO(i, 1:3) );
    twobody_J3_MEO(i) = norm( Y3_MEO(i, 1:3) - Y1_MEO(i, 1:3) );
    twobody_Drag_MEO(i) = norm( Y4_MEO(i, 1:3) - Y1_MEO(i, 1:3) );
    twobody_Sun_MEO(i) = norm( Y5_MEO(i, 1:3) - Y1_MEO(i, 1:3) );
    twobody_Moon_MEO(i) = norm( Y6_MEO(i, 1:3) - Y1_MEO(i, 1:3) );
    twobody_SRP_MEO(i) = norm( Y7_MEO(i, 1:3) - Y1_MEO(i, 1:3) );

end

figure(2)
hold on
abscissa = outTimes_mins - outTimes_mins(1);
set(gca, 'YScale', 'log')
grid on

plot(abscissa, twobody_J2_MEO, 'r', 'DisplayName', 'Two-Body + J2')
plot(abscissa, twobody_J3_MEO, 'b', 'DisplayName', 'Two-Body + J3')
plot(abscissa, twobody_Drag_MEO, 'g', 'DisplayName', 'Two-Body + Drag')
plot(abscissa, twobody_Sun_MEO, 'c', 'DisplayName', 'Two-Body + Sun')
plot(abscissa, twobody_Moon_MEO, 'k', 'DisplayName', 'Two-Body + Moon')
plot(abscissa, twobody_SRP_MEO, 'm', 'DisplayName', 'Two-Body + SRP')

xlabel("Minutes Since Epoch")
ylabel("Error in Distance [km]")
%xticks([0 1440 2880 4320 5760 7200])
legend

%% GEO Orbit:

Y0_GEO = COEstoRV(42164, 0.01, 0.5, -120, 0, 0, 0, 0, 0, MU_EARTH);

% Two-body only:
[T1_GEO, Y1_GEO] = ode45(@orbitPropagator, timespan, Y0_GEO, ode_options, MJD_TAI, MU_EARTH, 0, 0, 0, 0, 0, 0, 0, 0, 0);

% Two-body with J2:
[T2_GEO, Y2_GEO] = ode45(@orbitPropagator, timespan, Y0_GEO, ode_options, MJD_TAI, MU_EARTH, 0, 0, J2, 0, radiusEarth, 0, 0, 0, 0);

% Two-body with J3:
[T3_GEO, Y3_GEO] = ode45(@orbitPropagator, timespan, Y0_GEO, ode_options, MJD_TAI, MU_EARTH, 0, 0, 0, J3, radiusEarth, 0, 0, 0, 0);

% Two-body with Drag:
[T4_GEO, Y4_GEO] = ode45(@orbitPropagator, timespan, Y0_GEO, ode_options, MJD_TAI, MU_EARTH, 0, 0, 0, 0, radiusEarth, w_earth, Cd, 0, AreaMassRatio);

% Two-body with Sun third-body:
[T5_GEO, Y5_GEO] = ode45(@orbitPropagator, timespan, Y0_GEO, ode_options, MJD_TAI, MU_EARTH, MU_SUN, 0, 0, 0, 0, 0, 0, 0, 0);

% Two-body with Moon third-body:
[T6_GEO, Y6_GEO] = ode45(@orbitPropagator, timespan, Y0_GEO, ode_options, MJD_TAI, MU_EARTH, 0, MU_MOON, 0, 0, 0, 0, 0, 0, 0);

% Two-body with Solar Radiation Pressure:
[T7_GEO, Y7_GEO] = ode45(@orbitPropagator, timespan, Y0_GEO, ode_options, MJD_TAI, MU_EARTH, 0, 0, 0, 0, 0, 0, 0, Cr, AreaMassRatio);

twobody_J2_GEO = zeros(length(timespan), 1);
twobody_J3_GEO = zeros(length(timespan), 1);
twobody_Drag_GEO = zeros(length(timespan), 1);
twobody_Sun_GEO = zeros(length(timespan), 1);
twobody_Moon_GEO = zeros(length(timespan), 1);
twobody_SRP_GEO = zeros(length(timespan), 1);

for i = 1: length(timespan)
    
    twobody_J2_GEO(i) = norm( Y2_GEO(i, 1:3) - Y1_GEO(i, 1:3) );
    twobody_J3_GEO(i) = norm( Y3_GEO(i, 1:3) - Y1_GEO(i, 1:3) );
    twobody_Drag_GEO(i) = norm( Y4_GEO(i, 1:3) - Y1_GEO(i, 1:3) );
    twobody_Sun_GEO(i) = norm( Y5_GEO(i, 1:3) - Y1_GEO(i, 1:3) );
    twobody_Moon_GEO(i) = norm( Y6_GEO(i, 1:3) - Y1_GEO(i, 1:3) );
    twobody_SRP_GEO(i) = norm( Y7_GEO(i, 1:3) - Y1_GEO(i, 1:3) );

end

figure(3)
hold on
abscissa = outTimes_mins - outTimes_mins(1);
set(gca, 'YScale', 'log')
grid on

plot(abscissa, twobody_J2_GEO, 'r', 'DisplayName', 'Two-Body + J2')
plot(abscissa, twobody_J3_GEO, 'b', 'DisplayName', 'Two-Body + J3')
plot(abscissa, twobody_Drag_GEO, 'g', 'DisplayName', 'Two-Body + Drag')
plot(abscissa, twobody_Sun_GEO, 'c', 'DisplayName', 'Two-Body + Sun')
plot(abscissa, twobody_Moon_GEO, 'k', 'DisplayName', 'Two-Body + Moon')
plot(abscissa, twobody_SRP_GEO, 'm', 'DisplayName', 'Two-Body + SRP')

xlabel("Minutes Since Epoch")
ylabel("Error in Distance [km]")
%xticks([0 1440 2880 4320 5760 7200])
legend

%% Molniya Orbit:

Y0_Mol = COEstoRV(26000, 0.72, 75, 90, -90, 0, 0, 0, 0, MU_EARTH);

% Two-body only:
[T1_Mol, Y1_Mol] = ode45(@orbitPropagator, timespan, Y0_Mol, ode_options, MJD_TAI, MU_EARTH, 0, 0, 0, 0, 0, 0, 0, 0, 0);

% Two-body with J2:
[T2_Mol, Y2_Mol] = ode45(@orbitPropagator, timespan, Y0_Mol, ode_options, MJD_TAI, MU_EARTH, 0, 0, J2, 0, radiusEarth, 0, 0, 0, 0);

% Two-body with J3:
[T3_Mol, Y3_Mol] = ode45(@orbitPropagator, timespan, Y0_Mol, ode_options, MJD_TAI, MU_EARTH, 0, 0, 0, J3, radiusEarth, 0, 0, 0, 0);

% Two-body with Drag:
[T4_Mol, Y4_Mol] = ode45(@orbitPropagator, timespan, Y0_Mol, ode_options, MJD_TAI, MU_EARTH, 0, 0, 0, 0, radiusEarth, w_earth, Cd, 0, AreaMassRatio);

% Two-body with Sun third-body:
[T5_Mol, Y5_Mol] = ode45(@orbitPropagator, timespan, Y0_Mol, ode_options, MJD_TAI, MU_EARTH, MU_SUN, 0, 0, 0, 0, 0, 0, 0, 0);

% Two-body with Moon third-body:
[T6_Mol, Y6_Mol] = ode45(@orbitPropagator, timespan, Y0_Mol, ode_options, MJD_TAI, MU_EARTH, 0, MU_MOON, 0, 0, 0, 0, 0, 0, 0);

% Two-body with Solar Radiation Pressure:
[T7_Mol, Y7_Mol] = ode45(@orbitPropagator, timespan, Y0_Mol, ode_options, MJD_TAI, MU_EARTH, 0, 0, 0, 0, 0, 0, 0, Cr, AreaMassRatio);

twobody_J2_Mol = zeros(length(timespan), 1);
twobody_J3_Mol = zeros(length(timespan), 1);
twobody_Drag_Mol = zeros(length(timespan), 1);
twobody_Sun_Mol = zeros(length(timespan), 1);
twobody_Moon_Mol = zeros(length(timespan), 1);
twobody_SRP_Mol = zeros(length(timespan), 1);

for i = 1: length(timespan)
    
    twobody_J2_Mol(i) = norm( Y2_Mol(i, 1:3) - Y1_Mol(i, 1:3) );
    twobody_J3_Mol(i) = norm( Y3_Mol(i, 1:3) - Y1_Mol(i, 1:3) );
    twobody_Drag_Mol(i) = norm( Y4_Mol(i, 1:3) - Y1_Mol(i, 1:3) );
    twobody_Sun_Mol(i) = norm( Y5_Mol(i, 1:3) - Y1_Mol(i, 1:3) );
    twobody_Moon_Mol(i) = norm( Y6_Mol(i, 1:3) - Y1_Mol(i, 1:3) );
    twobody_SRP_Mol(i) = norm( Y7_Mol(i, 1:3) - Y1_Mol(i, 1:3) );

end

figure(4)
hold on
abscissa = outTimes_mins - outTimes_mins(1);
set(gca, 'YScale', 'log')
grid on

plot(abscissa, twobody_J2_Mol, 'r', 'DisplayName', 'Two-Body + J2')
plot(abscissa, twobody_J3_Mol, 'b', 'DisplayName', 'Two-Body + J3')
plot(abscissa, twobody_Drag_Mol, 'g', 'DisplayName', 'Two-Body + Drag')
plot(abscissa, twobody_Sun_Mol, 'c', 'DisplayName', 'Two-Body + Sun')
plot(abscissa, twobody_Moon_Mol, 'k', 'DisplayName', 'Two-Body + Moon')
plot(abscissa, twobody_SRP_Mol, 'm', 'DisplayName', 'Two-Body + SRP')

xlabel("Minutes Since Epoch")
ylabel("Error in Distance [km]")
%xticks([0 1440 2880 4320 5760 7200])
legend

%% Task 2:

% Burnout coordinates in degrees
lat1 = 30.3;
long1 = -120.6;
lat2 = 30.3;
long2 = -97.7;

% Known parameters at burnout
a = 6500; % km
e = 0.001;
nu1 = 20; % deg
w_earth_deg = 360 / 86164; % deg/s
theta_ERA_0 = 0;
T = 2 * pi * sqrt(a^3 / MU_EARTH); % period in seconds
n = 10; % number of orbit periods before flyby

% Algorithm - NASA Technical Note D-233 (1960)
t_nu1 = TimeFromPerigee(nu1, T, e);
long2e = long2 + n * w_earth_deg * T;

delta_long_1_2e = zeros(4, 1);
nu2e = zeros(4, 1);
t_nu2e = zeros(4, 1);
psi1 = zeros(4, 1);
inc = zeros(4, 1);
w = zeros(4, 1);

t_nu2e_nu1 = ((T / 360) * (long2e - long1));

for i = 1: 4

    if i == 1

        delta_long_1_2e(i) = long2e - long1 + w_earth_deg * t_nu2e_nu1;
        nu2e(i) = myacosd( sind(lat2) * sind(lat1) + cosd(lat2) * cosd(lat1) * cosd(delta_long_1_2e(i)) ) + nu1;
        t_nu2e(i) = TimeFromPerigee(nu2e(i), T, e);

        psi1(i) = myasind( (sind(delta_long_1_2e(i)) * cosd(lat2)) / sind(nu2e(i) - nu1) );

        inc(i) = myacosd(sind(psi1(i)) * cosd(lat1));

        w(i) = myasind(sind(lat1) / sind(inc(i))) - nu1;

        meanMotion = sqrt(MU_EARTH / a^3) * 180 / pi;
        p = a * (1 - e^2);
        delta_argP = ( (3 * meanMotion * radiusEarth^2 * J2 * (4 - 5 * sind(inc(i))^2)) / (4 * p^2) ) * (n * T + t_nu2e_nu1);
        delta_RAAN = ( (-3 * meanMotion * radiusEarth^2 * J2 * cosd(inc(i))) / (2 * p^2) ) * (n * T + t_nu2e_nu1);


        delta_lat2 = (sind(inc(i)) * cosd(w(i) + nu2e(i)) * delta_argP) / (cosd(lat2));
        delta_long2 = ( (cosd(inc(i)) * secd(w(i) + nu2e(i))^2 * delta_argP) / (1 + cosd(inc(i))^2 * tand(w(i) + nu2e(i))^2 ) ) + delta_RAAN;

        lat2 = lat2 - delta_lat2;
        long2 = long2 - delta_long2;
        long2e = long2 + n * w_earth_deg * T;

        t_nu2e_nu1 = t_nu2e(i) - t_nu1;

    else

        delta_long_1_2e(i) = long2e - long1 + w_earth_deg * t_nu2e_nu1;
        nu2e(i) = myacosd( sind(lat2) * sind(lat1) + cosd(lat2) * cosd(lat1) * cosd(delta_long_1_2e(i)) ) + nu1;
        t_nu2e(i) = TimeFromPerigee(nu2e(i), T, e);

        psi1(i) = myasind( (sind(delta_long_1_2e(i)) * cosd(lat2)) / sind(nu2e(i) - nu1) );

        inc(i) = myacosd(sind(psi1(i)) * cosd(lat1));

        t_nu2e_nu1 = t_nu2e(i) - t_nu1;

        w(i) = myasind(sind(lat1) / sind(inc(i))) - nu1;
    
    end

end

delta_longN1 = atand( sind(lat1) * tand(psi1(4)) );
longNref = long1 - delta_longN1;
RAAN = longNref + theta_ERA_0;

argP = w(4);
incl = inc(4);
true_anom = nu1;

timespanGT = 0: 30 : n * T + (0.3 * T - 240);

Y0 = COEstoRV(a, e, incl, RAAN, argP, true_anom, 0, 0, 0, MU_EARTH); % state vectors at burnout
[T_J2, Y_J2] = ode45(@orbitPropagator, timespanGT, Y0, ode_options, MJD_TAI, MU_EARTH, 0, 0, J2, 0, radiusEarth, 0, 0, 0, 0);
[T_full, Y_full] = ode45(@orbitPropagator, timespanGT, Y0, ode_options, MJD_TAI, MU_EARTH, MU_SUN, MU_MOON, J2, J3, radiusEarth, w_earth, Cd, Cr, AreaMassRatio);

theta_ERA = zeros(length(timespanGT), 1);

lat_J2 = zeros(length(timespanGT), 1);
long_J2 = zeros(length(timespanGT), 1);
r_ITRF_J2 = zeros(length(timespanGT), 3);
r_J2 = zeros(length(timespanGT), 1);

lat_full = zeros(length(timespanGT), 1);
long_full = zeros(length(timespanGT), 1);
r_ITRF_full = zeros(length(timespanGT), 3);
r_full = zeros(length(timespanGT), 1);

for i = 1: length(timespanGT)

    theta_ERA(i) = theta_ERA_0 + w_earth_deg * (timespanGT(i) - 0);
    
    % Coordinates for J2 propagation:
    r_ITRF_J2(i, :) =  r3Rotation(theta_ERA(i)) * transpose(Y_J2(i, 1:3));
    r_J2(i) = norm(r_ITRF_J2(i, :));

    lat_J2(i) = asind(r_ITRF_J2(i, 3) / r_J2(i));
    long_J2(i) = atan2d(r_ITRF_J2(i, 2), r_ITRF_J2(i, 1));

    % Coordinates for full-fidelity propagation:
    r_ITRF_full(i, :) =  r3Rotation(theta_ERA(i)) * transpose(Y_full(i, 1:3));
    r_full(i) = norm(r_ITRF_full(i, :));

    lat_full(i) = asind(r_ITRF_full(i, 3) / r_full(i));
    long_full(i) = atan2d(r_ITRF_full(i, 2), r_ITRF_full(i, 1));

end

load earth_coastline.mat

figure(5)

plot(earth_coastline(:, 1), earth_coastline(:, 2), 'k')
hold on
plot(long_J2(:), lat_J2(:), '.b')
plot(long_J2(1), lat_J2(1), 'k.', 'MarkerSize', 12)
plot(long_J2(end), lat_J2(end), 'r.', 'MarkerSize', 12)
grid on
axis equal
xlim([-180, 180])
ylim([-90, 90])
xticks([-180 -150 -120 -90 -60 -30 0 30 60 90 120 150 180])
yticks([-90 -60 -30 0 30 60 90])
xlabel("Longitude [deg]")
ylabel("Latitude [deg]")

figure(6)

plot(earth_coastline(:, 1), earth_coastline(:, 2), 'k')
hold on
plot(long_full(:), lat_full(:), '.b')
plot(long_full(1), lat_full(1), 'k.', 'MarkerSize', 12)
plot(long_full(end), lat_full(end), 'r.', 'MarkerSize', 12)
grid on
axis equal
xlim([-180, 180])
ylim([-90, 90])
xticks([-180 -150 -120 -90 -60 -30 0 30 60 90 120 150 180])
yticks([-90 -60 -30 0 30 60 90])
xlabel("Longitude [deg]")
ylabel("Latitude [deg]")


%% Functions:

% Local Function
function t = TimeFromPerigee(nu, T, e)
    
    nu = nu * pi / 180;

    E = 2 * atan2( sqrt(1 - e) * tan(nu / 2), sqrt(1 + e) );

    t = (T / (2 * pi)) * (E - e * sin(E));

end
