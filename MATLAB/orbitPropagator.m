% Copyright (c) 2022 Vrund Patel, USA
%

function dY = orbitPropagator(t, Y, MJD_TAI, MU_EARTH, MU_SUN, MU_MOON, J2, J3, radiusEarth, w_earth, Cd, Cr, AreaMassRatio)

    r_vec = Y(1:3);
    v_vec = Y(4:6);
    
    currentTime = MJD_TAI + t / 86400;
    
    % Assuming TAI = UT1
    JD_UT1 = currentTime;
    
    sun_pos_vec = Sun(JD_UT1);
    moon_pos_vec = Moon(JD_UT1);
    
    two_body_acc = TwoBody(r_vec, MU_EARTH);
    J2_acc = getJ2(r_vec, J2, radiusEarth, MU_EARTH);
    J3_acc = getJ3(r_vec, J3, radiusEarth, MU_EARTH);
    third_body_sun_acc = ThirdBodySun(r_vec, sun_pos_vec, MU_SUN);
    third_body_moon_acc = ThirdBodyMoon(r_vec, moon_pos_vec, MU_MOON);
    Drag_acc = getDrag(r_vec, v_vec, w_earth, radiusEarth, Cd, AreaMassRatio);
    SRP_acc = getSRP(r_vec, sun_pos_vec, Cr, AreaMassRatio);
    
    dY = zeros(size(Y));
    dY(1:3) = v_vec;
    
    dY(4:6) = two_body_acc + J2_acc + J3_acc + third_body_sun_acc + third_body_moon_acc + Drag_acc + SRP_acc;

end


function acc = TwoBody(r_vec, MU_EARTH)

    r = norm(r_vec);
    
    acc = (-MU_EARTH / r^3) * r_vec;

end

function acc = getJ2(r_vec, J2, radiusEarth, MU)

    r = norm(r_vec);
    ri = r_vec(1);
    rj = r_vec(2);
    rk = r_vec(3);
    
    ai = ( -3 * J2 * MU * radiusEarth^2 * ri / (2 * r^5) ) * (1 - 5 * rk^2 / r^2); 
    aj = ( -3 * J2 * MU * radiusEarth^2 * rj / (2 * r^5) ) * (1 - 5 * rk^2 / r^2);
    ak = ( -3 * J2 * MU * radiusEarth^2 * rk / (2 * r^5) ) * (3 - 5 * rk^2 / r^2);

    acc = [ai; aj; ak];

end

function acc = getJ3(r_vec, J3, radiusEarth, MU)

    r = norm(r_vec);
    ri = r_vec(1);
    rj = r_vec(2);
    rk = r_vec(3);

    ai = ( -5 * J3 * MU * radiusEarth^3 * ri / (2 * r^7) ) * (3 * rk - 7 * rk^3 / r^2);
    aj = ( -5 * J3 * MU * radiusEarth^3 * rj / (2 * r^7) ) * (3 * rk - 7 * rk^3 / r^2);
    ak = ( -5 * J3 * MU * radiusEarth^3 / (2 * r^7) ) * (6 * rk^2 - 7 * rk^4 / r^2 - 0.6 * r^2);

    acc = [ai; aj; ak];

end

function acc = ThirdBodySun(r_vec, sun_pos_vec, MU_SUN)

    sun_pos = norm(sun_pos_vec);
    
    sat_sun_pos_vec = sun_pos_vec - r_vec;
    sat_sun_pos = norm(sat_sun_pos_vec);
    
    acc = MU_SUN * ((sat_sun_pos_vec / sat_sun_pos^3) - (sun_pos_vec / sun_pos^3));

end

function acc = ThirdBodyMoon(r_vec, moon_pos_vec, MU_MOON)

    moon_pos = norm(moon_pos_vec);
    
    sat_moon_pos_vec = moon_pos_vec - r_vec;
    sat_moon_pos = norm(sat_moon_pos_vec);
    
    acc = MU_MOON * ((sat_moon_pos_vec / sat_moon_pos^3) - (moon_pos_vec / moon_pos^3));

end

function acc = getDrag(r_vec, v_vec, w_earth, radiusEarth, Cd, AreaMassRatio)

    v_rel_vec = v_vec - cross([0; 0; w_earth], r_vec);
    v_rel = norm(v_rel_vec);
    
    r = norm(r_vec);
    altitude = r - radiusEarth;
    params = getDensityParams(altitude);
    rho0 = params(1);
    h0 = params(2);
    H = params(3);
    density = rho0 * exp(-(altitude - h0) / H) * 1e9; 

    acc = -0.5 * Cd * AreaMassRatio * density * v_rel * v_rel_vec;

end

function acc = getSRP(r_vec, sun_pos_vec, Cr, AreaMassRatio)

    sat_sun_pos_vec = sun_pos_vec - r_vec;
    sat_sun_pos = norm(sat_sun_pos_vec);

    SF = 1367; %kg/s^3
    c_light = 299792.458; %km/s
    pressureSRP = SF / c_light; 
    gamma = 1;

    % Shadow Model
%     a = 695700 / norm(sun_pos_vec - r_vec); % radius of sun = 695700 km
%     b = myasin(6378.1363 / norm(r_vec)); % radius of earth = 6378.1363 km
%     c = myacos( dot(-r_vec, sun_pos_vec - r_vec) / (norm(r_vec) * norm(sun_pos_vec - r_vec)) );
%     x = (c^2 + a^2 - b^2) / (2 * c);
%     y = sqrt(a^2 - x^2);
%     area = b^2 * myacos((c - x) / b) + a^2 * myacos(x / a) - c * y; % eclipse area
% 
%     if c < abs(a - b)
% 
%         gamma = 0;
% 
%     elseif (a + b) <= c
% 
%         gamma = 1;
% 
%     else
% 
%         gamma = 1 - (area / (pi * a^2));
%         
%     end

    acc = -pressureSRP * gamma * Cr * AreaMassRatio * sat_sun_pos_vec / sat_sun_pos;

end


function pos_GCRF = Sun(JD_UT1)

    T_UT1 = (JD_UT1 - 2451545.0) / 36525;
    
    lambda_mean = 280.460 + 36000.771 * T_UT1;
    
    %Assuming T_TDB = T_UT1
    T_TDB = T_UT1;
    
    mean_anom = 357.52772333 + 35999.05034 * T_TDB;
    
    lambda_ecl = lambda_mean + 1.914666471 * sind(mean_anom) + 0.019994643 * sind(2 * mean_anom);
    
    r = 1.000140612 - 0.016708617 * cosd(mean_anom) - 0.000139589 * cosd(2 * mean_anom);
    
    oblq = 23.439291 - 0.0130042 * T_TDB;
    
    %TOD frame
    pos = [r * cosd(lambda_ecl); r * cosd(oblq) * sind(lambda_ecl); r * sind(oblq) * sind(lambda_ecl)]; % AU
    
    %Assuming TOD = MOD
    %Assuming T_TT = T_TDB = T_UT1
    %1AU = 149597870.7 km
    pos_GCRF = MODtoGCRF(pos * 149597870.7, T_TDB);

end


function Q = MODtoGCRF(r, T_TT)

    zeta = 2306.2181 * T_TT + 0.30188 * T_TT^2 + 0.017998 * T_TT^3;
    theta = 2004.3109 * T_TT - 0.42665 * T_TT^2 - 0.041833 * T_TT^3;
    z = 2306.2181 * T_TT + 1.09468 * T_TT^2 + 0.018203 * T_TT^3;
    
    arcsec_to_deg = @(x) x ./ 3600;
    
    zeta = arcsec_to_deg(zeta);
    theta = arcsec_to_deg(theta);
    z = arcsec_to_deg(z);
    
    M = r3Rotation(zeta) * r2Rotation(-theta) * r3Rotation(z);
    
    Q = M * r;

end


function pos_GCRF = Moon(JD_TDB)

    %Assuming T_TDB = T_UT1

    T_TDB = (JD_TDB - 2451545.0) / 36525;

    lambda_ecl = 218.32 + 481267.8813 * T_TDB + 6.29 * sind(134.9 + 477198.85 * T_TDB) - 1.27 * sind(259.2 - 413335.38 * T_TDB) + 0.66 * sind(235.7 + 890534.23 * T_TDB) + 0.21 * sind(269.9 + 954397.7 * T_TDB) - 0.19 * sind(357.5 + 35999.05 * T_TDB) - 0.11 * sind(186.6 + 966404.05 * T_TDB);

    phi_ecl = 5.13 * sind(93.3 + 483202.03 * T_TDB) + 0.28 * sind(228.2 + 960400.87 * T_TDB) - 0.28 * sind(318.3 + 6003.18 * T_TDB) - 0.17 * sind(217.6 - 407332.20 * T_TDB);

    p = 0.9508 + 0.0518 * cosd(134.9 + 477198.85 * T_TDB) + 0.0095 * cosd(259.2 - 413335.38 * T_TDB) + 0.0078 * cosd(235.7 + 890534.23 * T_TDB) + 0.0028 * cosd(269.9 + 954397.7 * T_TDB);

    oblq = 23.439291 - 0.0130042 * T_TDB - 1.64e-7 * T_TDB^2 + 5.04e-7 * T_TDB^3;

    %Radius of Earth = 6378.1363 km
    r = 6378.1363 / sind(p);

    %J2000.0 Frame
    pos = [r * cosd(phi_ecl) * cosd(lambda_ecl); r * cosd(oblq) * cosd(phi_ecl) * sind(lambda_ecl) - r * sind(oblq) * sind(phi_ecl); r * sind(oblq) * sind(lambda_ecl) + r * cosd(oblq) * sind(phi_ecl)];

    pos_GCRF = J2000toGCRF(pos);

end


function Q = J2000toGCRF(r)

    delta_alpha_0 = 0.0146;
    xi0 = -0.16617;
    eta0 = -0.0068192;
    
    arcsec_to_deg = @(x) x ./ 3600;
    
    delta_alpha_0 = arcsec_to_deg(delta_alpha_0);
    xi0 = arcsec_to_deg(xi0);
    eta0 = arcsec_to_deg(eta0);
    
    M = r3Rotation(-delta_alpha_0) * r2Rotation(-xi0) * r1Rotation(eta0);
    
    Q = M * r;

end


