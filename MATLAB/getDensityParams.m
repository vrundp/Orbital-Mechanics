
function params = getDensityParams( altitude )


% Density model parameters for the 1976 Standard Atmosphere exponential model
% 
%  rho0 - Nominal density in kg/m^3
%  h0   - Base altitude in km
%  H    - Scale height in km
%
% Assumptions/References:
%
%  Vallado and McClain, Fourth Edition.
%



    if altitude > 1000

        rho0 = 3.019e-15;
        h0 = 1000;
        H = 268;

        params = [rho0, h0, H];

    elseif altitude > 900

        rho0 = 5.245e-15;
        h0 = 900;
        H = 181.05;

        params = [rho0, h0, H];

    elseif altitude > 800

        rho0 = 1.170e-14;
        h0 = 800;
        H = 124.64;

        params = [rho0, h0, H];

    elseif altitude > 700

        rho0 = 3.614e-14;
        h0 = 700;
        H = 88.667;

        params = [rho0, h0, H];

    elseif altitude > 600

        rho0 = 1.454e-13;
        h0 = 600;
        H = 71.835;

        params = [rho0, h0, H];

    elseif altitude > 500

        rho0 = 6.967e-13;
        h0 = 500;
        H = 63.822;

        params = [rho0, h0, H];

    elseif altitude > 450

        rho0 = 1.585e-12;
        h0 = 450;
        H = 60.828;

        params = [rho0, h0, H];

    elseif altitude > 400

        rho0 = 3.725e-12;
        h0 = 400;
        H = 58.515;

        params = [rho0, h0, H];

    elseif altitude > 350

        rho0 = 9.518e-12;
        h0 = 350;
        H = 53.298;

        params = [rho0, h0, H];

    elseif altitude > 300

        rho0 = 2.418e-11;
        h0 = 300;
        H = 53.628;

        params = [rho0, h0, H];

    elseif altitude > 250

        rho0 = 7.248e-11;
        h0 = 250;
        H = 45.546;

        params = [rho0, h0, H];

    elseif altitude > 200

        rho0 = 2.789e-10;
        h0 = 200;
        H = 37.105;

        params = [rho0, h0, H];

    elseif altitude > 180

        rho0 = 5.464e-10;
        h0 = 180;
        H = 29.740;

        params = [rho0, h0, H];

    elseif altitude > 150

        rho0 = 2.070e-9;
        h0 = 150;
        H = 22.523;

        params = [rho0, h0, H];

    elseif altitude > 140

        rho0 = 3.845e-9;
        h0 = 140;
        H = 16.149;

        params = [rho0, h0, H];

    elseif altitude > 130

        rho0 = 8.484e-9;
        h0 = 130;
        H = 12.636;

        params = [rho0, h0, H];

    elseif altitude > 120
        rho0 = 2.438e-8;
        h0 = 120;
        H = 9.473;

        params = [rho0, h0, H];

    elseif altitude > 110

        rho0 = 9.661e-8;
        h0 = 110;
        H = 7.263;

        params = [rho0, h0, H];

    elseif altitude > 100

        rho0 = 5.297e-7;
        h0 = 100;
        H = 5.877;

        params = [rho0, h0, H];

    elseif altitude > 90

        rho0 = 3.396e-6;
        h0 = 90;
        H = 5.382;

        params = [rho0, h0, H];

    elseif altitude > 80

        rho0 = 1.905e-5;
        h0 = 80;
        H = 5.799;

        params = [rho0, h0, H];

    elseif altitude > 70

        rho0 = 8.770e-5;
        h0 = 70;
        H = 6.549;

        params = [rho0, h0, H];

    elseif altitude > 60

        rho0 = 3.206e-4;
        h0 = 60;
        H = 7.714;

        params = [rho0, h0, H];

    elseif altitude > 50

        rho0 = 1.057e-3;
        h0 = 50;
        H = 8.382;

        params = [rho0, h0, H];

    elseif altitude > 40

        rho0 = 3.972e-3;
        h0 = 40;
        H = 7.554;

        params = [rho0, h0, H];

    elseif altitude > 30

        rho0 = 1.774e-2;
        h0 = 30;
        H = 6.682;

        params = [rho0, h0, H];

    elseif altitude > 25

        rho0 = 3.899e-2;
        h0 = 25;
        H = 6.349;

        params = [rho0, h0, H];

    else

        rho0 = 1.225;
        h0 = 0;
        H = 7.249;

        params = [rho0, h0, H];

    end

end
