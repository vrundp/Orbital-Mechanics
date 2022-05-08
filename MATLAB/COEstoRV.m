
function RV = COEstoRV(a, e, i, RAAN, argP, true_anom, u, w, l, MU)

    if (e < 1e-12) && (i < 1e-12)

        argP = 0.0;
        RAAN = 0.0;
        true_anom = l;

    elseif (e < 1e-12) && (i ~= 0)

        argP = 0.0;
        true_anom = u;

    elseif (i < 1e-12) && (e ~= 0)

        RAAN = 0.0;
        argP = w;

    end

    p = a * (1 - e^2);

    r_vec_pqw = [p * cosd(true_anom) / (1 + e * cosd(true_anom)); p * sind(true_anom) / (1 + e * cosd(true_anom)); 0];
    v_vec_pqw = [-sqrt(MU / p) * sind(true_anom); sqrt(MU / p) * (e + cosd(true_anom)); 0];

    Q = r3Rotation(-RAAN) * r1Rotation(-i) * r3Rotation(-argP);

    r_vec_ijk = Q * r_vec_pqw;
    v_vec_ijk = Q * v_vec_pqw;

    RV = [r_vec_ijk; v_vec_ijk];

end