
function COEs = RVtoCOEs(r_vec, v_vec, MU)

    r = norm(r_vec);
    v = norm(v_vec);
    
    energy = (v^2 / 2) - MU / r;
    
    h_vec = cross(r_vec, v_vec);
    h = norm(h_vec);
    
    e_vec = ((v^2 - MU / r) * r_vec - dot(r_vec, v_vec) * v_vec) / MU;
    e = norm(e_vec);
    
    a = -MU / (2 * energy);
    
    i = myacosd(h_vec(3) / h);
    
    if (e < 1e-12) && (i ~= 0)
    
        argP = NaN;
        true_anom = NaN;
        w = NaN;
        l = NaN;
    
        n_vec = cross([0; 0; 1], h_vec);
        n = norm(n_vec);
        RAAN = myacosd(n_vec(1) / n);
        if n_vec(2) < 0
            RAAN = 360 - RAAN;
        end
    
        n_hat = n_vec ./ n;
        u = myacosd(dot(n_hat, r_vec) / r);
        if r_vec(3) < 0
            u = 360 - u;
        end
    
        COEs = [a, e, i, RAAN, argP, true_anom, u, w, l];
    
    
    elseif (i < 1e-12) && (e ~= 0)
    
        RAAN = NaN;
        argP = NaN;
        u = NaN;
        l = NaN;
    
        true_anom = myacosd(dot(r_vec, e_vec) / (r * e));
        if dot(r_vec, v_vec) < 0
            true_anom = 360 - true_anom;
        end
    
        w = myacosd(dot(e_vec, [1; 0; 0]) / e);
        if e_vec(2) < 0
            w = 360 - w;
        end
    
        COEs = [a, e, i, RAAN, argP, true_anom, u, w, l];
    
    elseif (e < 1e-12) && (i < 1e-12)
    
        true_anom = NaN;
        RAAN = NaN;
        argP = NaN; 
        u = NaN;
        w = NaN;
    
        l = myacosd(dot(r_vec, [1; 0; 0]) / r);
        if r_vec(2) < 0
            l = 360 - l;
        end
    
        COEs = [a, e, i, RAAN, argP, true_anom, u, w, l];
      
    else
    
        u = NaN;
        w = NaN;
        l = NaN;
    
        n_vec = cross([0; 0; 1], h_vec);
        n = norm(n_vec);
        RAAN = myacosd(n_vec(1) / n);
        if n_vec(2) < 0
            RAAN = 360 - RAAN;
        end
    
        argP = myacosd((dot(n_vec, e_vec)) / (n * e));
        if e_vec(3) < 0
            argP = 360 - argP;
        end
    
        true_anom = myacosd(dot(r_vec, e_vec) / (r * e));
        if dot(r_vec, v_vec) < 0
            true_anom = 360 - true_anom;
        end
    
        COEs = [a, e, i, RAAN, argP, true_anom, u, w, l];
    
    end

end
