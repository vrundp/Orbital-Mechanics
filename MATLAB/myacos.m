function angle = myacos(arg)

    if (abs(arg - 1.0) < 1e-12) && (arg < 1.0)

        arg = 1.0;

    elseif (abs(arg + 1.0) < 1e-12) && (arg < -1.0)

        arg = -1.0;
        
    end
    
    angle = acos(arg);

end