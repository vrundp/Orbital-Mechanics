function R2 = r2Rotation(x)

R2 = [cosd(x), 0., -sind(x); 0., 1., 0.; sind(x), 0., cosd(x)];

end