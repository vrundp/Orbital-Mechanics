function R1 = r1Rotation(x)

R1 = [1., 0., 0.; 0., cosd(x), sind(x); 0., -sind(x), cosd(x)];

end