centerPos = 103; %103 for not 32 km/s, 113 for 32 km/s
% 2.301533079552962e+02 nm for 8 km/s
% 4.986655005698084e+02 nm for 16 km/s
% 1.264222395810782e+03 nm for 32 km/s
wavelength = 4.986655005698084e+02; % nm
x_range = wavelength/x_spacing;
half_range = round(x_range/2);
left_x = centerPos-half_range;
right_x = centerPos+half_range;

density = 5.323/1000*(100)^3; % kg/m^3

Tot_Kin_Energ = 0;

Vol_Mult = (x_incr^2+y_incr^2)^(1/2)*(x_spacing/sin(atan(abs(y_incr)/abs(x_incr))))*1; % nm^3
Vol_Mult = Vol_Mult/(1e9)^2/1e3; % m^2*um

for curX = left_x:right_x
    
    curVel_Vec = vel_matr{curX}*1000; % pm/ps or m/s
    for dist_int = 1:size(curVel_Vec,1)
        Tot_Kin_Energ = Tot_Kin_Energ + density.*Vol_Mult.*(curVel_Vec(dist_int,1).^2 + curVel_Vec(dist_int,2).^2)./2;
    end
    
end