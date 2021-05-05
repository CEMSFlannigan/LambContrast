numFrames = 1;
centerPoint = 113;

% 103 for 8V and 16V
% 113 for 32V due to large wavelength

% 21 ps for 16 km/s
% 7 ps for 8 km/s
% 28 ps for 32 km/s

min_bf_diff_loc_ind = 1;

wavelength = 1.264222395810782e+03; % nm
% 2.301533079552962e+02 nm for 8 km/s
% 4.986655005698084e+02 nm for 16 km/s
% 1.264222395810782e+03 nm for 32 km/s
x_range = wavelength/x_spacing;
half_range = round(x_range/2);
left_x = centerPoint-half_range;
right_x = centerPoint+half_range;

C11 = 1.292e11; % Pa
C12 = 0.479e11; % Pa
C44 = 0.670e11; % Pa

Tot_Dilatational = 0;
Tot_Shear = 0;

Vol_Mult = (x_incr^2+y_incr^2)^(1/2)*(x_spacing/sin(atan(abs(y_incr)/abs(x_incr))))*1; % nm^3
Vol_Mult = Vol_Mult/(1e9)^2/1e3; % m^2*um

for curX = left_x:right_x
     
    StrainMatrTim = strain_matr{curX,min_bf_diff_loc_ind};
    for dist_int = 1:size(StrainMatrTim,1)
        Tot_Dilatational = Tot_Dilatational + Vol_Mult/2*(C11*StrainMatrTim(dist_int,1)^2 + C11*StrainMatrTim(dist_int,2)^2 + abs(2*C12*StrainMatrTim(dist_int,1)*StrainMatrTim(dist_int,2)));
        Tot_Shear = Tot_Shear + Vol_Mult/2*(C44*StrainMatrTim(dist_int,3)^2);
    end
    
end

Tot_Dilatational
Tot_Shear