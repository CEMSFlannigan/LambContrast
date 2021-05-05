%% Setting physical constants
lambda = 2.5079e-3*1e-9; % m, 200 keV
vel = 2.0844e8; % m/s, 200 keV
theta_0 = 0.117*32^(1/3)/(200^(1/2)); % rad
beta = 0.0038; % rad
m_0 = 9.11e-31; % kg
rel_mass = 1.2674e-30;% kg
a_0 = 5.29e-11; % m
Z = 32;
E_0 = 200*1000*1.602e-19; % J, 200 keV
light = 2.9979e8; % m/s
rho = 5.323/1000*(100)^3; % kg/m^3
N_0 = 6.022e23; % atoms/mol
A = 72.64/1000; % kg/mol
gam = 1/sqrt(1-vel^2/light^2); % nondim
sigma = (Z*lambda*(a_0/(Z^0.33))*(1 + E_0/m_0/light^2))^2/(pi*a_0^2*(1+(beta/theta_0)^2));

e_charge = 1.602e-19; % C
permitiv = 8.8541878e-12; % C/V/m
planck = 6.626e-34; % Js

%% Declare variables and constants relating to electron beam travel direction
k_mag = 398; % nm^-1
rel_acc = ((k_mag*1e9)*6.626e-34)^2/2/rel_mass/e_charge; % V
g_mags = 3.06; % nm^-1 (1-1-1)
g_mags_2 = 3.06; % nm^-1 (1-11)
alpha_tilt = -13.02;
gam_tilt = 0;
beta_tilt = 1.55;
beta_add = 0; 
alpha_add = 0.500; 
k_unit = [0, -1, 0];
k_unit = k_unit/sqrt(dot(k_unit,k_unit));
k_alph_rot = [1,0,0;0,cosd(alpha_tilt),-sind(alpha_tilt);0,sind(alpha_tilt),cosd(alpha_tilt)];
k_gam_rot = [cosd(gam_tilt),0,sind(gam_tilt);0,1,0;-sind(gam_tilt),0,cosd(gam_tilt)];
k_beta_rot = [cosd(beta_tilt),-sind(beta_tilt),0;sind(beta_tilt),cosd(beta_tilt),0;0,0,1];
k_alph_add = [1,0,0;0,cosd(alpha_add),-sind(alpha_add);0,sind(alpha_add),cosd(alpha_add)];
k_beta_add = [cosd(beta_add),-sind(beta_add),0;sind(beta_add),cosd(beta_add),0;0,0,1];
k_unit = k_beta_rot*k_gam_rot*k_alph_rot*k_beta_add*k_alph_add*k_unit';

% Formerly used equation for Rutherfordian scattering (not used in this code)
frac_scatter = @(t) exp(-(N_0./A).*sigma.*rho.*t);

% Atomic scattering factor constants with parameterization via Doyle-Turner
SFa1 = 0.2135;
SFa2 = 0.9761;
SFa3 = 1.6555;
SFa4 = 2.8938;
SFa5 = 1.6356;
SFb1 = 0.0989;
SFb2 = 0.9845;
SFb3 = 4.5527;
SFb4 = 21.5563;
SFb5 = 70.3903;
Ge_ASF = @(thet) (SFa1*exp(-SFb1*(sind(thet)/(lambda*1e10)).^2) + SFa2*exp(-SFb2*(sind(thet)/(lambda*1e10)).^2) + SFa3*exp(-SFb3*(sind(thet)/(lambda*1e10)).^2) + SFa4*exp(-SFb4*(sind(thet)/(lambda*1e10)).^2) + SFa5*exp(-SFb5*(sind(thet)/(lambda*1e10)).^2)); % A

% Absorption form factor constants with parameterization via Doyle-Turner
% For T = 293 K, 200keV
AFFa1 = -0.0008;
AFFa2 = -0.0281;
AFFa3 = 0.1559;
AFFa4 = 0.0873;
AFFa5 = 0.0110;
AFFb1 = 0.0002;
AFFb2 = 0.1071;
AFFb3 = 0.5417;
AFFb4 = 1.6461;
AFFb5 = 8.4617;
DW = 0.6000;
Beta100 = (1-(1+1.9579341e-3*100)^(-2))^(-1);
Beta200 = (1-(1+1.9579341e-3*200)^(-2))^(-1);
Ge_AFF = @(thet) (AFFa1*exp(-AFFb1*(sind(thet)/(lambda*1e10)).^2) + AFFa2*exp(-AFFb2*(sind(thet)/(lambda*1e10)).^2) + AFFa3*exp(-AFFb3*(sind(thet)/(lambda*1e10)).^2) + AFFa4*exp(-AFFb4*(sind(thet)/(lambda*1e10)).^2) + AFFa5*exp(-AFFb5*(sind(thet)/(lambda*1e10)).^2))*Beta100/Beta200./exp(-DW/2*(sind(thet)/(lambda*1e10)).^2);

% Crystal lattice contstants and probing dimensions used later in the code for probing strain or x-positions
latpar = 0.5658; % nm
strainpar = sqrt(dot(1/4*[1,1,1],1/4*[1,1,1]))*latpar; % nm
V_c = latpar^3;
ExtDist = @(thet,h,k,l) pi*(V_c./1e9./1e9./1e9)*cosd(thet)/lambda/(Ge_ASF(thet)/1e10)/real(Ge_SF(h,k,l));
x_spacing = 10*latpar; % nm

k_mag = sqrt(2*rel_mass*e_charge*(rel_acc + 47.87801/(V_c*10^3)*(Ge_ASF(0))))/planck*1e-9; % nm^-1
k_actual = k_unit*k_mag;

% Unit cell positions
r1 = [0, 0, 0];
r2 = [1/4, 1/4, 1/4];
r3 = [0, 1/2, 1/2];
r4 = [1/2, 0, 1/2];
r5 = [1/2, 1/2, 0];
r6 = [1/4, 3/4, 3/4];
r7 = [3/4, 1/4, 3/4];
r8 = [3/4, 3/4, 1/4];

% Diamond cubic structure factor
Ge_SF = @(ha,ka,la,r1a,r2a,r3a,r4a,r5a,r6a,r7a,r8a) exp(-2*pi*sqrt(-1)*(ha*r1a(1) + ka*r1a(2) + la*r1a(3))) + exp(-2*pi*sqrt(-1)*(ha*r2a(1) + ka*r2a(2) + la*r2a(3))) + exp(-2*pi*sqrt(-1)*(ha*r3a(1) + ka*r3a(2) + la*r3a(3))) + exp(-2*pi*sqrt(-1)*(ha*r4a(1) + ka*r4a(2) + la*r4a(3))) + exp(-2*pi*sqrt(-1)*(ha*r5a(1) + ka*r5a(2) + la*r5a(3))) + exp(-2*pi*sqrt(-1)*(ha*r6a(1) + ka*r6a(2) + la*r6a(3))) + exp(-2*pi*sqrt(-1)*(ha*r7a(1) + ka*r7a(2) + la*r7a(3))) + exp(-2*pi*sqrt(-1)*(ha*r8a(1) + ka*r8a(2) + la*r8a(3)));

% Real Fourier coefficients of the crystal lattice (elastic Bragg-Bragg scattering)
Ge_RPOT = @(thet,gx,gy,gz) 47.87801/(V_c*10^3)*(Ge_ASF(thet))*2*rel_mass*e_charge/(planck)^2./(sqrt((k_actual(1)*ones(size(gx))+gx).^2 + (k_actual(2)*ones(size(gy))+gy).^2 + (k_actual(3)*ones(size(gz))+gz).^2)*1e9);
% Imaginary Fourier coefficients of the crystal lattice (diffuse scattering (anomalous), absorption)
Ge_IPOT = @(thet,gx,gy,gz) 47.87801/(V_c*10^3)*(Ge_AFF(thet))*2*rel_mass*e_charge/(planck)^2./(sqrt((k_actual(1)*ones(size(gx))+gx).^2 + (k_actual(2)*ones(size(gy))+gy).^2 + (k_actual(3)*ones(size(gz))+gz).^2)*1e9);

% Distance by which to increment along the path of the electron, assuming column approximation
thick_incr = 0.1; % nm
xy_dist_incr = sqrt(thick_incr^2*(1-k_unit(3)^2));
%xy_dist_incr = 0.0001; % nm
%dist_const = xy_dist_incr^2/(k_unit(1)^2+k_unit(2)^2);
%thick_incr = sqrt(xy_dist_incr^2 + dist_const^2*k_unit(3)^2); % nm

% Setting the initial shape of the wedge
yinit = 25; %[25 65]; % nm
xvals = 0:x_spacing:1000; % nm

% Setting the timepoints
timevals = 7; % ps

% Setting the shape of the wedge
angle = 2/2; % degrees!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
slope = sind(angle); % nm/nm

% Setting the shape of the undisturbed surfaces of the wedge
y_top = @(x) slope*x+yinit;
y_bot = @(x) -slope*x-yinit;

% Setting distortion parameters
LambAmp = 3;%0.005;%3;%6;%4e2;%3.75e1;%3.75; % nm
k = 0.0273; %.0049693; %0.010 %0.0270; % 1/nm
omega = (35.7400e9/1e12)*2*pi; % 1/ps
VL = 5.350; % nm/ps
VT = 3.570; % nm/ps
V = 8; % nm/ps
h = 36; % nm
alpha = 0; % Set to pi/2 for asymmetric mode

p = sqrt(omega^2*(1/VL^2 - 1/V^2));
q = sqrt(omega^2*(1/VT^2 - 1/V^2));

syms x_dist(y) y_dist(y) time_space(x,t) combined_x(x,c,t,d) combined_y(x,y,t)

% Distortion in x and y
x_dist(y) = q*LambAmp*(cos(q*y+alpha) - 2*k^2/(k^2-q^2)*cos(q*h+alpha)/cos(p*h+alpha)*cos(p*y+alpha));
y_dist(y) = k*LambAmp*(sin(q*y+alpha) + 2*p*q/(k^2-q^2)*cos(q*h+alpha)/cos(p*h+alpha)*sin(p*y+alpha));
time_space(x,t) = cos(k*x-omega*t);
combined_x(x,y,t) = x_dist(y).*time_space(x,t);
combined_y(x,y,t) = y_dist(y).*time_space(x,t);

strain_matr = cell(length(xvals), length(timevals));

% setting up the Bragg spots
simspots = 400;
spot_angles = zeros(1,simspots);
ext_dist_arr = zeros(1,simspots);
SF_arr = zeros(1,simspots);
gvec_list = cell(1,simspots);
hvec_list = cell(1,simspots);
Vg_mag_list = zeros(1,simspots);

% Bragg spot basis
g_basis_1 = [1,-1,-1];
g_basis_2 = [1,-1,1];
or_z_basis = [1, g_basis_1(2), 0];
or_true_basis = [1, 0, 0];
basis_angle = acosd(dot(g_basis_1,g_basis_2)/sqrt(sum(g_basis_1.^2)*sum(g_basis_2.^2)));
orientation_z_angle = acosd(dot(g_basis_1,or_z_basis)/sqrt(sum(g_basis_1.^2)*sum(or_z_basis.^2)));
orientation_x_angle = acosd(dot(or_z_basis,or_true_basis)/sqrt(sum(or_z_basis.^2)*sum(or_true_basis.^2)));

% Orienting the first Bragg spot
g_alpha = 0;
g_beta = 1;%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
g_unit = [1, 0, 0]; % 0.0353 straight perpendicular ---> ZOLZ Weiss Zone Rule
g_unit = g_unit/sqrt(dot(g_unit,g_unit));
g_alph_rot = [cosd(g_alpha),0,sind(g_alpha);0,1,0;-sind(g_alpha),0,cosd(g_alpha)];
g_beta_rot = [cosd(g_beta),-sind(g_beta),0;sind(g_beta),cosd(g_beta),0;0,0,1];
g_unit = (g_beta_rot*g_alph_rot*g_unit')';

or_z_rot = [cosd(orientation_z_angle),-sind(orientation_z_angle),0;sind(orientation_z_angle),cosd(orientation_z_angle),0;0,0,1];
or_x_rot = [cosd(-orientation_x_angle),0,sind(-orientation_x_angle);0,1,0;-sind(-orientation_x_angle),0,cosd(-orientation_x_angle)];
true_x_dir = (or_z_rot*or_x_rot*g_unit')';
true_y_dir = ([cosd(90), -sind(90), 0; sind(90), cosd(90), 0; 0, 0, 1]*true_x_dir')';
true_z_dir = ([cosd(-90), 0, sind(-90); 0, 1, 0; -sind(-90), 0, cosd(-90)]*true_x_dir')';

% Finding the second Bragg spot basis
g_unit_2 = ([cosd(basis_angle),0,sind(basis_angle);0,1,0;-sind(basis_angle),0,cosd(basis_angle)]*g_unit')';
g_unit_2 = g_unit_2/sqrt(dot(g_unit_2,g_unit_2));

gvec_list{1} = g_unit.*g_mags;
gvec_list{2} = g_unit_2.*g_mags_2;
gvec_list{3} = -g_unit.*g_mags;
gvec_list{4} = -g_unit_2.*g_mags_2;

hvec_list{1} = g_basis_1;
hvec_list{2} = g_basis_2;
hvec_list{3} = -g_basis_1;
hvec_list{4} = -g_basis_2;

% Generate Bragg spots
[gvec_list,hvec_list] = RecipSpotGenerator(4,gvec_list,simspots,hvec_list);

% Function to calculate the excitation error
s_g = @(g)  (k_mag^2 - dot(k_actual' + g,k_actual'+g))/(2*sqrt(dot(k_actual'+g,k_actual'+g)));

% Bragg angles
angles = @(x,y,z) asind(lambda/2*sqrt(x.*x+y.*y+z.*z)*1e9);

% Determine which surface is the first surface of contact based on electron path travel
if k_unit(2) > 0
    start_y = @(x) y_bot(x);
    end_y = @(x) y_top(x);
    top_bot = 1;
elseif k_unit(2) < 0
    start_y = @(x) y_top(x);
    end_y = @(x) y_bot(x);
    top_bot = 0;
else
    throw('Stop using weird geometries. This will not propagate through the specimen.');
end

% Calculating slopes of the surfaces and associated inverse slopes
top_slope = slope;
bot_slope = -slope;
top_inv_slope = -1/top_slope;
bot_inv_slope = -1/bot_slope;
top_inv_vec = [1, top_inv_slope, 0];
top_inv_vec = top_inv_vec/sqrt(dot(top_inv_vec,top_inv_vec));
bot_inv_vec = [1, bot_inv_slope, 0];
bot_inv_vec = bot_inv_vec/sqrt(dot(bot_inv_vec,bot_inv_vec));

% Setting up intensities of each Bragg spot
summed_intensities = gpuArray(zeros(1,simspots+1));
summed_intensities(1) = 1;

% Setting up contrast
bf_contrast = gpuArray(zeros(length(xvals),length(timevals)));

% Setting up the travel distance of each increment
k_slope = k_unit(2)/k_unit(1);
x_incr = sign(k_unit(1))*xy_dist_incr/sqrt(1+k_slope^2);
y_incr = sign(k_unit(2))*xy_dist_incr/sqrt(1+(1/k_slope)^2);
dist_travel = sqrt(abs(x_incr)^2 + abs(y_incr)^2 + abs(k_unit(3)/k_unit(2)*y_incr)^2);

tic;
for k = 1:length(timevals)
    parfor i = 1:length(xvals)
        
        summed_intensities = zeros(1,simspots+1);
        summed_intensities(1) = 1;
        
        curx = xvals(i);
        cur_start_y = start_y(curx);
        cury = cur_start_y;
        
        tot_thick = 0;
        num_incr = 0;
        strain_vec = [];
        
        prevx_lamb = curx + combined_x(curx, cury, timevals(k));
        prevy_lamb = cury + combined_y(curx, cury, timevals(k));
        while abs(cury) <= abs(end_y(curx))
            
            s_g_vec = zeros(simspots+1,3);
            s_g_mag = zeros(simspots+1,1);
            mod_gvec_list = zeros(simspots+1,3);
            mod_hvec_list = zeros(simspots+1,3);
            g_angles = zeros(simspots+1,simspots+1);
            gx_dif_matr = zeros(simspots+1,simspots+1);
            gy_dif_matr = zeros(simspots+1,simspots+1);
            gz_dif_matr = zeros(simspots+1,simspots+1);
            hx_dif_matr = zeros(simspots+1,simspots+1);
            hy_dif_matr = zeros(simspots+1,simspots+1);
            hz_dif_matr = zeros(simspots+1,simspots+1);
            pot_phase = zeros(simspots+1, simspots+1);
            pot_mag = zeros(simspots+1, simspots+1);
            pot_mag_imag = zeros(simspots+1, simspots+1);
            pot_abs = zeros(simspots+1, simspots+1);
            pot_phase_fact_r = zeros(simspots+1, simspots+1);
            pot_phase_fact_i = zeros(simspots+1, simspots+1);
            U_matr = zeros(simspots+1, simspots+1);
            
            curx_lamb = curx + combined_x(curx, cury, timevals(k));
            cury_lamb = cury + combined_y(curx, cury, timevals(k));
            
            %% CURRENT ISSUE: NEED TO FIND s_g
            
            %% CALCULATE THE STRAIN IN THIS REGION VIA 6 POINTS
            % SEPARATION DISTANCE WILL BE ARBITRARY (NN BOND DISTANCE)
            
            x1 = curx - strainpar/2;
            x2 = curx + strainpar/2;
            y1 = cury - strainpar/2;
            y2 = cury + strainpar/2;
            %z1
            %z2
            
            ux1 = double(combined_x(x1, cury, timevals(k)));
            ux2 = double(combined_x(x2, cury, timevals(k)));
            uy1 = double(combined_y(curx, y1, timevals(k)));
            uy2 = double(combined_y(curx, y2, timevals(k)));
            %uz1
            %uz2
            
            strain_x = (ux2-ux1)/(x2-x1).*[1, 0, 0];
            strain_x_mag = (ux2-ux1)/(x2-x1);
            strain_y = (uy2-uy1)/(y2-y1).*[0, 1, 0];
            strain_y_mag = (uy2-uy1)/(y2-y1);
            strain_z = 0.*[0, 0, 1];
            strain_z_mag = 0;
            
            strain_true_x = dot(strain_x,true_x_dir)*true_x_dir + dot(strain_y,true_x_dir)*true_x_dir + dot(strain_z,true_x_dir)*true_x_dir;
            strain_true_x_mag = dot(strain_x,true_x_dir) + dot(strain_y,true_x_dir) + dot(strain_z,true_x_dir);
            strain_true_y = dot(strain_x,true_y_dir)*true_y_dir + dot(strain_y,true_y_dir)*true_y_dir + dot(strain_z,true_y_dir)*true_y_dir;
            strain_true_y_mag = dot(strain_x,true_y_dir) + dot(strain_y,true_y_dir) + dot(strain_z,true_y_dir);
            strain_true_z = dot(strain_x,true_z_dir)*true_z_dir + dot(strain_y,true_z_dir)*true_z_dir + dot(strain_z,true_z_dir)*true_z_dir;
            strain_true_z_mag = dot(strain_x,true_z_dir) + dot(strain_y,true_z_dir) + dot(strain_z,true_z_dir);
            
            strain_vec = [strain_vec; strain_true_x_mag, strain_true_y_mag, strain_true_z_mag];
            
            strain_tot_vec = strain_x+strain_y+strain_z;
            strain_tot_mag = sqrt(dot(strain_tot_vec,strain_tot_vec));
            
            r1n = [r1(1)*(1+strain_true_x_mag),r1(2)*(1+strain_true_y_mag),r1(3)*(1+strain_true_z_mag)];
            r2n = [r2(1)*(1+strain_true_x_mag),r2(2)*(1+strain_true_y_mag),r2(3)*(1+strain_true_z_mag)];
            r3n = [r3(1)*(1+strain_true_x_mag),r3(2)*(1+strain_true_y_mag),r3(3)*(1+strain_true_z_mag)];
            r4n = [r4(1)*(1+strain_true_x_mag),r4(2)*(1+strain_true_y_mag),r4(3)*(1+strain_true_z_mag)];
            r5n = [r5(1)*(1+strain_true_x_mag),r5(2)*(1+strain_true_y_mag),r5(3)*(1+strain_true_z_mag)];
            r6n = [r6(1)*(1+strain_true_x_mag),r6(2)*(1+strain_true_y_mag),r6(3)*(1+strain_true_z_mag)];
            r7n = [r7(1)*(1+strain_true_x_mag),r7(2)*(1+strain_true_y_mag),r7(3)*(1+strain_true_z_mag)];
            r8n = [r8(1)*(1+strain_true_x_mag),r8(2)*(1+strain_true_y_mag),r8(3)*(1+strain_true_z_mag)];
            
            %% THEN NEED TO FIND s_g via change in g
            
            for hvec_count = 1:simspots
                
                curg_vec = gvec_list{hvec_count};
                curh_vec = hvec_list{hvec_count};

                mod_gvec_list(hvec_count+1,:) = [curg_vec(1)/(1+strain_x_mag), curg_vec(2)/(1+strain_y_mag), curg_vec(3)];%curg_vec+g_distort*curg_vec./curg_mag;
                
                mod_hvec_list(hvec_count+1,:) = [curh_vec(1)/(1+strain_true_x_mag), curh_vec(2)/(1+strain_true_y_mag), curh_vec(3)/(1+strain_true_z_mag)];
                
                curg_vec_str = mod_gvec_list(hvec_count+1,:);
                
                dir_vec = (k_actual' + curg_vec_str)/sqrt(dot(k_actual'+curg_vec_str,k_actual'+curg_vec_str));
                
                s_g_vec(hvec_count+1,:) = dir_vec*s_g(curg_vec_str);
                
                s_g_mag(hvec_count+1,:) = s_g(curg_vec_str)*1e9;
            
            end
            
%             for g_count = 1:simspots+1
%                 for g2_count = g_count:simspots+1
%                     
%                     g_dif_vec = mod_gvec_list(g2_count,:) - mod_gvec_list(g_count,:);
%                     h_dif_vec = mod_hvec_list(g2_count,:) - mod_hvec_list(g_count,:);
%                     
%                     g_angle = asind(lambda/2*sqrt(dot(g_dif_vec,g_dif_vec))*1e9);
%                     
%                     if sqrt(dot(h_dif_vec,h_dif_vec)) ~= 0
%                         pot_phase(g_count,g2_count) = exp(-2*pi*sqrt(-1)*dot(h_dif_vec,r1n)) + exp(-2*pi*sqrt(-1)*dot(h_dif_vec,r2n)) + exp(-2*pi*sqrt(-1)*dot(h_dif_vec,r3n)) + exp(-2*pi*sqrt(-1)*dot(h_dif_vec,r4n)) + exp(-2*pi*sqrt(-1)*dot(h_dif_vec,r5n)) + exp(-2*pi*sqrt(-1)*dot(h_dif_vec,r6n)) + exp(-2*pi*sqrt(-1)*dot(h_dif_vec,r7n)) + exp(-2*pi*sqrt(-1)*dot(h_dif_vec,r8n));
%                         pot_phase(g2_count,g_count) = pot_phase(g_count,g2_count);
%                         pot_mag(g_count,g2_count) = 47.87801/(V_c*10^3)*(Ge_ASF(g_angle))*2*rel_mass*e_charge/(planck)^2./(sqrt(dot(k_actual'+g_dif_vec,k_actual'+g_dif_vec))*1e9);
%                         pot_mag(g2_count,g_count) = pot_mag(g_count,g2_count);
%                         
%                         pot_mag_imag(g_count,g2_count) = 47.87801/(V_c*10^3)*(Ge_AFF(g_angle))*2*rel_mass*e_charge/(planck)^2./(sqrt(dot(k_actual'+g_dif_vec,k_actual'+g_dif_vec))*1e9);
%                         pot_mag_imag(g2_count,g_count) = pot_mag_imag(g_count,g2_count);
%                     else
%                         pot_phase(g_count,g2_count) = exp(-2*pi*sqrt(-1)*dot(h_dif_vec,r1n)) + exp(-2*pi*sqrt(-1)*dot(h_dif_vec,r2n)) + exp(-2*pi*sqrt(-1)*dot(h_dif_vec,r3n)) + exp(-2*pi*sqrt(-1)*dot(h_dif_vec,r4n)) + exp(-2*pi*sqrt(-1)*dot(h_dif_vec,r5n)) + exp(-2*pi*sqrt(-1)*dot(h_dif_vec,r6n)) + exp(-2*pi*sqrt(-1)*dot(h_dif_vec,r7n)) + exp(-2*pi*sqrt(-1)*dot(h_dif_vec,r8n));
%                         pot_phase(g2_count,g_count) = pot_phase(g_count,g2_count);
%                     end
%                 end
%             end

            D = bsxfun(@minus, mod_gvec_list(:,1), mod_gvec_list(:,1)');                                  % Subtract To Create Matrix
            gx_dif_matr = D;
            
            D = bsxfun(@minus, mod_gvec_list(:,2), mod_gvec_list(:,2)');                                  % Subtract To Create Matrix
            gy_dif_matr = D;
            
            D = bsxfun(@minus, mod_gvec_list(:,3), mod_gvec_list(:,3)');                                  % Subtract To Create Matrix
            gz_dif_matr = D;
            
            D = bsxfun(@minus, mod_hvec_list(:,1), mod_hvec_list(:,1)');                                  % Subtract To Create Matrix
            hx_dif_matr = D;
            
            D = bsxfun(@minus, mod_hvec_list(:,2), mod_hvec_list(:,2)');                                  % Subtract To Create Matrix
            hy_dif_matr = D;
            
            D = bsxfun(@minus, mod_hvec_list(:,3), mod_hvec_list(:,3)');                                  % Subtract To Create Matrix
            hz_dif_matr = D;
            
            g_angles = gpuArray(angles(gx_dif_matr,gy_dif_matr,gz_dif_matr));
            
            pot_phase = gpuArray(Ge_SF(hx_dif_matr,hy_dif_matr,hz_dif_matr,r1n,r2n,r3n,r4n,r5n,r6n,r7n,r8n));
            
            pot_mag = gpuArray(Ge_RPOT(g_angles,gx_dif_matr,gy_dif_matr,gz_dif_matr));
            
            pot_mag_imag = gpuArray(Ge_IPOT(g_angles,gx_dif_matr,gy_dif_matr,gz_dif_matr));
            
            for strack = 1:simspots+1
                pot_mag(strack,strack) = 0;
                pot_mag_imag(strack,strack) = 0;
            end
            
            U_matr_imag = gpuArray(abs(pot_phase.*pot_mag_imag));
            U_matr_real = gpuArray(abs(pot_phase.*pot_mag));
            U_matr = gpuArray(exp(sqrt(-1)*pot_phase_fact_r).*U_matr_real);% + exp(sqrt(-1)*pot_phase_fact_i)*sqrt(-1).*U_matr_imag; %don't want TDS yet
            
            if k_unit(1) ~= 0
                true_thick_incr = double(sqrt((curx_lamb - prevx_lamb)^2 + (cury_lamb - prevy_lamb)^2 + (k_unit(3)/k_unit(1)*(curx_lamb - prevx_lamb))^2));
            elseif k_unit(2) ~= 0
                true_thick_incr = double(sqrt((curx_lamb - prevx_lamb)^2 + (cury_lamb - prevy_lamb)^2 + (k_unit(3)/k_unit(2)*(cury_lamb - prevy_lamb))^2));
            else
                true_thick_incr = thick_incr;
            end
            
            exc_err_matr = gpuArray(2*pi*diag(s_g_mag));
            interact_matr = gpuArray(pi*U_matr);
            
            DHW_matr = gpuArray(exc_err_matr + interact_matr);
            
            exp_matr = gpuArray(sqrt(-1)*(DHW_matr)*true_thick_incr*1e-9);
            
%             assignin('base','pot_phase_old',pot_phase);
%             assignin('base','g_angles_old',g_angles);
%             assignin('base','s_g_old',s_g_mag);
%             assignin('base','U_matr_old',U_matr);
%             assignin('base','U_matr_imag_old',abs(pot_phase.*pot_mag_imag));
%             assignin('base','pot_phase_fact_r_old',pot_phase_fact_r);
%             assignin('base','pot_phase_fact_i_old',pot_phase_fact_i);
%             assignin('base','U_matr_real_old',abs(pot_phase.*pot_mag));
%             assignin('base','pot_phase_old',pot_phase);
%             assignin('base','DHW_matr_old',DHW_matr);
%             assignin('base','exp_matr_old',exp_matr);
            
            summed_intensities = expm(exp_matr)*summed_intensities';
            summed_intensities = summed_intensities';
            
%             change_summed_intensities_exc_err = pi*sqrt(-1)*(true_thick_incr*1e-9)*2*s_g_mag.*summed_intensities';
%             change_summed_intensities_interact =  pi*sqrt(-1)*(true_thick_incr*1e-9)*U_matr*summed_intensities';
%             
%             summed_intensities = summed_intensities+change_summed_intensities_exc_err' + change_summed_intensities_interact';
            summed_intensities = summed_intensities./sqrt(sum(abs(summed_intensities).^2));
            %abs(summed_intensities(1))
            
            curx = curx + x_incr;
            cury = cury + y_incr;
            
            tot_thick = true_thick_incr+tot_thick;
            
            prevx_lamb = curx_lamb;
            prevy_lamb = cury_lamb;
            
        end
        
        strain_matr{i,k} = strain_vec;
        
        absorb = Ge_SF(0,0,0,r1,r2,r3,r4,r5,r6,r7,r8)*47.87801/(V_c*10^3)*(Ge_AFF(0))*2*rel_mass*e_charge/(planck)^2./(sqrt(dot(k_actual',k_actual'))*1e9);
        summed_intensities = summed_intensities*exp(-2*pi*(tot_thick*1e-9)*absorb);
        
        bf_contrast(i,k) = abs(summed_intensities(1));
    end
    k
end
toc