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
k_unit = [0, -1, 0];
k_unit = k_unit/sqrt(dot(k_unit,k_unit));
k_alph_rot = [1,0,0;0,cosd(alpha_tilt),-sind(alpha_tilt);0,sind(alpha_tilt),cosd(alpha_tilt)];
k_gam_rot = [cosd(gam_tilt),0,sind(gam_tilt);0,1,0;-sind(gam_tilt),0,cosd(gam_tilt)];
k_beta_rot = [cosd(beta_tilt),-sind(beta_tilt),0;sind(beta_tilt),cosd(beta_tilt),0;0,0,1];
k_beta_add = [cosd(beta_add),-sind(beta_add),0;sind(beta_add),cosd(beta_add),0;0,0,1];
k_unit = k_beta_rot*k_gam_rot*k_alph_rot*k_beta_add*k_unit';

% Crystal lattice contstants and probing dimensions used later in the code for probing strain or x-positions
latpar = 0.5658; % nm
strainpar = sqrt(dot(1/4*[1,1,1],1/4*[1,1,1]))*latpar; % nm
V_c = latpar^3;
ExtDist = @(thet,h,k,l) pi*(V_c./1e9./1e9./1e9)*cosd(thet)/lambda/(Ge_ASF(thet)/1e10)/real(Ge_SF(h,k,l));
x_spacing = 10*latpar; % nm

% Unit cell positions
r1 = [0, 0, 0];
r2 = [1/4, 1/4, 1/4];
r3 = [0, 1/2, 1/2];
r4 = [1/2, 0, 1/2];
r5 = [1/2, 1/2, 0];
r6 = [1/4, 3/4, 3/4];
r7 = [3/4, 1/4, 3/4];
r8 = [3/4, 3/4, 1/4];

% Distance by which to increment along the path of the electron, assuming column approximation
thick_incr = 0.1; % nm
xy_dist_incr = sqrt(thick_incr^2*(1-k_unit(3)^2));
%xy_dist_incr = 0.0001; % nm
%dist_const = xy_dist_incr^2/(k_unit(1)^2+k_unit(2)^2);
%thick_incr = sqrt(xy_dist_incr^2 + dist_const^2*k_unit(3)^2); % nm

% Setting the initial shape of the wedge
yinit = 25; %[25 65]; % nm
xvals = 0:x_spacing:1200; % nm

% Setting the timepoints
timevals = 0:1:17; % ps

% Setting the shape of the wedge
angle = 2/2; % degrees!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
slope = sind(angle); % nm/nm

% Setting the shape of the undisturbed surfaces of the wedge
y_top = @(x) slope*x+yinit;
y_bot = @(x) -slope*x-yinit;

% Setting distortion parameters
LambAmp = 6;%0.005;%3;%6;%4e2;%3.75e1;%3.75; % nm
k = 0.027286040502180; %.0049693; %0.010 %0.0270; % 1/nm
omega = (35.739999668189228e9/1e12)*2*pi; % 1/ps
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

% Setting up the travel distance of each increment
k_slope = k_unit(2)/k_unit(1);
x_incr = sign(k_unit(1))*xy_dist_incr/sqrt(1+k_slope^2);
y_incr = sign(k_unit(2))*xy_dist_incr/sqrt(1+(1/k_slope)^2);
dist_travel = sqrt(abs(x_incr)^2 + abs(y_incr)^2 + abs(k_unit(3)/k_unit(2)*y_incr)^2);

tic;
i = 111;
parfor k = 1:length(timevals)
    
        
        curx = xvals(i);
        cur_start_y = start_y(curx);
        cury = cur_start_y;
        
        strain_vec = [];
        
        prevx_lamb = curx + combined_x(curx, cury, timevals(k));
        prevy_lamb = cury + combined_y(curx, cury, timevals(k));
        while abs(cury) <= abs(end_y(curx))
            
            curx_lamb = curx + combined_x(curx, cury, timevals(k));
            cury_lamb = cury + combined_y(curx, cury, timevals(k));
            
            x1 = curx - strainpar/2;
            x2 = curx + strainpar/2;
            y1 = cury - strainpar/2;
            y2 = cury + strainpar/2;
            
            uxxn1n = double(combined_x(x1, cury, timevals(k)));
            uxxn2n = double(combined_x(x2, cury, timevals(k)));
            uxyn1n = double(combined_x(curx, y1, timevals(k)));
            uxyn2n = double(combined_x(curx, y2, timevals(k)));
            uxxy11 = double(combined_x(x1, y1, timevals(k)));
            uxxy12 = double(combined_x(x1, y2, timevals(k)));
            uxxy21 = double(combined_x(x2, y1, timevals(k)));
            uxxy22 = double(combined_x(x2, y2, timevals(k)));
            
            uyyny1 = double(combined_y(curx, y1, timevals(k)));
            uyyny2 = double(combined_y(curx, y2, timevals(k)));
            uyxn1n = double(combined_y(x1, cury, timevals(k)));
            uyxn2n = double(combined_y(x2, cury, timevals(k)));
            uyxy11 = double(combined_y(x1, y1, timevals(k)));
            uyxy12 = double(combined_y(x1, y2, timevals(k)));
            uyxy21 = double(combined_y(x2, y1, timevals(k)));
            uyxy22 = double(combined_y(x2, y2, timevals(k)));
            
            strain_xx = (uxxn2n-uxxn1n)/(x2-x1);
            strain_yy = (uyyny2-uyyny1)/(y2-y1);
            strain_xyyz = ((uyxn2n - uyxn1n)/(x2-x1) + (uxyn2n - uxyn1n)/(y2-y1))/2;
            
            strain_vec = [strain_vec; strain_xx, strain_yy, strain_xyyz];
                                 
            curx = curx + x_incr;
            cury = cury + y_incr;
            
            prevx_lamb = curx_lamb;
            prevy_lamb = cury_lamb;
            
        end
        
        strain_matr{i,k} = strain_vec;
end
toc