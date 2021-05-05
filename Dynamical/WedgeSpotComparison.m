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
alpha_add = 0; 
k_unit = [0, -1, 0];
k_unit = k_unit/sqrt(dot(k_unit,k_unit));
k_alph_rot = [1,0,0;0,cosd(alpha_tilt),-sind(alpha_tilt);0,sind(alpha_tilt),cosd(alpha_tilt)];
k_gam_rot = [cosd(gam_tilt),0,sind(gam_tilt);0,1,0;-sind(gam_tilt),0,cosd(gam_tilt)];
k_beta_rot = [cosd(beta_tilt),-sind(beta_tilt),0;sind(beta_tilt),cosd(beta_tilt),0;0,0,1];
k_alph_add = [1,0,0;0,cosd(alpha_add),-sind(alpha_add);0,sind(alpha_add),cosd(alpha_add)];
k_beta_add = [cosd(beta_add),-sind(beta_add),0;sind(beta_add),cosd(beta_add),0;0,0,1];
k_unit = k_beta_rot*k_gam_rot*k_alph_rot*k_beta_add*k_alph_add*k_unit';

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
g_basis_3 = [1,1,0];
or_z_basis = [1, -1, 0];
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

slope = sind(1);
xR = -15:15;
zR = -15:15;
top_y = @(x) slope*x + 25;
bot_y = @(x) -slope*x - 25;

figure; hold on; 
for i = 1:length(gvec_list)
curgvec = gvec_list{i};
scatter3(curgvec(1),curgvec(2),curgvec(3));
end

plot3(linspace(-10*k_unit(1),10*k_unit(1),2), linspace(-10*k_unit(2),10*k_unit(2),2), linspace(-10*k_unit(3),10*k_unit(3),2));
plot(xR, top_y(xR));
plot(xR, bot_y(xR));