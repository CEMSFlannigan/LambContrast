%% Setting physical constants
lambda = 2.5079e-3*1e-9; % m, 200 keV

g_mags = 3.06; % nm^-1 (1-1-1)
g_mags_2 = 3.06; % nm^-1 (1-11)

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

% Finding the second Bragg spot basis
g_unit_2 = ([cosd(basis_angle),0,sind(basis_angle);0,1,0;-sind(basis_angle),0,cosd(basis_angle)]*g_unit')';
g_unit_2 = g_unit_2/sqrt(dot(g_unit_2,g_unit_2));

zone_axis_h = [1,1,0];
zone_axis_g = g_beta_rot*g_alph_rot*zone_axis_h'*g_mags;

gvec_list_h_init = [g_unit(1)*g_mags, g_unit_2(1)*g_mags_2, -g_unit(1)*g_mags, -g_unit_2(1)*g_mags_2];
gvec_list_k_init = [g_unit(2)*g_mags, g_unit_2(2)*g_mags_2, -g_unit(2)*g_mags, -g_unit_2(2)*g_mags_2];
gvec_list_l_init = [g_unit(3)*g_mags, g_unit_2(3)*g_mags_2, -g_unit(3)*g_mags, -g_unit_2(3)*g_mags_2];

hvec_list_h_init = [g_basis_1(1), g_basis_2(1), -g_basis_1(1), -g_basis_2(1)];
hvec_list_k_init = [g_basis_1(2), g_basis_2(2), -g_basis_1(2), -g_basis_2(2)];
hvec_list_l_init = [g_basis_1(3), g_basis_2(3), -g_basis_1(3), -g_basis_2(3)];

gvec_list_h_master = [g_unit(1)*g_mags, g_unit_2(1)*g_mags_2, -g_unit(1)*g_mags, -g_unit_2(1)*g_mags_2];
gvec_list_k_master = [g_unit(2)*g_mags, g_unit_2(2)*g_mags_2, -g_unit(2)*g_mags, -g_unit_2(2)*g_mags_2];
gvec_list_l_master = [g_unit(3)*g_mags, g_unit_2(3)*g_mags_2, -g_unit(3)*g_mags, -g_unit_2(3)*g_mags_2];

hvec_list_h_master = [g_basis_1(1), g_basis_2(1), -g_basis_1(1), -g_basis_2(1)];
hvec_list_k_master = [g_basis_1(2), g_basis_2(2), -g_basis_1(2), -g_basis_2(2)];
hvec_list_l_master = [g_basis_1(3), g_basis_2(3), -g_basis_1(3), -g_basis_2(3)];

zero_list_h = [];
zero_list_k = [];
zero_list_l = [];

angled = 3.1;
%angled = 2.5;

angles = @(x,y,z) asind(lambda*1e9/2*sqrt(x.*x+y.*y+z.*z));
numLaue = 0;
curangle = 0;
while(curangle < angled)
    numLaue = numLaue + 1;
    cur_Laue_g = numLaue*zone_axis_g;
    curangle = angles(cur_Laue_g(1), cur_Laue_g(2), cur_Laue_g(3));
end

totcount = 0;

for i = 0:1
    
    gvec_list_h = gvec_list_h_init + zone_axis_g(1)*i;
    gvec_list_k = gvec_list_k_init + zone_axis_g(2)*i;
    gvec_list_l = gvec_list_l_init + zone_axis_g(3)*i;
    
    hvec_list_h = hvec_list_h_init + zone_axis_h(1)*i;
    hvec_list_k = hvec_list_k_init + zone_axis_h(2)*i;
    hvec_list_l = hvec_list_l_init + zone_axis_h(3)*i;
    
    [gvec_list_h,gvec_list_k,gvec_list_l,hvec_list_h,hvec_list_k,hvec_list_l,count] = RecipSpotGeneratorAngle(angled,4,gvec_list_h,gvec_list_k,gvec_list_l,hvec_list_h,hvec_list_k,hvec_list_l, g_unit, g_unit_2, g_mags, g_mags_2, g_basis_1, g_basis_2);

    for j = 1:count
        gvec_list_h_master = [gvec_list_h_master gvec_list_h(j)];
        gvec_list_k_master = [gvec_list_k_master gvec_list_k(j)];
        gvec_list_l_master = [gvec_list_l_master gvec_list_l(j)];
        
        hvec_list_h_master = [hvec_list_h_master hvec_list_h(j)];
        hvec_list_k_master = [hvec_list_k_master hvec_list_k(j)];
        hvec_list_l_master = [hvec_list_l_master hvec_list_l(j)];
    end
    
    totcount = totcount + count;
    
end

hvec_list = cell(1,length(hvec_list_h_master));
gvec_list = cell(1,length(hvec_list_h_master));
zero_h_list = cell(1,length(zero_list_h));

for i = 1:length(hvec_list_h_master)
    
    if (hvec_list_h_master(i)+hvec_list_k_master(i)) == 0
        zero_list_h = [zero_list_h hvec_list_h_master(i)];
        zero_list_k = [zero_list_k hvec_list_k_master(i)];
        zero_list_l = [zero_list_l hvec_list_l_master(i)];
    end
    hvec_list{i} = [hvec_list_h_master(i) hvec_list_k_master(i) hvec_list_l_master(i)];
    gvec_list{i} = [gvec_list_h_master(i) gvec_list_k_master(i) gvec_list_l_master(i)];
    
end

for i = 1:length(zero_list_h)
    zero_h_list{i} = [zero_list_h(i) zero_list_k(i) zero_list_l(i)];
end

scatter3(hvec_list_h_master, hvec_list_k_master, hvec_list_l_master);
pbaspect([1 1 1]);
xlabel('h');
ylabel('k');
zlabel('l');