% %% Proceed to contrast analysis here:
% 
% %% Declare variables and constants
% k_mag = 398; % nm^-1
% g_mags = 3.06; % nm^-1 (1-1-1)
% g_mags_2 = 3.06; % nm^-1 (1-11)
% 
% simspots = 400;
% spot_angles = zeros(1,simspots);
% ext_dist_arr = zeros(1,simspots);
% SF_arr = zeros(1,simspots);
% gvec_list = cell(1,simspots);
% hvec_list = cell(1,simspots);
% 
% g_basis_1 = [1,-1,-1];
% g_basis_2 = [1,-1,1];
% basis_angle = acosd(dot(g_basis_1,g_basis_2)/sqrt(sum(g_basis_1.^2)*sum(g_basis_2.^2)));
% 
% % alpha = 10;%-4.33;
% % k_unit = [0.0349, -0.999, 0]; % tilted -4.33 degrees in alpha direction (to mimic UEM)
% % k_unit = k_unit/sqrt(dot(k_unit,k_unit));
% % k_flat_rot = [cosd(-2),-sind(-2),0;sind(-2),cosd(-2),0;0,0,1];
% % k_alph_rot = [cosd(alpha),0,sind(alpha);0,1,0;-sind(alpha),0,cosd(alpha)];
% % k_flat_rerot = [cosd(2),-sind(2),0;sind(2),cosd(2),0;0,0,1];
% % k_unit = k_flat_rerot*k_alph_rot*k_flat_rot*k_unit';
% % k_actual = k_unit*k_mag;
% 
% % k_unit = [0.0353, -1, -0.0699]; % tilted -4.33 degrees in alpha direction (to mimic UEM)
% % k_unit = k_unit/sqrt(dot(k_unit,k_unit));
% % k_actual = k_unit*k_mag;
% 
% g_unit = [0.999847695156391	0.0174524064372835	0]; % 0.0353 straight perpendicular ---> ZOLZ Weiss Zone Rule
% g_unit = g_unit/sqrt(dot(g_unit,g_unit));
% 
% g_unit_2 = ([cosd(basis_angle),0,sind(basis_angle);0,1,0;-sind(basis_angle),0,cosd(basis_angle)]*g_unit')';
% g_unit_2 = g_unit_2/sqrt(dot(g_unit_2,g_unit_2));
% g = g_unit.*g_mags;
% 
% gvec_list{1} = g_unit.*g_mags;
% gvec_list{2} = g_unit_2.*g_mags_2;
% gvec_list{3} = -g_unit.*g_mags;
% gvec_list{4} = -g_unit_2.*g_mags_2;
% 
% hvec_list{1} = g_basis_1;
% hvec_list{2} = g_basis_2;
% hvec_list{3} = -g_basis_1;
% hvec_list{4} = -g_basis_2;
% 
% [gvec_list,hvec_list] = RecipSpotGenerator(4,gvec_list,simspots,hvec_list);
% 
% lambda = 2.5079e-3*1e-9; % m, 200 keV
% vel = 2.0844e8; % m/s, 200 keV
% theta_0 = 0.117*32^(1/3)/(200^(1/2)); % rad
% rel_mass = 1.2674e-30;% kg
% beta = 0.0038; % rad
% m_0 = 9.11e-31; % kg
% a_0 = 5.29e-11; % m
% Z = 32;
% E_0 = 200*1000*1.602e-19; % J, 200 keV
% light = 2.9979e8; % m/s
% rho = 5.323/1000*(100)^3; % kg/m^3
% N_0 = 6.022e23; % atoms/mol
% A = 72.64/1000; % kg/mol
% sigma = (Z*lambda*(a_0/(Z^0.33))*(1 + E_0/m_0/light^2))^2/(pi*a_0^2*(1+(beta/theta_0)^2));
% %frac_scatter = @(t) exp(-(N_0./A).*sigma.*rho.*t);
% 
% e_charge = 1.602e-19; % C
% permitiv = 8.8541878e-12; % C/V/m
% planck = 6.626e-34; % Js
% 
% % For T = 293 K, 200keV
% AFFa1 = -0.0008;
% AFFa2 = -0.0281;
% AFFa3 = 0.1559;
% AFFa4 = 0.0873;
% AFFa5 = 0.0110;
% AFFb1 = 0.0002;
% AFFb2 = 0.1071;
% AFFb3 = 0.5417;
% AFFb4 = 1.6461;
% AFFb5 = 8.4617;
% DW = 0.6000;
% Beta100 = (1-(1+1.9579341e-3*100)^(-2))^(-1);
% Beta200 = (1-(1+1.9579341e-3*200)^(-2))^(-1);
% Ge_AFF = @(thet) (AFFa1*exp(-AFFb1*(sind(thet)/(lambda*1e10)).^2) + AFFa2*exp(-AFFb2*(sind(thet)/(lambda*1e10)).^2) + AFFa3*exp(-AFFb3*(sind(thet)/(lambda*1e10)).^2) + AFFa4*exp(-AFFb4*(sind(thet)/(lambda*1e10)).^2) + AFFa5*exp(-AFFb5*(sind(thet)/(lambda*1e10)).^2))*Beta100/Beta200;
% 
% r1 = [0, 0, 0];
% r2 = [1/4, 1/4, 1/4];
% r3 = [0, 1/2, 1/2];
% r4 = [1/2, 0, 1/2];
% r5 = [1/2, 1/2, 0];
% r6 = [1/4, 3/4, 3/4];
% r7 = [3/4, 1/4, 3/4];
% r8 = [3/4, 3/4, 1/4];
% 
% Ge_SF = @(ha,ka,la,r1a,r2a,r3a,r4a,r5a,r6a,r7a,r8a) exp(-2*pi*sqrt(-1)*(ha*r1a(1) + ka*r1a(2) + la*r1a(3))) + exp(-2*pi*sqrt(-1)*(ha*r2a(1) + ka*r2a(2) + la*r2a(3))) + exp(-2*pi*sqrt(-1)*(ha*r3a(1) + ka*r3a(2) + la*r3a(3))) + exp(-2*pi*sqrt(-1)*(ha*r4a(1) + ka*r4a(2) + la*r4a(3))) + exp(-2*pi*sqrt(-1)*(ha*r5a(1) + ka*r5a(2) + la*r5a(3))) + exp(-2*pi*sqrt(-1)*(ha*r6a(1) + ka*r6a(2) + la*r6a(3))) + exp(-2*pi*sqrt(-1)*(ha*r7a(1) + ka*r7a(2) + la*r7a(3))) + exp(-2*pi*sqrt(-1)*(ha*r8a(1) + ka*r8a(2) + la*r8a(3)));
% 
% absorb = Ge_SF(0,0,0,r1,r2,r3,r4,r5,r6,r7,r8)*47.87801/(V_c*10^3)*(Ge_AFF(0))*2*rel_mass*e_charge/(planck)^2./(sqrt(dot(k_actual',k_actual'))*1e9);
% frac_scatter = @(t) exp(-2*pi*(t*1e-9)*absorb);
% 
% s_g = @(vec,g)  -dot(g,2*k_actual'+g)/(2*sqrt(dot(k_actual'+g,k_actual'+g))*(dot(vec,k_actual'+g)/(sqrt(dot(vec,vec))*sqrt(dot(k_actual'+g,k_actual'+g)))));
% 
% D = @(z0,vec,g) sin(pi*s_g(vec,g)*z0)/(pi*s_g(vec,g));
% 
% SFa1 = 0.2135;
% SFa2 = 0.9761;
% SFa3 = 1.6555;
% SFa4 = 2.8938;
% SFa5 = 1.6356;
% SFb1 = 0.0989;
% SFb2 = 0.9845;
% SFb3 = 4.5527;
% SFb4 = 21.5563;
% SFb5 = 70.3903;
% Ge_ASF = @(thet) (1-(vel/light)^2)^(-1/2)*(SFa1*exp(-SFb1*(sind(thet)/(lambda*1e10)).^2) + SFa2*exp(-SFb2*(sind(thet)/(lambda*1e10)).^2) + SFa3*exp(-SFb3*(sind(thet)/(lambda*1e10)).^2) + SFa4*exp(-SFb4*(sind(thet)/(lambda*1e10)).^2) + SFa5*exp(-SFb5*(sind(thet)/(lambda*1e10)).^2));
% Ge_SF = @(h,k,l) (1+exp(sqrt(-1)*pi/2*(h+k+l)))*(1+exp(-sqrt(-1)*pi*(h+k)) + exp(-sqrt(-1)*pi*(h+l)) + exp(-sqrt(-1)*pi*(k+l)));
% ExtDist = @(thet,h,k,l) pi*(V_c./1e9./1e9./1e9)*cosd(thet)/lambda/(Ge_ASF(thet)/1e10)/abs(Ge_SF(h,k,l));
% 
% I_S = @(thet,h,k,l,z0,vec,g) (D(z0,vec,g).^2)./((ExtDist(thet,h,k,l)*1e9).^2);
% 
% for i = 1:simspots
%     curhvec = hvec_list{i};
%     spot_angles(i) = asind(lambda/2*sqrt(dot(gvec_list{i},gvec_list{i}))*1e9);
%     scatter3(curhvec(1),curhvec(2),curhvec(3));
%     hold on;
%     SF_arr(i) = Ge_SF(curhvec(1),curhvec(2),curhvec(3));
%     ASF_arr(i) = Ge_ASF(spot_angles(i));
%     ext_dist_arr(i) = ExtDist(spot_angles(i),curhvec(1),curhvec(2),curhvec(3));
% end
% 
% g_ave_arr = zeros(simspots,length(xvals),length(timevals));
% D_arr_1 = zeros(simspots,length(xvals),length(timevals));
% D_arr_2 = zeros(simspots,length(xvals),length(timevals));
% D_sq_total = zeros(simspots,length(xvals),length(timevals));
% s_g_arr_1 = zeros(simspots,length(xvals),length(timevals));
% s_g_arr_2 = zeros(simspots,length(xvals),length(timevals));
% I_S_arr = zeros(simspots,length(xvals),length(timevals));
% 
% example_depth = 5000; % nm
% ex_image = cell(1,length(timevals));
% 
% %% Now actually calculating the contrast
% 
% strain_x_dir = [1 0 0];
% strain_y_dir = [0 1 0];
% 
% w2 = waitbar(0,'Processing image contrast.');
% for k = 1:length(timevals)
%     parfor i = 1:length(xvals)
%         m_top = [-1 -top_inv_slopes(i,k) 0];
%         m_top = m_top/sqrt(dot(m_top,m_top));
%         m_bot = [-1 -bot_inv_slopes(i,k) 0];
%         m_bot = m_bot/sqrt(dot(m_bot,m_bot));
%         
%         cur_strain = strain_x_dir.*ave_prop_strain_x(i,k) + strain_y_dir.*ave_prop_strain_y(i,k);
%         tot_strain = sqrt(dot(cur_strain,cur_strain));
%         
%         for p = 1:simspots
%             curg_vec = gvec_list{p};
%             curh_vec = hvec_list{p};
%             
%             curg_mag = sqrt(dot(curg_vec,curg_vec));
%             g_distort = curg_mag - curg_mag./(1+tot_strain.*dot(cur_strain,curg_vec)/(sqrt(dot(cur_strain,cur_strain))*curg_mag));
%             
%             g_actual = curg_vec+g_distort*curg_vec./curg_mag;
%             
%             g_ave_arr(p,i,k) = g_distort;
%             
%             % top
%             s_g_arr_1(p,i,k) = s_g(m_top,g_actual);
%             D_arr_1(p,i,k) = D(top_surf_thickness(i,k),m_top,g_actual);
%             
%             % bot
%             s_g_arr_2(p,i,k) = s_g(m_bot,g_actual);
%             D_arr_2(p,i,k) = D(bot_surf_thickness(i,k),m_bot,g_actual);
%             I_S_arr(p,i,k) = I_S(spot_angles(p),curh_vec(1),curh_vec(2),curh_vec(3),top_surf_thickness(i,k),m_top,g_actual) + I_S(spot_angles(p),curh_vec(1),curh_vec(2),curh_vec(3),bot_surf_thickness(i,k),m_bot,g_actual);
%         end
%     end
%     waitbar(k/length(timevals));
% end
% close(w2);

look_length = round(length(xvals)/2)+2;
bf_contrast = zeros(look_length,length(timevals));
mag_I_sq_arr = zeros(look_length,length(timevals));


for k = 1:length(timevals)
    %     if simspots == 1
    %         sum_D_arr_sq = D_arr_1(1,:,k).^2+D_arr_2(1,:,k).^2;
    %
    %         min_D_sq = 0;
    %         max_D_sq = max(prop_thickness(:,k))^2;
    %
    %         sum_D_arr_sq = (sum_D_arr_sq-2*min_D_sq)./(max_D_sq*2 - min_D_sq*2);
    %
    %         bf_contrast(:,k) = (1-sum_D_arr_sq).*frac_scatter(prop_thickness(:,k)'./1e9);
    %     else
    sum_D_arr_sq = sum(D_arr_1(1:simspots,1:look_length,k).^2+D_arr_2(1:simspots,1:look_length,k).^2,1);
    
    mag_I_sq_arr(:,k) = sum(I_S_arr(1:simspots,1:look_length,k),1) + ones(size(I_S_arr(1,1:look_length,k)));
    
    min_D_sq = 0;
    max_D_sq = max(prop_thickness(1:look_length,k))^2;
    
    sum_D_arr_sq = (sum_D_arr_sq-2*min_D_sq)./(max_D_sq*2 - min_D_sq*2);
    
    bf_contrast(:,k) = ones(size(mag_I_sq_arr(:,k)))./mag_I_sq_arr(:,k).*frac_scatter(prop_thickness(1:look_length,k));
end

%bf_contrast = (bf_contrast-min(min(bf_contrast)))./(max(max(bf_contrast)) - min(min(bf_contrast)));

max_bf_contrast = max(max(bf_contrast));
min_bf_contrast = min(min(bf_contrast));

for k = 1:length(timevals)
    cur_img = zeros(example_depth,round((look_length)*x_spacing)+3);%);
    append_img = ones(size(cur_img(1,4:end)));
    
     for i = 1:look_length%length(xvals)
         append_img(1,round((i-1)*x_spacing)+1:(round((i)*x_spacing))+1) = bf_contrast(i,k);
     end
    
    for l = 1:example_depth
        cur_img(l,1) = min_bf_contrast;
        cur_img(l,2) = max_bf_contrast;
        cur_img(l,3:end) = append_img(1,:);
    end
    cur_img(cur_img < 0) = 0;
    cur_img(cur_img > 1) = 1;
    ex_image{k} = cur_img;
end