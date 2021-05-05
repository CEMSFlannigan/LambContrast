k_mag = 398; % nm^-1
g_mags = 3.06.*ones(1,1); % nm^-1 (1-1-1)

k_unit = [0.0347794, -0.996559, 0.075241]; % tilted -4.33 degrees in alpha direction (to mimic UEM)
k_actual = k_unit*k_mag;
g_unit = [0, 0.0349736, 0.999388]; % straight perpendicular ---> ZOLZ Weiss Zone Rule
%g = g_unit.*g_mags;

angle = 4.0448; % degrees
c = sind(angle); % nm/nm

lambda = 2.5079e-3*1e-9; % m, 200 keV
theta_0 = 0.117*32^(1/3)/(200^(1/2)); % rad
beta = 0.0038; % rad
m_0 = 9.11e-31; % kg
a_0 = 5.29e-11; % m
Z = 32;
E_0 = 200*1000*1.602e-19; % J, 200 keV
light = 2.9979e8; % m/s
rho = 5.323/1000*(100)^3; % kg/m^3
N_0 = 6.022e23; % atoms/mol
A = 72.64/1000; % kg/mol
sigma = (Z*lambda*(a_0/(Z^0.33))*(1 + E_0/m_0/light^2))^2/(pi*a_0^2*(1+(beta/theta_0)^2));
frac_scatter = @(t) exp(-(N_0./A).*sigma.*rho.*t);

% Want 220 nm width

example_depth = 5000; % nm

d = 50/2;

distortion_amps = [0 1 5 10];% 1 2];

syms f(x,mag_distort) f2(x,mag_distort) curline(x) %
f(x,mag_distort) = (c/2*x+d)*(1 + mag_distort/100*sin(2*pi*x/200));%(0.005*mag_distort*(c/2*250+d)*2)*exp(-(x-250)^2/50^2) + (0.005*mag_distort*(c/2*500+d)*2)*exp(-(x-500)^2/50^2) + (0.005*mag_distort*(c/2*750+d)*2)*exp(-(x-750)^2/50^2);% + (0.005*mag_distort*(c/2*1000+d)*2)*exp(-(x-1000)^2/50^2) + (0.005*mag_distort*(c/2*1250+d)*2)*exp(-(x-1250)^2/50^2) + (0.005*mag_distort*(c/2*1500+d)*2)*exp(-(x-1500)^2/50^2) + (0.005*mag_distort*(c/2*1750+d)*2)*exp(-(x-1750)^2/50^2);% - (0.005*mag_distort*(c/2*350+50)*2)*exp(-(x-350)^2/50^2);% - 1*sin(x/10/pi); % nm
f2(x,mag_distort) = (-c/2*x-d)*(1 + mag_distort/100*sin(2*pi*x/200));%(0.005*mag_distort*(c/2*250+d)*2)*exp(-(x-250)^2/50^2) - (0.005*mag_distort*(c/2*500+d)*2)*exp(-(x-500)^2/50^2) - (0.005*mag_distort*(c/2*750+d)*2)*exp(-(x-750)^2/50^2);% - (0.005*mag_distort*(c/2*1000+d)*2)*exp(-(x-1000)^2/50^2) - (0.005*mag_distort*(c/2*1250+d)*2)*exp(-(x-1250)^2/50^2) - (0.005*mag_distort*(c/2*1500+d)*2)*exp(-(x-1500)^2/50^2) - (0.005*mag_distort*(c/2*1750+d)*2)*exp(-(x-1750)^2/50^2);% + (0.005*mag_distort*(c/2*350+50)*2)*exp(-(x-350)^2/50^2);

df = diff(f,x);
dff = diff(df,x);
df2 = diff(f2,x);
dff2 = diff(f2,x);

s_g = @(vec,g)  -dot(g,2*k_actual+g)/(2*sqrt(dot(k_actual+g,k_actual+g))*(dot(vec,k_actual+g)/(sqrt(dot(vec,vec))*sqrt(dot(k_actual+g,k_actual+g)))));

D = @(z0,vec,g) -(2*dot(vec,k_actual+g))/(pi*z0*dot(g,g)*sqrt(dot(vec,vec)))*sin(-(pi*z0*dot(g,g)*sqrt(dot(vec,vec)))/(2*dot(vec,k_actual+g)));

xvals = 0:1:900; % nm
D_arr_1 = zeros(length(distortion_amps),length(xvals));
D_arr_2 = zeros(length(distortion_amps),length(xvals));
D_sq_total = zeros(length(distortion_amps),length(xvals));
s_g_arr_1 = zeros(length(distortion_amps),length(xvals));
s_g_arr_2 = zeros(length(distortion_amps),length(xvals));
tot_strain = zeros(length(distortion_amps),length(xvals));
ex_image = cell(1,length(distortion_amps));
thick_arr = zeros(length(distortion_amps),length(xvals));
thick_arr_2 = zeros(length(distortion_amps),length(xvals));

h1 = figure;

for k = 1:length(distortion_amps)
    cur_distort = distortion_amps(k);
    
    figure(h1);
    subplot(2,4,1);hold on;
    plot(xvals,f(xvals,cur_distort)); 
    plot(xvals,f2(xvals,cur_distort));
    hold off;
    
    pause(0.001);
    
    thickness_wedge = double(f(xvals,0) - f2(xvals,0));
    thickness_distort = double(f(xvals,cur_distort) - f2(xvals,cur_distort));
    thickness_beam = zeros(length(distortion_amps),length(thickness_wedge));
    
    for i = 1:length(xvals)
        beam_slope = k_unit(2)/k_unit(1);
        curx = xvals(i);
        q = double(f(curx,cur_distort));
        curline(x) = beam_slope*(x-curx)+q;
        x2 = vpasolve(curline(x) - f2(x,cur_distort) == 0,x);%(q-invslope*curx)/invslope;
        q2 = double(f2(x2,cur_distort));
        y2 = (q2 - q)/k_unit(2)*k_unit(3);
        thickness_beam(k,i) = sqrt((x2 - curx)^2 + (q2 - q)^2 + y2^2);
    end
    
    strain = (thickness_distort-thickness_wedge)./thickness_wedge;
    
    strain_dir = [0 1 0]; % Poisson's ratio for Ge is 0.28
    
    for j = 1:length(g_mags)
        g_cur = g_unit*g_mags(j);
        
        for i = 1:length(xvals)
            curx = xvals(i);
            df_multip = (double(df(curx,cur_distort)) < 0)*1 + (double(df(curx,cur_distort)) > 0)*-1 + (double(df(curx,cur_distort)) == 0)*0;
            strain_dir(1) = abs(strain_dir(1))*df_multip;
            g_distort = g_mags(j) - g_mags(j)./(1+strain(i).*dot(strain_dir,g_unit)/(dot(strain_dir,strain_dir)*dot(g_unit,g_unit)));
            g_actual = g_cur+g_distort*g_unit;
            
            tot_strain(k,i) = sqrt(dot(strain(i).*strain_dir,strain(i).*strain_dir));
            
            %% First relrod
            q = double(f(curx,cur_distort));
            slope = double(df(curx,cur_distort));
            invslope = -1/slope;
            curline(x) = invslope*(x-curx)+q;
            x2 = vpasolve(curline(x) - f2(x,cur_distort) == 0,x);%(q-invslope*curx)/invslope;
            q2 = double(f2(x2,cur_distort));
            thickness = sqrt((x2-curx)^2+(q2-q)^2);
            thick_arr(k,i) = thickness;
            m = [-(x2-curx), q-q2, 0];
            s_g_arr_1(k,i) = s_g(m,g_actual);
            D_arr_1(k,i) = D(thickness,m,g_actual);
            
            %% Second relrod
            
            q = double(f2(curx,cur_distort));
            slope = double(df2(curx,cur_distort));
            invslope = -1/slope;
            curline(x) = invslope*(x-curx)+q;
            x2 = vpasolve(curline(x) - f(x,cur_distort) == 0,x);%(q-invslope*curx)/invslope;
            q2 = double(f(x2,cur_distort));
            thickness = sqrt((x2-curx)^2+(q2-q)^2);
            thick_arr_2(k,i) = thickness;
            m = [-(x2-curx), q-q2, 0];
            s_g_arr_2(k,i) = s_g(m,g_actual);
            D_arr_2(k,i) = D(thickness,m,g_actual);
        end
    end
    
    figure(h1);
    subplot(2,4,5);hold on;
    plot(xvals,strain,'-');
    plot(xvals,tot_strain(k,:),'--');
    hold off;

    subplot(2,4,2);hold on;
    plot(xvals,s_g_arr_1(k,:));
    plot(xvals,s_g_arr_2(k,:));
    hold off;
    
    subplot(2,4,3);hold on;
    plot(xvals,D_arr_1(k,:));
    plot(xvals,D_arr_2(k,:));
    hold off;
    
    subplot(2,4,4);hold on;
    plot(xvals,D_arr_1(k,:).^2+D_arr_2(k,:).^2);
    hold off;
    
    init_D_arr_sq = D_arr_1(1,:).^2+D_arr_2(1,:).^2;
    sum_D_arr_sq = D_arr_1(k,:).^2+D_arr_2(k,:).^2;
    
    subplot(2,4,6);hold on;
    plot(xvals,thick_arr(k,:));
    plot(xvals,thickness_beam(k,:));
    hold off;

    subplot(2,4,7:8);hold on;
    AmplitudeChange = (sum_D_arr_sq-init_D_arr_sq);
    plot(xvals,AmplitudeChange);
    ylim([-0.75 0.75]);
    hold off;

    sum_D_arr_sq = (sum_D_arr_sq-min(sum_D_arr_sq))./(max(sum_D_arr_sq) - min(sum_D_arr_sq));
    
    bf_contrast = (1-sum_D_arr_sq).*frac_scatter(thickness_beam(k,:)./1e9);
    
    cur_img = zeros(example_depth,length(xvals)+101);
    for l = 1:example_depth
        cur_img(l,1) = 0;
        cur_img(l,2:101) = ones(size(cur_img(l,2:101)));
        cur_img(l,102:end) = ones(size(cur_img(l,102:end))).*bf_contrast;
    end
    cur_img(cur_img < 0) = 0;
    cur_img(cur_img > 1) = 1;
    ex_image{k} = cur_img;
    
    figure;
    imagesc(ex_image{k});
    
    pause(0.001);
end

D_sq_total = D_arr_1.^2+D_arr_2.^2;

