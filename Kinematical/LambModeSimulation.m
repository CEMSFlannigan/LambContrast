% add script to later calculate the location of the front?

latpar = 0.5658; % nm

y_point_num = 1000; % points
yvals = [25 65]; % nm
xvals = 0:5*latpar:400; % nm
timevals = 0:1:70; % ps

yval_arr = zeros(length(xvals),y_point_num);
slope = (yvals(end) - yvals(1))/(xvals(end) - xvals(1)); % degrees
for i = 1:length(xvals)
    cury = slope*xvals(i)+yvals(1);
    yval_arr(i,:) = linspace(0,cury,y_point_num);
end

x_res = zeros(length(xvals),length(yvals),length(timevals));
y_res = zeros(length(xvals),length(yvals),length(timevals));
x_strain = zeros(length(xvals),length(yvals),length(timevals));
y_strain = zeros(length(xvals),length(yvals),length(timevals));
x_adj = zeros(size(x_res));
y_adj = zeros(size(y_res));

A = 1e-2; % nm
k = 2*pi/200; % 1/nm
omega = (5e9/1e12)*2*pi; % 1/ps
VL = 5.350; % nm/ps
VT = 3.570; % nm/ps
V = 10; % nm/ps
h = 25; % nm
alpha = 0;

p = sqrt(omega^2*(1/VL^2 - 1/V^2));
q = sqrt(omega^2*(1/VT^2 - 1/V^2));

syms x_dist(y) y_dist(y) time_space(x,t) combined_x(x,y,t) combined_y(x,y,t)

x_dist(y) = q*A*(cos(q*y+alpha) - 2*k^2/(k^2-q^2)*cos(q*h+alpha)/cos(p*h+alpha)*cos(p*y+alpha));
y_dist(y) = k*A*(sin(q*y+alpha) + 2*p*q/(k^2-q^2)*cos(q*h+alpha)/cos(p*h+alpha)*sin(p*y+alpha));
time_space(x,t) = cos(k*x-omega*t);
combined_x(x,y,t) = x_dist(y)*time_space(x,t);
combined_y(x,y,t) = y_dist(y)*time_space(x,t);

% Let's first test what happens if we do not account for the current
% y-state of the wave and only base our calculations off of the original
% specimen

for j = 1:length(yvals)
    for k = 1:length(timevals)
        cury = yvals(j);
        curt = timevals(k);
        x_res(:,j,k) = combined_x(xvals,cury,curt);
        y_res(:,j,k) = combined_y(xvals,cury,curt);
        
        x_adj(:,j,k) = xvals' + x_res(:,j,k);
        y_adj(:,j,k) = yvals(j).*ones(size(y_res(:,j,k))) + y_res(:,j,k);
    end
    disp(yvals(j));
end

for i = 1:length(xvals)
    curx = xvals(i);
    x_strain(i,:,:) = (x_adj(i,:,:) - curx.*ones(size(x_adj(i,:,:))))./latpar.*ones(size(x_adj(i,:,:)));
end

for j = 1:length(yvals)
    cury = yvals(j);
    y_strain(:,j,:) = (y_adj(:,j,:) - cury.*ones(size(y_adj(:,j,:))))./latpar.*ones(size(y_adj(:,j,:)));
end

figure; subplot(1,2,1);
for k = 1:length(timevals)
    hold on;
    plot(x_adj(:,end,k),y_adj(:,end,k));
    plot(x_adj(:,1,k),y_adj(:,1,k));
    pause(0.05);
end

subplot(1,2,2);
for k = 1:length(timevals)
    hold on;
    plot(xvals,x_res(:,end,k));
    pause(0.05);
end

figure;
for k = 1:length(timevals)
    subplot(1,2,1); hold on;
    plot(xvals,x_strain(:,end,k));
    subplot(1,2,2); hold on;
    plot(xvals,y_strain(:,end,k));
end

figure;
for k = 1:length(timevals)
    subplot(1,2,1); hold on;
    plot(yvals,x_strain(end,:,k));
    subplot(1,2,2); hold on;
    plot(yvals,y_strain(end,:,k));
end

% Need to add in width increase along the wedge!!!