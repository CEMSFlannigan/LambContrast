% add script to later calculate the location of the front?

latpar = 0.5658; % nm
V_c = latpar^3;
x_spacing = 10*latpar; % nm

y_point_num = 65; % points
yvals = [25 65]; %[25 65]; % nm
xvals = 0:x_spacing:2000; % nm
timevals = 0:1:100; % ps

angle = 4/2; % degrees!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
slope = sind(angle); % nm/nm

yval_arr = zeros(length(xvals),y_point_num);
for i = 1:length(xvals)
    cury = slope*xvals(i)+yvals(1);
    yval_arr(i,:) = linspace(-1*cury,cury,y_point_num);
end

x_res = zeros(length(xvals),y_point_num,length(timevals));
y_res = zeros(length(xvals),y_point_num,length(timevals));
x_strain = zeros(length(xvals)-1,y_point_num,length(timevals));
y_strain = zeros(length(xvals),y_point_num-1,length(timevals));
x_adj = zeros(size(x_res));
y_adj = zeros(size(y_res));

LambAmp = 37.5;%6;%4e2;%3.75e1;%3.75; % nm
k = 0.0270; %0.010 %0.0270; % 1/nm
omega = (35e9/1e12)*2*pi; % 1/ps
VL = 5.350; % nm/ps
VT = 3.570; % nm/ps
V = 8; % nm/ps
h = 25; % nm
alpha = 0;

p = sqrt(omega^2*(1/VL^2 - 1/V^2));
q = sqrt(omega^2*(1/VT^2 - 1/V^2));

syms x_dist(y) y_dist(y) time_space(x,t) combined_x(x,c,t,d) combined_y(x,y,t)

x_dist(y) = q*LambAmp*(cos(q*y+alpha) - 2*k^2/(k^2-q^2)*cos(q*h+alpha)/cos(p*h+alpha)*cos(p*y+alpha));
y_dist(y) = k*LambAmp*(sin(q*y+alpha) + 2*p*q/(k^2-q^2)*cos(q*h+alpha)/cos(p*h+alpha)*sin(p*y+alpha));
time_space(x,t) = cos(k*x-omega*t);
combined_x(x,y,t) = x_dist(y).*time_space(x,t);
combined_y(x,y,t) = y_dist(y).*time_space(x,t);

w1 = waitbar(0,'Calculating adjustments.');
for k = 1:length(timevals)
    parfor j = 1:y_point_num
        curt = timevals(k);
        x_res(:,j,k) = double(combined_x(xvals,yval_arr(:,j)',curt));
        y_res(:,j,k) = double(combined_y(xvals,yval_arr(:,j)',curt));
        
        x_adj(:,j,k) = xvals' + x_res(:,j,k);
        y_adj(:,j,k) = yval_arr(:,j) + y_res(:,j,k);
    end
    waitbar(k/length(timevals));
end
close(w1);

xvalsm = (xvals(2:end) - xvals(1:end-1))./2+xvals(1:end-1);
yval_arrm = zeros(length(yval_arr(:,1))-1,length(yval_arr(1,:))-1);

for i = 1:length(xvals)
    for j = 1:y_point_num - 1
        yval_arrm(i,j) = (yval_arr(i,j+1) + yval_arr(i,j))/2;
    end
end

w2 = waitbar(0,'Calculating strains.');
for k = 1:length(timevals)
    for i = 1:length(xvals)-1
        for j = 1:y_point_num
            xadj_dif = (x_adj(i+1,j,k) - x_adj(i,j,k));
            xorig_dif = (xvals(i+1) - xvals(i));
            x_strain(i,j,k) = (xadj_dif - xorig_dif)./(xorig_dif);
        end
    end
    
    for i = 1:length(xvals)
        for j = 1:y_point_num-1
            yadj_dif = (y_adj(i,j+1,k) - y_adj(i,j,k));
            yorig_dif = (yval_arr(i,j+1) - yval_arr(i,j));
            y_strain(i,j,k) = (yadj_dif - yorig_dif)./(yorig_dif);
        end
    end
    waitbar(k/length(timevals));
end
close(w2);