% Strain and Contrast Analysis Tool

k_arr = [0.00497, 0.0126, 0.0273];
k_arr_fit = fliplr(k_arr);
k_arr_ext = k_arr(1):0.0001:k_arr(end);
LA_arr = [0.005, 0.05, 0.5, 3, 6];
log_LA_arr = log(LA_arr)/log(10);
LA_arr_ext = LA_arr(1):0.001:LA_arr(end);
log_LA_arr_ext = log(LA_arr_ext);

intens_arr = 1/2*[6.07978262049946e-05,9.48144108329774e-05,0.000139429497182042;0.000607977775656021,0.000950850545870019,0.00139428776040600;0.00607913708402400,0.00951055305222400,0.0139358660196960;0.0361306349486370,0.0556381722348260,0.0785859325216220;0.0692628647621960,0.0996080926979300,0.124454466162722];
str_xx_arr = abs([9.84200000000000e-07,3.20050000000000e-07,-5.70420000000000e-07;9.84200000000000e-06,3.20050000000000e-06,-5.70420000000000e-06;9.84200000000000e-05,3.20050000000000e-05,-5.70420000000000e-05;0.000590520000000000,0.000192030000000000,-0.000342250000000000;0.00120000000000000,0.000384060000000000,-0.000684500000000000]);
str_yy_arr = abs([-1.37710000000000e-05,-2.14510000000000e-05,-2.77030000000000e-05;-0.000137710000000000,-0.000214510000000000,-0.000277030000000000;-0.00140000000000000,-0.00210000000000000,-0.00280000000000000;-0.00830000000000000,-0.0129000000000000,-0.0166000000000000;-0.0165000000000000,-0.0257000000000000,-0.0332000000000000]);
str_xy_arr = abs([-5.27770000000000e-06,-7.13590000000000e-06,-7.15870000000000e-06;-5.27770000000000e-05,-7.13590000000000e-05,-7.15870000000000e-05;-0.000527770000000000,-0.000713590000000000,-0.000715870000000000;-0.00320000000000000,-0.00430000000000000,-0.00430000000000000;-0.00630000000000000,-0.00860000000000000,-0.00860000000000000]);

% figure;
% subplot(2,2,1);
% surf(intens_arr);
% subplot(2,2,2);
% surf(str_x_arr);
% subplot(2,2,3);
% surf(str_y_arr);
% subplot(2,2,4);
% surf(str_z_arr);

syms a b x;

f_atan_ft = fittype(@(a,b,x) 1/pi*atan(a * (x - b)) + 1/2);
f_atan = @(a,b,x) 1/pi*atan(a * (x - b)) + 1/2;
f_hyper_ft = fittype(@(a,b,x) 1/2*tanh(a * (x - b)) + 1/2);
f_hyper = @(a,b,x) 1/2*tanh(a * (x - b)) + 1/2;
f_erf_ft = fittype(@(a,b,x) 1/2*erf(a * (x - b)) + 1/2);
f_erf = @(a,b,x) 1/2*erf(a * (x - b)) + 1/2;
f_exp_ft = fittype(@(a,b,x) a*exp(b*x));
f_ln_ft = fittype(@(a,b,x) a*x + b);

[fit_erf_1, gof_erf_1] = fit(log_LA_arr', intens_arr(:,1), f_erf_ft, 'Lower',[0,0],'Upper',[Inf,Inf], 'StartPoint',[1 1]);
[fit_erf_2, gof_erf_2] = fit(log_LA_arr', intens_arr(:,2), f_erf_ft, 'Lower',[0,0],'Upper',[Inf,Inf], 'StartPoint',[1 1]);
[fit_erf_3, gof_erf_3] = fit(log_LA_arr', intens_arr(:,3), f_erf_ft, 'Lower',[0,0],'Upper',[Inf,Inf], 'StartPoint',[1 1]);

[s_x_fit_exp_1, s_x_gof_exp_1] = fit(LA_arr', str_xx_arr(:,1), f_ln_ft, 'StartPoint',[1 1]);
[s_x_fit_exp_2, s_x_gof_exp_2] = fit(LA_arr', str_xx_arr(:,2), f_ln_ft, 'StartPoint',[1 1]);
[s_x_fit_exp_3, s_x_gof_exp_3] = fit(LA_arr', str_xx_arr(:,3), f_ln_ft, 'StartPoint',[1 1]);

[s_y_fit_exp_1, s_y_gof_exp_1] = fit(LA_arr', str_yy_arr(:,1), f_ln_ft, 'StartPoint',[1 1]);
[s_y_fit_exp_2, s_y_gof_exp_2] = fit(LA_arr', str_yy_arr(:,2), f_ln_ft, 'StartPoint',[1 1]);
[s_y_fit_exp_3, s_y_gof_exp_3] = fit(LA_arr', str_yy_arr(:,3), f_ln_ft, 'StartPoint',[1 1]);

[s_z_fit_exp_1, s_z_gof_exp_1] = fit(LA_arr', str_xy_arr(:,1), f_ln_ft, 'StartPoint',[1 1]);
[s_z_fit_exp_2, s_z_gof_exp_2] = fit(LA_arr', str_xy_arr(:,2), f_ln_ft, 'StartPoint',[1 1]);
[s_z_fit_exp_3, s_z_gof_exp_3] = fit(LA_arr', str_xy_arr(:,3), f_ln_ft, 'StartPoint',[1 1]);

% figure;
% subplot(2,2,1);
% scatter(log_LA_arr, intens_arr(:,1),'b');
% hold on;
% scatter(log_LA_arr, intens_arr(:,2),'g');
% scatter(log_LA_arr, intens_arr(:,3),'r');
% plot(-5:0.001:10, fit_erf_1(-5:0.001:10),'b');
% plot(-5:0.001:10, fit_erf_2(-5:0.001:10),'g');
% plot(-5:0.001:10, fit_erf_3(-5:0.001:10),'r');
% 
% subplot(2,2,2);
% scatter(log_LA_arr, str_x_arr(:,1),'b');
% hold on;
% scatter(log_LA_arr, str_x_arr(:,2),'g');
% scatter(log_LA_arr, str_x_arr(:,3),'r');
% plot(-5:0.001:2, s_x_fit_exp_1(-5:0.001:2),'b');
% plot(-5:0.001:2, s_x_fit_exp_2(-5:0.001:2),'g');
% plot(-5:0.001:2, s_x_fit_exp_3(-5:0.001:2),'r');
% 
% subplot(2,2,3);
% scatter(log_LA_arr, str_y_arr(:,1),'b');
% hold on;
% scatter(log_LA_arr, str_y_arr(:,2),'g');
% scatter(log_LA_arr, str_y_arr(:,3),'r');
% plot(-5:0.001:2, s_y_fit_exp_1(-5:0.001:2),'b');
% plot(-5:0.001:2, s_y_fit_exp_2(-5:0.001:2),'g');
% plot(-5:0.001:2, s_y_fit_exp_3(-5:0.001:2),'r');
% 
% subplot(2,2,4);
% scatter(log_LA_arr, str_z_arr(:,1),'b');
% hold on;
% scatter(log_LA_arr, str_z_arr(:,2),'g');
% scatter(log_LA_arr, str_z_arr(:,3),'r');
% plot(-5:0.001:2, s_z_fit_exp_1(-5:0.001:2),'b');
% plot(-5:0.001:2, s_z_fit_exp_2(-5:0.001:2),'g');
% plot(-5:0.001:2, s_z_fit_exp_3(-5:0.001:2),'r');

final_fit_intens = @(k,log_LA) 1/2*erf(0.8169*(k.^0.06642).*(log_LA'*ones(1,length(k)) - 2.709363243731673)) + 1/2;

final_fit_str_xx = @(k,LA) 1.259346411068840e-04.*LA'.*ones(1,length(k)); % (0.02465*ones(length(LA),1).*k.^(0.8708)).*LA'.*ones(1,length(k))
final_fit_str_yy = @(k,LA) 0.004189919548932.*LA'.*ones(1,length(k));
final_fit_str_xy = @(k,LA) 0.001306249464878.*LA'.*ones(1,length(k));

final_data = final_fit_intens(k_arr_ext,log_LA_arr_ext);
final_data_str_xx = final_fit_str_xx(k_arr_ext,LA_arr_ext);
final_data_str_yy = final_fit_str_yy(k_arr_ext,LA_arr_ext);
final_data_str_xy = final_fit_str_xy(k_arr_ext,LA_arr_ext);

figure;
subplot(2,2,1);
surf(ones(length(log_LA_arr),1)*k_arr_fit,log_LA_arr'*ones(1,length(k_arr_fit)), intens_arr);
hold on;
surf(ones(length(log_LA_arr_ext),1)*k_arr_ext,log_LA_arr_ext'*ones(1,length(k_arr_ext)),final_data,'edgecolor', 'none');
subplot(2,2,2);
surf(ones(length(LA_arr_ext),1)*k_arr_ext,LA_arr_ext'*ones(1,length(k_arr_ext)),final_data_str_xx,'edgecolor', 'none');
subplot(2,2,3);
surf(ones(length(LA_arr_ext),1)*k_arr_ext,LA_arr_ext'*ones(1,length(k_arr_ext)),final_data_str_yy,'edgecolor', 'none');
subplot(2,2,4);
surf(ones(length(LA_arr_ext),1)*k_arr_ext,LA_arr_ext'*ones(1,length(k_arr_ext)),final_data_str_xy,'edgecolor', 'none');