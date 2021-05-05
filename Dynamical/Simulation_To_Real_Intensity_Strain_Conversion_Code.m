final_fit_intens = @(k,log_LA) 1/2*erf(0.8169*(k.^0.06642).*(log_LA' - 2.709363243731673)) + 1/2;
% @(k,log_LA) 1/2*erf((0.8877*ones(length(log_LA),1)*k.^(0.09331)).*(log_LA'*ones(1,length(k)) - (2.919*ones(length(log_LA),1)*k.^(0.01667)))) + 1/2;

final_fit_str_xx = @(LA) 1.259346411068840e-04.*LA'*ones(1,length(k)); % (0.02465*ones(length(LA),1)*k.^(0.8708)).*LA'*ones(1,length(k))
final_fit_str_yy = @(LA) 0.004189919548932.*LA'*ones(1,length(k)); % (0.0006043*ones(length(LA),1)*k.^(-0.4419)).*LA'*ones(1,length(k))
final_fit_str_xy = @(LA) 0.001306249464878.*LA'*ones(1,length(k)); % (0.0008365*ones(length(LA),1)*k.^(-0.1452)).*LA'*ones(1,length(k))

init_guess = 0;
LA_wave_01 = fzero(@(logLA) final_fit_intens(kArrReal(1),logLA) - intensArr(1),init_guess); % final_fit_intens(kArrReal(1),logLA) - intensArr(1),init_guess)
LA_wave_02 = fzero(@(logLA) final_fit_intens(kArrReal(2),logLA) - intensArr(2),init_guess);
LA_wave_03 = fzero(@(logLA) final_fit_intens(kArrReal(3),logLA) - intensArr(3),init_guess);
LA_wave_04 = fzero(@(logLA) final_fit_intens(kArrReal(4),logLA) - intensArr(4),init_guess);
LA_wave_05 = fzero(@(logLA) final_fit_intens(kArrReal(5),logLA) - intensArr(5),init_guess);
LA_wave_06 = fzero(@(logLA) final_fit_intens(kArrReal(6),logLA) - intensArr(6),init_guess);
LA_wave_07 = fzero(@(logLA) final_fit_intens(kArrReal(7),logLA) - intensArr(7),init_guess);
LA_wave_08 = fzero(@(logLA) final_fit_intens(kArrReal(8),logLA) - intensArr(8),init_guess);
LA_wave_09 = fzero(@(logLA) final_fit_intens(kArrReal(9),logLA) - intensArr(9),init_guess);
LA_wave_10 = fzero(@(logLA) final_fit_intens(kArrReal(10),logLA) - intensArr(10),init_guess);
LA_wave_11 = fzero(@(logLA) final_fit_intens(kArrReal(11),logLA) - intensArr(11),init_guess);

LA_Arr = [LA_wave_01,LA_wave_02,LA_wave_03,LA_wave_04,LA_wave_05,LA_wave_06,LA_wave_07,LA_wave_08,LA_wave_09,LA_wave_10,LA_wave_11];
LA_Arr = 10.^(LA_Arr);
str_xx_vals = [final_fit_str_xx(LA_Arr(1)),final_fit_str_xx(LA_Arr(2)),final_fit_str_xx(LA_Arr(3)),final_fit_str_xx(LA_Arr(4)),final_fit_str_xx(LA_Arr(5)),final_fit_str_xx(LA_Arr(6)),final_fit_str_xx(LA_Arr(7)),final_fit_str_xx(LA_Arr(8)),final_fit_str_xx(LA_Arr(9)),final_fit_str_xx(LA_Arr(10)),final_fit_str_xx(LA_Arr(11))]; % 
% [final_fit_str_xx(kArrReal(1),LA_Arr(1)),final_fit_str_xx(kArrReal(2),LA_Arr(2)),final_fit_str_xx(kArrReal(3),LA_Arr(3)),final_fit_str_xx(kArrReal(4),LA_Arr(4)),final_fit_str_xx(kArrReal(5),LA_Arr(5)),final_fit_str_xx(kArrReal(6),LA_Arr(6)),final_fit_str_xx(kArrReal(7),LA_Arr(7)),final_fit_str_xx(kArrReal(8),LA_Arr(8)),final_fit_str_xx(kArrReal(9),LA_Arr(9)),final_fit_str_xx(kArrReal(10),LA_Arr(10)),final_fit_str_xx(kArrReal(11),LA_Arr(11))];
str_yy_vals = [final_fit_str_yy(LA_Arr(1)),final_fit_str_yy(LA_Arr(2)),final_fit_str_yy(LA_Arr(3)),final_fit_str_yy(LA_Arr(4)),final_fit_str_yy(LA_Arr(5)),final_fit_str_yy(LA_Arr(6)),final_fit_str_yy(LA_Arr(7)),final_fit_str_yy(LA_Arr(8)),final_fit_str_yy(LA_Arr(9)),final_fit_str_yy(LA_Arr(10)),final_fit_str_yy(LA_Arr(11))];
str_xy_vals = [final_fit_str_xy(LA_Arr(1)),final_fit_str_xy(LA_Arr(2)),final_fit_str_xy(LA_Arr(3)),final_fit_str_xy(LA_Arr(4)),final_fit_str_xy(LA_Arr(5)),final_fit_str_xy(LA_Arr(6)),final_fit_str_xy(LA_Arr(7)),final_fit_str_xy(LA_Arr(8)),final_fit_str_xy(LA_Arr(9)),final_fit_str_xy(LA_Arr(10)),final_fit_str_xy(LA_Arr(11))];

figure;
subplot(1,3,1);
plot(str_xx_vals);
subplot(1,3,2);
plot(str_yy_vals);
subplot(1,3,3);
plot(str_xy_vals);