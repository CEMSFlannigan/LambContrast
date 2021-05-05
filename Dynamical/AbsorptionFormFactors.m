A1 = [-0.0099, -0.0099, -0.0099, -0.0133, -0.0171, -0.0199, -0.0007, -0.0008, -0.0008, -0.0008];
A2 = [0.0506, 0.0514, 0.0523, 0.0560, 0.0790, 0.1043, -0.0245, -0.0278, -0.0281, -0.0283];
A3 = [0.0343, 0.0351, 0.0355, 0.0403, 0.0515, 0.0341, 0.1386, 0.1551, 0.1559, 0.1558];
A4 = [0.0244, 0.0238, 0.0233, 0.0324, 0.0356, 0.0547, 0.0829, 0.0868, 0.0873, 0.0879];
A5 = [0.0049, 0.0044, 0.0044, 0.0066, 0.0066, 0.0100, 0.0110, 0.0109, 0.0110, 0.0115];
B1 = [0.0267, 0.0267, 0.0268, 0.0369, 0.0517, 0.0647, 0.0002, 0.0004, 0.0002, 0.0003];
B2 = [0.1513, 0.1536, 0.1566, 0.1776, 0.2715, 0.3758, 0.0903, 0.1058, 0.1071, 0.1078];
B3 = [0.4689, 0.4845, 0.4993, 0.4749, 0.6888, 0.7795, 0.4719, 0.5377, 0.5417, 0.5433];
B4 = [1.3187, 1.3795, 1.4088, 1.2870, 1.6476, 1.5093, 1.4956, 1.6436, 1.6461, 1.6336];
B5 = [5.0758, 5.4966, 5.4880, 5.2629, 7.0253, 6.4112, 7.6944, 8.5279, 8.4617, 8.1604];
DW = [0.1597, 0.1611, 0.1632, 0.2049, 0.3073, 0.4097, 0.5120, 0.5938, 0.6000, 0.6041];
T = [75, 77, 80, 100, 150, 200, 250, 290, 293, 295];

s_range = 0:0.0001:6;

Beta100 = (1-(1+1.9579341e-3*100)^(-2))^(-1);
Beta200 = (1-(1+1.9579341e-3*200)^(-2))^(-1);

plot(s_range, A1(1)*exp(-B1(1)*s_range.^2) + A2(1)*exp(-B2(1)*s_range.^2) + A3(1)*exp(-B3(1)*s_range.^2) + A4(1)*exp(-B4(1)*s_range.^2) + A5(1)*exp(-B5(1)*s_range.^2));
hold on;
plot(s_range, A1(2)*exp(-B1(2)*s_range.^2) + A2(2)*exp(-B2(2)*s_range.^2) + A3(2)*exp(-B3(2)*s_range.^2) + A4(2)*exp(-B4(2)*s_range.^2) + A5(2)*exp(-B5(2)*s_range.^2));
plot(s_range, A1(3)*exp(-B1(3)*s_range.^2) + A2(3)*exp(-B2(3)*s_range.^2) + A3(3)*exp(-B3(3)*s_range.^2) + A4(3)*exp(-B4(3)*s_range.^2) + A5(3)*exp(-B5(3)*s_range.^2));
plot(s_range, A1(4)*exp(-B1(4)*s_range.^2) + A2(1)*exp(-B2(4)*s_range.^2) + A3(1)*exp(-B3(4)*s_range.^2) + A4(1)*exp(-B4(4)*s_range.^2) + A5(1)*exp(-B5(4)*s_range.^2));
plot(s_range, A1(5)*exp(-B1(5)*s_range.^2) + A2(1)*exp(-B2(5)*s_range.^2) + A3(1)*exp(-B3(5)*s_range.^2) + A4(1)*exp(-B4(5)*s_range.^2) + A5(1)*exp(-B5(5)*s_range.^2));
plot(s_range, A1(6)*exp(-B1(6)*s_range.^2) + A2(1)*exp(-B2(6)*s_range.^2) + A3(1)*exp(-B3(6)*s_range.^2) + A4(1)*exp(-B4(6)*s_range.^2) + A5(1)*exp(-B5(6)*s_range.^2));
plot(s_range, A1(7)*exp(-B1(7)*s_range.^2) + A2(1)*exp(-B2(7)*s_range.^2) + A3(1)*exp(-B3(7)*s_range.^2) + A4(1)*exp(-B4(7)*s_range.^2) + A5(1)*exp(-B5(7)*s_range.^2));
plot(s_range, A1(8)*exp(-B1(8)*s_range.^2) + A2(1)*exp(-B2(8)*s_range.^2) + A3(1)*exp(-B3(8)*s_range.^2) + A4(1)*exp(-B4(8)*s_range.^2) + A5(1)*exp(-B5(8)*s_range.^2));
plot(s_range, A1(9)*exp(-B1(9)*s_range.^2) + A2(1)*exp(-B2(9)*s_range.^2) + A3(1)*exp(-B3(9)*s_range.^2) + A4(1)*exp(-B4(9)*s_range.^2) + A5(1)*exp(-B5(9)*s_range.^2));
plot(s_range, A1(10)*exp(-B1(10)*s_range.^2) + A2(1)*exp(-B2(10)*s_range.^2) + A3(1)*exp(-B3(10)*s_range.^2) + A4(1)*exp(-B4(10)*s_range.^2) + A5(1)*exp(-B5(10)*s_range.^2));

figure;
plot(s_range, (A1(1)*exp(-B1(1)*s_range.^2) + A2(1)*exp(-B2(1)*s_range.^2) + A3(1)*exp(-B3(1)*s_range.^2) + A4(1)*exp(-B4(1)*s_range.^2) + A5(1)*exp(-B5(1)*s_range.^2))./exp(-DW(1)/2*s_range.^2));
hold on;
plot(s_range, (A1(2)*exp(-B1(2)*s_range.^2) + A2(2)*exp(-B2(2)*s_range.^2) + A3(2)*exp(-B3(2)*s_range.^2) + A4(2)*exp(-B4(2)*s_range.^2) + A5(2)*exp(-B5(2)*s_range.^2))./exp(-DW(2)/2*s_range.^2));
plot(s_range, (A1(3)*exp(-B1(3)*s_range.^2) + A2(3)*exp(-B2(3)*s_range.^2) + A3(3)*exp(-B3(3)*s_range.^2) + A4(3)*exp(-B4(3)*s_range.^2) + A5(3)*exp(-B5(3)*s_range.^2))./exp(-DW(3)/2*s_range.^2));
plot(s_range, (A1(4)*exp(-B1(4)*s_range.^2) + A2(1)*exp(-B2(4)*s_range.^2) + A3(1)*exp(-B3(4)*s_range.^2) + A4(1)*exp(-B4(4)*s_range.^2) + A5(1)*exp(-B5(4)*s_range.^2))./exp(-DW(4)/2*s_range.^2));
plot(s_range, (A1(5)*exp(-B1(5)*s_range.^2) + A2(1)*exp(-B2(5)*s_range.^2) + A3(1)*exp(-B3(5)*s_range.^2) + A4(1)*exp(-B4(5)*s_range.^2) + A5(1)*exp(-B5(5)*s_range.^2))./exp(-DW(5)/2*s_range.^2));
plot(s_range, (A1(6)*exp(-B1(6)*s_range.^2) + A2(1)*exp(-B2(6)*s_range.^2) + A3(1)*exp(-B3(6)*s_range.^2) + A4(1)*exp(-B4(6)*s_range.^2) + A5(1)*exp(-B5(6)*s_range.^2))./exp(-DW(6)/2*s_range.^2));
plot(s_range, (A1(7)*exp(-B1(7)*s_range.^2) + A2(1)*exp(-B2(7)*s_range.^2) + A3(1)*exp(-B3(7)*s_range.^2) + A4(1)*exp(-B4(7)*s_range.^2) + A5(1)*exp(-B5(7)*s_range.^2))./exp(-DW(7)/2*s_range.^2));
plot(s_range, (A1(8)*exp(-B1(8)*s_range.^2) + A2(1)*exp(-B2(8)*s_range.^2) + A3(1)*exp(-B3(8)*s_range.^2) + A4(1)*exp(-B4(8)*s_range.^2) + A5(1)*exp(-B5(8)*s_range.^2))./exp(-DW(8)/2*s_range.^2));
plot(s_range, (A1(9)*exp(-B1(9)*s_range.^2) + A2(1)*exp(-B2(9)*s_range.^2) + A3(1)*exp(-B3(9)*s_range.^2) + A4(1)*exp(-B4(9)*s_range.^2) + A5(1)*exp(-B5(9)*s_range.^2))./exp(-DW(9)/2*s_range.^2));
plot(s_range, (A1(10)*exp(-B1(10)*s_range.^2) + A2(1)*exp(-B2(10)*s_range.^2) + A3(1)*exp(-B3(10)*s_range.^2) + A4(1)*exp(-B4(10)*s_range.^2) + A5(1)*exp(-B5(10)*s_range.^2))./exp(-DW(10)/2*s_range.^2));

AgA = [-0.0222, 0.1318, 0.0837, 0.0372, 0.0019];
AgB = [0.0298, 0.1874, 0.6389, 2.4460, 8.6206];
AgDW = 0.1933;

figure;
%plot(s_range,(AgA(1)*exp(-AgB(1)*s_range.^2) + AgA(2)*exp(-AgB(2)*s_range.^2) + AgA(3)*exp(-AgB(3)*s_range.^2) + AgA(4)*exp(-AgB(4)*s_range.^2) + AgA(5)*exp(-AgB(5)*s_range.^2)).*exp(-AgDW/2*s_range.^2));
%hold on;
%plot(s_range,(AgA(1)*exp(-AgB(1)*s_range.^2) + AgA(2)*exp(-AgB(2)*s_range.^2) + AgA(3)*exp(-AgB(3)*s_range.^2) + AgA(4)*exp(-AgB(4)*s_range.^2) + AgA(5)*exp(-AgB(5)*s_range.^2)));
%hold on;
plot(s_range,(AgA(1)*exp(-AgB(1)*s_range.^2) + AgA(2)*exp(-AgB(2)*s_range.^2) + AgA(3)*exp(-AgB(3)*s_range.^2) + AgA(4)*exp(-AgB(4)*s_range.^2) + AgA(5)*exp(-AgB(5)*s_range.^2))*Beta100/Beta200);