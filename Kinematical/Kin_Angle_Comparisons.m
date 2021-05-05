Alph_Tilts = [0,0,0,0,0,0,0,2,5,5,5,5,5,5,5,10,10,10,10,10,10,10,20,20,20,20,20,20,20];
Bet_Tilts = [0,2,4,7,12,-3,-8,0,0,2,4,7,12,-3,-8,0,2,4,7,12,-3,-8,0,2,-4,7,12,-3,-8];
Max_Oscillations = 100.*[0.0051,0.0031,0.0154,0.0044,0.0040,0.0029,0.0037,0.0063,0.0106,0.0162,0.0064,0.0033,0.0039,0.0052,0.0022,0.0065,0.0069,0.0038,0.0116,0.0035,0.0049,0.0042,0.0039,0.0054,0.0034,0.0034,0.0033,0.0065,0.0034];

figure;
surf(Alph_Tilts_Interp,Bet_Tilts_Interp,medfilt2(Max_Osc_Interp,[25 25]),'edgecolor','none');
shading interp

xlabel('Alpha Tilt');
ylabel('Beta Tilt');
zlabel('Intensity Change (%)');

figure;
surf(Alph_Tilts_Interp,Bet_Tilts_Interp,imfilter(Max_Osc_Interp,fspecial('average',[20 20])),'edgecolor','none');
colorMap = gray(512);
colormap(colorMap);
shading interp

pbaspect([1 1 1]);
xlim([0 20]);
ylim([-8 12]);
set(gcf,'Position',[25 25 2500 2500]);
set(gca,'FontSize', 48);
set(gca,'LineWidth',10);
set(gca,'TickDir','out');
box on;
colorbar('AxisLocation','out','LineWidth',5,'Position',[0.591493054996762,0.360078277886497,0.011111111111111,0.565557729941291]);
xlabel('Alpha Tilt (degrees)','FontSize',50);
ylabel('Beta Tilt (degrees)','FontSize',50);
zlabel('Intensity Change (%)');
view(0,90);