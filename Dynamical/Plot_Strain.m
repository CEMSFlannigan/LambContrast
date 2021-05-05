numFrames = 17;
x_strain_time = zeros(1,numFrames);
y_strain_time = zeros(1,numFrames);
xy_strain_time = zeros(1,numFrames);
curmaxx = 0;
curmaxy = 0;
curmaxxy = 0;
for k = 1:numFrames%length(timevals)
    curstrain = strain_matr{111,k};
    
    [Mx,I] = max(abs(curstrain(:,1)));
    x_strain_time(k) = curstrain(I,1);
    [My,I] = max(abs(curstrain(:,2)));
    y_strain_time(k) = curstrain(I,2);
    [Mxy,I] = max(abs(curstrain(:,3)));
    xy_strain_time(k) = curstrain(I,3);
    if Mx > curmaxx
        curmaxx = Mx;
    end
    
    if My > curmaxy
        curmaxy = My;
    end
    
    if Mxy > curmaxxy
        curmaxxy = Mxy;
    end
end

curmaxx
curmaxy
curmaxxy

figure;
subplot(1,3,1);
plot(x_strain_time,'b');
subplot(1,3,2);
plot(y_strain_time,'r');
subplot(1,3,3);
plot(xy_strain_time,'k');
% 
% xdil = zeros(size(xvals));
% ydil = zeros(size(xvals));
% xysh = zeros(size(xvals));
% for i = 1:length(xvals)
%     curstrain = strain_matr{i,1};
%     xdil(i) = curstrain(1,1);
%     ydil(i) = curstrain(1,2);
%     xysh(i) = curstrain(1,3);
% end
% figure;
% hold on;
% plot(xvals,xdil*100,'b','LineWidth',9);
% plot(xvals,ydil*100,'r','LineWidth',9);
% plot(xvals,xysh*100,'k','LineWidth',9);
% pbaspect([1 1 1]);
% box on;
% set(gca,'FontSize',55,'LineWidth',12,'FontWeight','bold');
% set(gca,'TickDir','out');
% xlim([0 1000]);
% ylim([-0.6 0.6]);
% yticks([-0.6 -0.3 0 0.3 0.6]);
% xticks([0 250 500 750 1000]);
% xlabel("Position {\it x}-dir (nm)",'FontSize',65,'FontWeight','bold');
% ylabel("Strain (%)","FontSize",65,'FontWeight','bold');
% set(gcf,'Position',[-25 -25 2000 2000]);
% %legend({" {\it x}-dilatational"," {\it y}-dilatational"," {\it xy}-shear"},'FontSize',30,'Location','se','box','off');
% % print('Strain_x_dir','-dpng','-r300');
% 
% figure;
% 
% curstrain = strain_matr{end,1};
% xdil = zeros(1,length(curstrain));
% ydil = zeros(1,length(curstrain));
% xysh = zeros(1,length(curstrain));
% for i = 1:length(curstrain)
%     xdil(i) = curstrain(i,1);
%     ydil(i) = curstrain(i,2);
%     xysh(i) = curstrain(i,3);
% end
% yvals = linspace(0,100,length(curstrain));
% hold on;
% plot(yvals-50,xdil*100,'b','LineWidth',9);
% plot(yvals-50,ydil*100,'r','LineWidth',9);
% plot(yvals-50,xysh*100,'k','LineWidth',9);
% pbaspect([1 1 1]);
% box on;
% set(gca,'FontSize',55,'LineWidth',12,'FontWeight','bold');
% set(gca,'TickDir','out');
% xticks([-50 -25 0 25 50]);
% xlim([-50 50]);
% ylim([-0.25 0.5]);
% yticks([-0.25 0 0.25 0.5]);
% xlabel("Position {\it y}-dir (nm)",'FontSize',65,'FontWeight','bold');
% ylabel("Strain (%)","FontSize",65,'FontWeight','bold');
% set(gcf,'Position',[-25 -25 2000 2000]);
%legend({" {\it x}-dilatational"," {\it y}-dilatational"," {\it xy}-shear"},'FontSize',30,'Location','se','box','off');
% print('Strain_y_dir','-dpng','-r300');
% pause(0.05);

% hold on;
% for k = 1:length(timevals)
%     curstrain = strain_matr{round(end/2),k};
%     
%     plot(curstrain(:,1));
% end
% 
% subplot(1,3,2);
% hold on;
% for k = 1:length(timevals)
%     curstrain = strain_matr{round(end),k};
%     
%     plot(curstrain(:,2));
% end
% 
% subplot(1,3,3);
% hold on;
% for k = 1:length(timevals)
%     curstrain = strain_matr{round(end),k};
%     
%     plot(curstrain(:,3));
% end

% k0p0126_LA3_str_x = [0.00550475732675804];
% k0p0126_LA3_str_y = [0.00733251194437042];
% k0p0126_LA3_str_z = [0.00555172328352643];
% 
% k0p0126_LA6_str_x = [0.0110095146535161];
% k0p0126_LA6_str_y = [0.0146650238887408];
% k0p0126_LA6_str_z = [0.0111034465670529];
% 
% k0p0126_LA12_str_x = [0.0220190293070321];
% k0p0126_LA12_str_y = [0.0293300477774817];
% k0p0126_LA12_str_z = [0.0222068931341057];
% 
% k0p0273_LA3_str_x = [0.00367277111025167];
% k0p0273_LA3_str_y = [0.00480612480044810];
% k0p0273_LA3_str_z = [0.00373317570869010];
% 
% k0p0273_LA6_str_x = [0.00734554222050335];
% k0p0273_LA6_str_y = [0.00961224960089620];
% k0p0273_LA6_str_z = [0.00746635141738019];
% 
% k0p0273_LA12_str_x = [0.0146910844410067];
% k0p0273_LA12_str_y = [0.0192244992017924];
% k0p0273_LA12_str_z = [0.0149327028347604];
% 
% k0p00497_LA3_str_x = [0.00795364804816097];
% k0p00497_LA3_str_y = [0.0107045436021833];
% k0p00497_LA3_str_z = [0.00799608426562321];
% 
% k0p00497_LA6_str_x = [0.0159072960963219];
% k0p00497_LA6_str_y = [0.0214090872043667];
% k0p00497_LA6_str_z = [0.0159921685312464];
% 
% k0p00497_LA12_str_x = [0.0318145921926439];
% k0p00497_LA12_str_y = [0.0428181744087333];
% k0p00497_LA12_str_z = [0.0319843370624928];
% 
% figure;
% k_vals = [0.00497, 0.0126, 0.0273];
% LA_vals = [3, 6, 12];
% scatter(k_vals,100.*[k0p00497_LA3_str_x, k0p0126_LA3_str_x, k0p0273_LA3_str_x],100,'k','MarkerFaceColor','k');
% hold on;
% scatter(k_vals,100.*[k0p00497_LA6_str_x, k0p0126_LA6_str_x, k0p0273_LA6_str_x],100,'r','MarkerFaceColor','r');
% scatter(k_vals,100.*[k0p00497_LA12_str_x, k0p0126_LA12_str_x, k0p0273_LA12_str_x],100,'b','MarkerFaceColor','b');
% pbaspect([1 1 1]);
% set(gcf,'Position',[25 25 2500 2500]);
% set(gca,'FontSize', 48);
% set(gca,'LineWidth',5);
% set(gca,'TickDir','out');
% box on;
% xlabel('Wavenumber (nm^{-1})', 'FontSize', 50);
% ylabel('Strain (%)','FontSize', 50);
% legend({'3 nm^2', '6 nm^2', '12 nm^2'},'FontSize',35,'box','off');

% k_arr = [0.0126, 0.0126, 0.0126, 0.0273, 0.0273, 0.0273, 0.00497, 0.00497, 0.00497];
% LA_arr = [3, 6, 12, 3, 6, 12, 3, 6, 12];
% str_x_arr = [k0p0126_LA3_str_x, k0p0126_LA6_str_x, k0p0126_LA12_str_x, k0p0273_LA3_str_x, k0p0273_LA6_str_x, k0p0273_LA12_str_x, k0p00497_LA3_str_x, k0p00497_LA6_str_x, k0p00497_LA12_str_x];
% str_y_arr = [k0p0126_LA3_str_y, k0p0126_LA6_str_y, k0p0126_LA12_str_y, k0p0273_LA3_str_y, k0p0273_LA6_str_y, k0p0273_LA12_str_y, k0p00497_LA3_str_y, k0p00497_LA6_str_y, k0p00497_LA12_str_y];
% str_z_arr = [k0p0126_LA3_str_z, k0p0126_LA6_str_z, k0p0126_LA12_str_z, k0p0273_LA3_str_z, k0p0273_LA6_str_z, k0p0273_LA12_str_z, k0p00497_LA3_str_z, k0p00497_LA6_str_z, k0p00497_LA12_str_z];
% figure;
% hold on;
% scatter3(k_arr,LA_arr,str_x_arr,'k');
% scatter3(k_arr,LA_arr,str_y_arr,'b');
% scatter3(k_arr,LA_arr,str_z_arr,'r');