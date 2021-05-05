% h1 = figure; hold on;
% plot(xvals,yval_arr(:,1),'k', 'LineWidth',8);
% 
% plot(xvals,yval_arr(:,end),'k', 'LineWidth',8);
% 
% 
% y13 = round(length(y_adj(1,:,1))/3);
% y23 = 2*y13;
% x13 = round(length(y_adj(:,1,1))/3);
% x23 = 2*x13;
% 
% plot(xvals,yval_arr(:,y23,1),'k', 'LineWidth',8);
% plot(xvals,yval_arr(:,y13,1),'k', 'LineWidth',8);
% 
% plot(xvals(1).*ones(size(yval_arr(1,:,1))),yval_arr(1,:,1),'k','LineWidth',8);
% plot(xvals(end).*ones(size(yval_arr(end,:,1))),yval_arr(end,:,1),'k','LineWidth',8);
% plot(xvals(x13).*ones(size(yval_arr(x13,:,1))),yval_arr(x13,:,1),'k', 'LineWidth',8);
% plot(xvals(x23).*ones(size(yval_arr(x13,:,1))),yval_arr(x23,:,1),'k', 'LineWidth',8);
% 
% set(gcf,'Position',[0 0 2000 2000]);
% box off; axis off;
% pbaspect([1 1 1]);
% 
% set(gca,'FontSize',20);
% set(gca,'YTickLabel',[]);
% set(gca,'XTickLabel',[]);
% 
% print('Undistorted','-dpng','-r300');
%  
% h2 = figure;
% C = linspecer(5,'qualitative');
% filebase = 'S:\DanielD\Contrast Analysis\Distortion_Movies\Distort_';
% for k = 1:1%length(timevals)
%     plot(x_adj(:,end,k),y_adj(:,end,k),'k', 'LineWidth',10, 'Color',C(3,:));
%     hold on;
%     
%     plot(x_adj(:,y23,k),y_adj(:,y23,k),'k', 'LineWidth',10);
%     plot(x_adj(:,y13,k),y_adj(:,y13,k),'k', 'LineWidth',10);
%     
%     plot(x_adj(1,:,k),y_adj(1,:,k),'k','LineWidth',10);
%     plot(x_adj(x13,:,k),y_adj(x13,:,k),'k', 'LineWidth',10);
%     plot(x_adj(x23,:,k),y_adj(x23,:,k),'k', 'LineWidth',10);
%     plot(x_adj(end,:,k),y_adj(end,:,k),'k','LineWidth',10, 'Color',C(5,:));
%     
%     plot(x_adj(:,1,k),y_adj(:,1,k),'k', 'LineWidth',10);
%     
%     set(gcf,'Position',[0 0 2000 2000]);
%     box off;%box on;
%     axis off;
%     pbaspect([1 1 1]);
%     
%     set(gca,'FontSize',20);
%     set(gca,'YTickLabel',[]);
%     set(gca,'XTickLabel',[]);
%     xlim([-5 2000]);
%     
%     numdig = length(num2str(k));
%     needadd = 4-numdig;
%     add = num2str(10^needadd);
%     add = add(2:end);
%     filename = strcat(filebase,add);
%     filename = strcat(filename,num2str(k));
%     filename = strcat(filename,'.png');
%     
%     saveas(h2,filename);
%     hold off;
% end
% 
% print('Distorted_Actual','-dpng','-r300');

% figure; 
% 
% plot(xvals(x23:end),yval_arr(x23:end,end),'k', 'LineWidth',6);
% 
% set(gcf,'Position',[0 0 2000 2000]);
% box off;%box on;
% axis off;
% pbaspect([1 1 1]);
% 
% set(gca,'FontSize',20);
% set(gca,'YTickLabel',[]);
% set(gca,'XTickLabel',[]);
% 
% print('Undistorted_segment','-dpng','-r300');
% 
% figure;
% 
% plot(x_adj(x23:end,end,k),y_adj(x23:end,end,k),'k', 'LineWidth',6);
% 
% set(gcf,'Position',[0 0 2000 2000]);
% box off;%box on;
% axis off;
% pbaspect([1 1 1]);
% 
% set(gca,'FontSize',20);
% set(gca,'YTickLabel',[]);
% set(gca,'XTickLabel',[]);
% 
% print('Distorted_segment','-dpng','-r300');

% figure;
% for k = 1:length(timevals)
%     hold on;
%     plot(xvals,x_res(:,end,k));
%     pause(0.05);
% end
% 
% figure;
% for k = 1:length(timevals)
%     hold on;
%     plot(xvals,y_res(:,end,k));
%     pause(0.05);
% end

figure;
% for k = 1:length(timevals)
%     subplot(1,2,1); hold on;
hold on;
    plot(xvalsm,x_strain(:,end,k)*100,'b','LineWidth',6);
    plot(xvals,y_strain(:,end,k)*100,'r','LineWidth',6);
    pbaspect([1 1 1]);
    box on;
    set(gca,'FontSize',48,'LineWidth',10,'FontWeight','bold');
    set(gca,'TickDir','out');
    xticks([0 250 500 750 1000]);
    xlabel("Position x-dir (nm)",'FontSize',50,'FontWeight','bold');
    ylabel("Strain (%)","FontSize",50,'FontWeight','bold');
    set(gcf,'Position',[-25 -25 2000 2000]);
    %legend({" x-strain"," y-strain"},'FontSize',36,'Location','se','box','off');
    print('Strain_x_dir','-dpng','-r300');
%     subplot(1,2,2); hold on;
%     
%     pause(0.05);
%end

figure;
hold on;
% for k = 1:length(timevals)
%     subplot(1,2,1); hold on;
%     plot(yval_arr(end,:),x_strain(end,:,k));
%     subplot(1,2,2); hold on;
    plot(yval_arr(end,:),x_strain(end,:,k)*100,'b','LineWidth',6);
    plot(yval_arrm(end,:),y_strain(end,:,k)*100,'r','LineWidth',6);
    pbaspect([1 1 1]);
    box on;
    set(gca,'FontSize',48,'LineWidth',10);
    set(gca,'TickDir','out');
    xticks([0 250 500 750 1000]);
    ylim([-7 6]);
    xlabel("Position y-dir (nm)",'FontSize',50,'FontWeight','bold');
    ylabel("Strain (%)","FontSize",50,'FontWeight','bold');
    set(gcf,'Position',[-25 -25 2000 2000]);
    %legend({" x-strain"," y-strain"},'FontSize',36,'Location','se','box','off');   
    print('Strain_y_dir','-dpng','-r300');
    pause(0.05);
%end