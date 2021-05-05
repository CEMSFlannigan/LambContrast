% h2 = figure;
% 
% mkdir('S:\DanielD\Contrast Analysis\Kinematic\Movie_400beam_5A_0G_2B_Clean_no_vac_colorbar');
% filebase = 'S:\DanielD\Contrast Analysis\Kinematic\Movie_400beam_5A_0G_2B_Clean_no_vac_colorbar\Distort_';
% for k = 1:length(timevals)
% 
%     imagesc(ex_image{k});
%     colorMap = gray(512);
%     colormap(colorMap);
%     %colorbar;
%     
%     set(gcf,'Position', [0 0 1500 1000]);
%     set(gca,'FontSize', 50, 'LineWidth', 1.5);
%     pbaspect([1 1 1]);
%     %xlabel('x direction (nm)', 'FontSize', 20);
%     %ylabel('z direction (nm)', 'FontSize', 20);
%     ticks = ((1000*linspace(min_bf_contrast,max_bf_contrast,5)) - mod((1000*linspace(min_bf_contrast,max_bf_contrast,5)),1) + round(mod((1000*linspace(min_bf_contrast,max_bf_contrast,5)),1)))/1000;
%     colorbar('Ticks',ticks,'LineWidth',5,'TickDirection','out','Position',[0.8, .118, .03, .8]);
%     yticks([]);
%     xticks([]);
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

% figure;
% for k = 1:length(timevals)
%     subplot(2,2,1);hold on;
%     plot(xvals,s_g_arr_1(1,:,k));
%     plot(xvals,s_g_arr_2(1,:,k));
%     hold off;
%     
%     subplot(2,2,2);hold on;
%     plot(xvals,D_arr_1(1,:,k));
%     plot(xvals,D_arr_2(1,:,k));
%     hold off;
%     
%     subplot(2,2,3);hold on;
%     plot(xvals,sum_D_arr_sq(1,:,k));
%     hold off;
% end
% 
% figure;
% for k = 1:length(timevals)
%     plot(xvals, bf_contrast(:,k));%,'k','LineWidth',2);
%     hold on;
% end
% set(gcf,'Position', [0 0 1000 1000]);
% set(gca,'FontSize', 20);
% box on;
% xlabel('Distance (nm)', 'FontSize', 20);
% ylabel('Normalized Intensity', 'FontSize', 20);

% figure;
% C = linspecer(7,'sequential');
% h1 = plot(xvals(1:look_length),bf_contrast(1:look_length,1)-0.04,'LineWidth',6,'Color',C(1,:));
% hold on; 
% h2 = plot(xvals(1:look_length),bf_contrast(1:look_length,07) + 0.01,'LineWidth',6,'Color',C(2,:)); % 0.06
% h3 = plot(xvals(1:look_length),bf_contrast(1:look_length,13) + 0.06,'LineWidth',6,'Color',C(3,:)); % 0.06
% h4 = plot(xvals(1:look_length),bf_contrast(1:look_length,18) + 0.11,'LineWidth',6,'Color',C(4,:)); % 0.06
% h5 = plot(xvals(1:look_length),bf_contrast(1:look_length,24) + 0.16,'LineWidth',6,'Color',C(5,:)); % 0.12
% h6 = plot(xvals(1:look_length),bf_contrast(1:look_length,30) + 0.21,'LineWidth',6,'Color',C(6,:)); % 0.12
% pbaspect([1 1 1]);
% xlim([0 1000]);
% %ylim([0.1 0.6]);
% %legend([h1, h3, h5, h7], {" t_i", " t_i + 10 ps", " t_i + 20 ps", " t_i + 30 ps"},'FontSize',25,'Location','ne','Box','off');
% set(gcf,'Position',[25 25 2500 2500]);
% set(gca,'FontSize', 48);
% set(gca,'LineWidth',5);
% set(gca,'TickDir','out');
% yticks([]);
% xticks([0 250 500 750 1000]);
% xlabel('Position (nm)', 'FontSize', 50);
% ylabel('Intensity','FontSize', 50);
% print('BF_Contrast_5A_0G_2B_Plot','-dpng','-r300');
% 
