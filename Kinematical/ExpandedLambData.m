tic;
k_mag = 398; % nm^-1
alpha_tilt = 5;%-15.22;%-4.33 mimics UEM;
gam_tilt = 0;
beta_tilt = 2;%-3.64;
k_unit = [0, -1, 0];
k_unit = k_unit/sqrt(dot(k_unit,k_unit));
k_alph_rot = [1,0,0;0,cosd(alpha_tilt),-sind(alpha_tilt);0,sind(alpha_tilt),cosd(alpha_tilt)];
k_gam_rot = [cosd(gam_tilt),0,sind(gam_tilt);0,1,0;-sind(gam_tilt),0,cosd(gam_tilt)];
k_beta_rot = [cosd(beta_tilt),-sind(beta_tilt),0;sind(beta_tilt),cosd(beta_tilt),0;0,0,1];
k_unit = k_beta_rot*k_gam_rot*k_alph_rot*k_unit';
k_actual = k_unit*k_mag;

%%
% Calculate thicknesses and strain
w0 = waitbar(0,'Generating additional data.');

strain_sample = floor(yvals(2)*2/latpar);
xvals_exp = -1*xvals(end):x_spacing:xvals(end)*2;

yval_exp_arr = zeros(length(xvals),y_point_num);
for i = 1:length(xvals_exp)
    cury = slope*(xvals_exp(i))+yvals(1);
    yval_exp_arr(i,:) = linspace(-1*cury,cury,y_point_num);
end

waitbar(1/3);

x_res_exp = zeros(length(xvals_exp),y_point_num,length(timevals));
y_res_exp = zeros(length(xvals_exp),y_point_num,length(timevals));
x_strain_exp = zeros(length(xvals_exp)-1,y_point_num,length(timevals));
y_strain_exp = zeros(length(xvals_exp),y_point_num-1,length(timevals));
x_adj_exp = zeros(size(x_res_exp));
y_adj_exp = zeros(size(y_res_exp));

for k = 1:length(timevals)
    parfor j = 1:y_point_num
        curt = timevals(k);
        x_res_exp(:,j,k) = double(combined_x(xvals_exp,yval_exp_arr(:,j)',curt));
        y_res_exp(:,j,k) = double(combined_y(xvals_exp,yval_exp_arr(:,j)',curt));
        
        x_adj_exp(:,j,k) = xvals_exp' + x_res_exp(:,j,k);
        y_adj_exp(:,j,k) = yval_exp_arr(:,j) + y_res_exp(:,j,k);
    end
end

waitbar(2/3);

for i = 1:length(xvals_exp)
    for j = 1:y_point_num - 1
        yval_arrm_exp(i,j) = (yval_exp_arr(i,j+1) + yval_exp_arr(i,j))/2;
    end
end

for k = 1:length(timevals)
    for i = 1:length(xvals_exp)-1
        for j = 1:y_point_num
            x_adj_dif = (x_adj_exp(i+1,j,k) - x_adj_exp(i,j,k));
            xorig_dif = (xvals_exp(i+1) - xvals_exp(i));
            x_strain_exp(i,j,k) = (x_adj_dif - xorig_dif)./(xorig_dif);
        end
    end
    
    for i = 1:length(xvals_exp)
        for j = 1:y_point_num-1
            y_adj_dif = (y_adj_exp(i,j+1,k) - y_adj_exp(i,j,k));
            yorig_dif = (yval_exp_arr(i,j+1) - yval_exp_arr(i,j));
            y_strain_exp(i,j,k) = (y_adj_dif - yorig_dif)./(yorig_dif);
        end
    end
end

waitbar(3/3);
close(w0);
%%

xvals_top = x_adj_exp(:,end,1);
xvals_bot = x_adj_exp(:,1,1);

yvals_top = y_adj_exp(:,end,1);
yvals_bot = y_adj_exp(:,1,1);

top_inv_slopes = zeros(length(xvals),length(timevals));
bot_inv_slopes = zeros(length(xvals),length(timevals));
top_surf_thickness = zeros(length(xvals),length(timevals));
bot_surf_thickness = zeros(length(xvals),length(timevals));
prop_thickness = zeros(length(xvals),length(timevals));
ave_prop_strain_x = zeros(length(xvals),length(timevals));
ave_prop_strain_y = zeros(length(xvals),length(timevals));
ave_prop_strain = zeros(length(xvals),length(timevals));

w4 = waitbar(0,'Processing time.');
for k = 1:length(timevals)
    w1 = waitbar(0,'Processing strain and thicknesses.');
    curx_bot = x_adj_exp(:,1,k);
    cury_bot = y_adj_exp(:,1,k);
    
    curx_top = x_adj_exp(:,end,k);
    cury_top = y_adj_exp(:,end,k);
 
    %%
    for i = 1:length(xvals)
        ai = i+length(xvals);
        top_slope = -1/((cury_top(ai + 1) - cury_top(ai-1))./(curx_top(ai+1) - curx_top(ai-1)));
        top_inv_slopes(i,k) = top_slope;
        bot_slope = -1/((cury_bot(ai + 1) - cury_bot(ai-1))./(curx_bot(ai+1) - curx_bot(ai-1)));
        bot_inv_slopes(i,k) = bot_slope;
        top_intersect = cury_top(ai) - top_slope*curx_top(ai);
        bot_intersect = cury_bot(ai) - bot_slope*curx_bot(ai);
        
        top_line_y = top_slope.*curx_bot+top_intersect;
        bot_line_y = bot_slope.*curx_top+bot_intersect;
        
        % Find intersect on top (bottom thickness)
        dif_top = bot_line_y - cury_top;
        top_indx = find(diff(sign(dif_top)))+1;
        bot_line = @(x) top_slope*x+top_intersect;
        
        top_curve_slope = (cury_top(top_indx) - cury_top(top_indx-1))./(curx_top(top_indx) - curx_top(top_indx-1));
        top_curve_intersect = cury_top(top_indx) - top_curve_slope*curx_top(top_indx);
        top_curve = @(x) top_curve_slope*x+top_curve_intersect;
        
        top_lc_inter_x = fzero(@(x) top_curve(x) - bot_line(x),curx_top(top_indx));
        top_lc_inter_y = top_curve(top_lc_inter_x);
        
        bot_surf_thickness(i,k) = sqrt((curx_bot(ai) - top_lc_inter_x)^2 + (cury_bot(ai) - top_lc_inter_y)^2);
        
        % Find intersect on bot (top thickness)
        dif_bot = top_line_y - cury_bot;
        bot_indx = find(diff(sign(dif_bot)))+1;
        top_line = @(x) top_slope*x+top_intersect;
        
        bot_curve_slope = (cury_bot(bot_indx) - cury_bot(bot_indx-1))./(curx_bot(bot_indx) - curx_bot(bot_indx-1));
        bot_curve_intersect = cury_bot(bot_indx) - bot_curve_slope*curx_bot(bot_indx);
        bot_curve = @(x) bot_curve_slope*x+bot_curve_intersect;
        
        bot_lc_inter_x = fzero(@(x) bot_curve(x) - top_line(x),curx_bot(bot_indx));
        bot_lc_inter_y = bot_curve(bot_lc_inter_x);
        
        top_surf_thickness(i,k) = sqrt((curx_top(ai) - bot_lc_inter_x)^2 + (cury_top(ai) - bot_lc_inter_y)^2);
    end
    %%
    
    for i = 1:length(xvals)
        ai = i+length(xvals);
        
        prop_slope = k_unit(2)/k_unit(1);
        
        
        
        if prop_slope > 0 % intersection with the top
            prop_intersect = cury_bot(ai) - prop_slope*curx_bot(ai);
            xy_line = @(x) prop_slope*x + prop_intersect;
            
            prop_y = prop_slope*curx_top + prop_intersect;
            
            prop_dif = prop_y - cury_top;
            prop_indx = find(diff(sign(prop_dif)))+1;
            
            prop_curve_slope = (cury_top(prop_indx) - cury_top(prop_indx-1))./(curx_top(prop_indx) - curx_top(prop_indx-1));
            prop_curve_intersect = cury_top(prop_indx) - prop_curve_slope*curx_top(prop_indx);
            
            prop_curve = @(x) prop_curve_slope*x+prop_curve_intersect;
            
            prop_lc_inter_x = fzero(@(x) prop_curve(x) - xy_line(x),curx_top(prop_indx));
            prop_lc_inter_y = prop_curve(prop_lc_inter_x);
            
            lateral_dist = (prop_lc_inter_y - cury_bot(ai))/k_unit(2)*k_unit(3);
            prop_thickness(i,k) = sqrt(lateral_dist^2 + (prop_lc_inter_x - curx_bot(ai))^2 + (prop_lc_inter_y - cury_bot(ai))^2);
            
            xythickness = sqrt((prop_lc_inter_x - curx_bot(ai))^2 + (prop_lc_inter_y - cury_bot(ai))^2);
            
            line_points_y = linspace(prop_lc_inter_y,cury_bot(ai),strain_sample);
            line_points_x = linspace(prop_lc_inter_x,curx_bot(ai),strain_sample);
            
            tot_strain_x = 0;
            tot_strain_y = 0;
            for curidx = 1:length(line_points_y)
                idxy1 = 0;
                idxy2 = 0;
                idxx1 = 0;
                idxx2 = 0;
                strain1 = 0; % top left
                strain2 = 0; % top right
                strain3 = 0; % bottom left
                strain4 = 0; % bottom right
                cur_str_y = line_points_y(curidx);
                cur_str_x = line_points_x(curidx);
                [strain_indexes closest_dif_y closest_dif_x] = findClosest(x_adj_exp(:,:,k),y_adj_exp(:,:,k),cur_str_x,cur_str_y);
                
                if closest_dif_y > 0 && strain_indexes(2) > 1 && strain_indexes(2) < size(y_strain_exp,2)-1
                    idxy1 = strain_indexes(2)-1;
                    idxy2 = strain_indexes(2);
                elseif closest_dif_y < 0 && strain_indexes(2) > 1 && strain_indexes(2) < size(y_strain_exp,2)-1
                    idxy1 = strain_indexes(2);
                    idxy2 = strain_indexes(2)+1;
                elseif closest_dif_y == 0 && strain_indexes(2) > 1 && strain_indexes(2) < size(y_strain_exp,2)-1
                    idxy1 = strain_indexes(2);
                    idxy2 = strain_indexes(2);
                elseif strain_indexes(2) <= 1
                    idxy1 = 1;
                    idxy2 = 1;
                elseif strain_indexes(2) >= size(y_strain_exp,2)-1
                    idxy1 = size(y_strain_exp,2);
                    idxy2 = size(y_strain_exp,2);
                end
                
                if closest_dif_x > 0 && strain_indexes(1) > 1 && strain_indexes(1) < size(x_strain_exp,1)
                    idxx1 = strain_indexes(1)-1;
                    idxx2 = strain_indexes(1);
                elseif closest_dif_x < 0 && strain_indexes(1) > 1 && strain_indexes(1) < size(x_strain_exp,1)
                    idxx1 = strain_indexes(1);
                    idxx2 = strain_indexes(1)+1;
                elseif closest_dif_x == 0 && strain_indexes(1) > 1 && strain_indexes(1) < size(x_strain_exp,1)
                    idxx1 = strain_indexes(1);
                    idxx2 = strain_indexes(1);
                elseif strain_indexes(1) <= 1
                    idxx1 = 1;
                    idxx2 = 1;
                elseif strain_indexes(1) >= size(x_strain_exp,1)
                    idxx1 = size(x_strain_exp,1);
                    idxx2 = size(x_strain_exp,1);
                end
                
                tot_strain_x = tot_strain_x + (x_strain_exp(idxx1,idxy1,k) + x_strain_exp(idxx1,idxy2,k) + x_strain_exp(idxx2,idxy1,k) + x_strain_exp(idxx2,idxy2,k))/4;
                tot_strain_y = tot_strain_y + (y_strain_exp(idxx1,idxy1,k) + y_strain_exp(idxx1,idxy2,k) + y_strain_exp(idxx2,idxy1,k) + y_strain_exp(idxx2,idxy2,k))/4;
                
            end
            ave_prop_strain_x(i,k) = tot_strain_x/strain_sample;
            ave_prop_strain_y(i,k) = tot_strain_y/strain_sample;
            
            ave_prop_strain(i,k) = sqrt(ave_prop_strain_x(i,k)^2 + ave_prop_strain_y(i,k)^2);
            
        elseif prop_slope < 0 % intersection with the bottom
            
            %% NEED TO THROW IN EDGE CASE FOR INFINITE PROPAGATION SLOPE
            if ~isinf(prop_slope)
                prop_intersect = cury_top(ai) - prop_slope*curx_top(ai);
                xy_line = @(x) prop_slope*x + prop_intersect;
                
                prop_y = prop_slope*curx_bot + prop_intersect;
                
                prop_dif = prop_y - cury_bot;
                prop_indx = find(diff(sign(prop_dif)))+1;
                
                prop_curve_slope = (cury_bot(prop_indx) - cury_bot(prop_indx-1))./(curx_bot(prop_indx) - curx_bot(prop_indx-1));
                prop_curve_intersect = cury_bot(prop_indx) - prop_curve_slope*curx_bot(prop_indx);
                
                prop_curve = @(x) prop_curve_slope*x+prop_curve_intersect;
                
                prop_lc_inter_x = fzero(@(x) prop_curve(x) - xy_line(x),curx_bot(prop_indx));
                prop_lc_inter_y = prop_curve(prop_lc_inter_x);
            else
                xy_line = @(x) cury_top(ai);
                
                prop_lc_inter_x = curx_top(ai);
                prop_lc_inter_y = cury_bot(ai);
            end
            
            lateral_dist = (prop_lc_inter_y - cury_top(ai))/k_unit(2)*k_unit(3);
            prop_thickness(i,k) = sqrt(lateral_dist^2 + (prop_lc_inter_x - curx_top(ai))^2 + (prop_lc_inter_y - cury_top(ai))^2);
            
            xythickness = sqrt((prop_lc_inter_x - curx_top(ai))^2 + (prop_lc_inter_y - cury_top(ai))^2);
            
            line_points_y = linspace(prop_lc_inter_y,cury_top(ai),strain_sample);
            line_points_x = linspace(prop_lc_inter_x,curx_top(ai),strain_sample);
            
            tot_strain_x = 0;
            tot_strain_y = 0;
            
            for curidx = 1:length(line_points_y)
                idxy1 = 0;
                idxy2 = 0;
                idxx1 = 0;
                idxx2 = 0;
                strain1 = 0; % top left
                strain2 = 0; % top right
                strain3 = 0; % bottom left
                strain4 = 0; % bottom right
                cur_str_y = line_points_y(curidx);
                cur_str_x = line_points_x(curidx);

                [strain_indexes closest_dif_y closest_dif_x] = findClosest(x_adj_exp(:,:,k),y_adj_exp(:,:,k),cur_str_x,cur_str_y);
                
                if closest_dif_y > 0 && strain_indexes(2) > 1 && strain_indexes(2) < size(y_strain_exp,2)-1
                    idxy1 = strain_indexes(2)-1;
                    idxy2 = strain_indexes(2);
                elseif closest_dif_y < 0 && strain_indexes(2) > 1 && strain_indexes(2) < size(y_strain_exp,2)-1
                    idxy1 = strain_indexes(2);
                    idxy2 = strain_indexes(2)+1;
                elseif closest_dif_y == 0 && strain_indexes(2) > 1 && strain_indexes(2) < size(y_strain_exp,2)-1
                    idxy1 = strain_indexes(2);
                    idxy2 = strain_indexes(2);
                elseif strain_indexes(2) <= 1
                    idxy1 = 1;
                    idxy2 = 1;
                elseif strain_indexes(2) >= size(y_strain_exp,2)-1
                    idxy1 = size(y_strain_exp,2);
                    idxy2 = size(y_strain_exp,2);
                end
                
                if closest_dif_x > 0 && strain_indexes(1) > 1 && strain_indexes(1) < size(x_strain_exp,1)
                    idxx1 = strain_indexes(1)-1;
                    idxx2 = strain_indexes(1);
                elseif closest_dif_x < 0 && strain_indexes(1) > 1 && strain_indexes(1) < size(x_strain_exp,1)
                    idxx1 = strain_indexes(1);
                    idxx2 = strain_indexes(1)+1;
                elseif closest_dif_x == 0 && strain_indexes(1) > 1 && strain_indexes(1) < size(x_strain_exp,1)
                    idxx1 = strain_indexes(1);
                    idxx2 = strain_indexes(1);
                elseif strain_indexes(1) <= 1
                    idxx1 = 1;
                    idxx2 = 1;
                elseif strain_indexes(1) >= size(x_strain_exp,1)
                    idxx1 = size(x_strain_exp,1);
                    idxx2 = size(x_strain_exp,1);
                end
                
                tot_strain_x = tot_strain_x + (x_strain_exp(idxx1,idxy1,k) + x_strain_exp(idxx1,idxy2,k) + x_strain_exp(idxx2,idxy1,k) + x_strain_exp(idxx2,idxy2,k))/4;
                tot_strain_y = tot_strain_y + (y_strain_exp(idxx1,idxy1,k) + y_strain_exp(idxx1,idxy2,k) + y_strain_exp(idxx2,idxy1,k) + y_strain_exp(idxx2,idxy2,k))/4;
                
            end
            ave_prop_strain_x(i,k) = tot_strain_x/strain_sample;
            ave_prop_strain_y(i,k) = tot_strain_y/strain_sample;
            
            ave_prop_strain(i,k) = sqrt(ave_prop_strain_x(i,k)^2 + ave_prop_strain_y(i,k)^2);
            
        end
        waitbar(i/length(xvals));
    end
    close(w1);
    waitbar(k/length(timevals));
end
close(w4);

toc