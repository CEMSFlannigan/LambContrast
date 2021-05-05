max_pos_xx_index = 0;
max_pos_yy_index = 0;
max_pos_xy_index = 0;

max_thick_xx_index = 0;
max_thick_yy_index = 0;
max_thick_xy_index = 0;

max_strain_xx = 0;
max_strain_yy = 0;
max_strain_xy = 0;

for i = 1:length(strain_matr)
    cur_strain_set = cell2mat(strain_matr(i));
    
    [cur_set_xx cur_set_xx_thick] = max(abs(cur_strain_set(:,1)));
    [cur_set_yy cur_set_yy_thick] = max(abs(cur_strain_set(:,2)));
    [cur_set_xy cur_set_xy_thick] = max(abs(cur_strain_set(:,3)));
    
    if (cur_set_xx > abs(max_strain_xx))
        max_strain_xx = cur_set_xx;
        max_thick_xx_index = cur_set_xx_thick;
        max_pos_xx_index = i;
    end
    
    if (cur_set_yy > abs(max_strain_yy))
        max_strain_yy = cur_set_yy;
        max_thick_yy_index = cur_set_yy_thick;
        max_pos_yy_index = i;
    end
    
    if (cur_set_xy > abs(max_strain_xy))
        max_strain_xy = cur_set_xy;
        max_thick_xy_index = cur_set_xy_thick;
        max_pos_xy_index = i;
    end
end

max_strains = [max_strain_xx, max_strain_yy, max_strain_xy];
max_strain_pos = [max_pos_xx_index, max_pos_yy_index, max_pos_xy_index];
max_strain_thick = [max_thick_xx_index, max_thick_yy_index, max_thick_xy_index];