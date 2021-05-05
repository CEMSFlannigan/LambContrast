function [gvec_list,hvec_list] = RecipRasterGenerator(rec_angle_limit,LZ_Limit,g2h_matr,resolution, kvec)

kvec_norm = kvec'./sqrt(sum(kvec.*kvec));

gvec_1 = [1, 0, (1-kvec_norm(1))/kvec_norm(3)] - kvec_norm;
gvec_1 = gvec_1./sqrt(sum(gvec_1.*gvec_1));
gvec_2 = [0, 1, (1-kvec_norm(2))/kvec_norm(3)] - kvec_norm;
gvec_2 = gvec_2./sqrt(sum(gvec_2.*gvec_2));

tot = 0;
if (mod(rec_angle_limit,resolution) ~= 0)
    gvec_list = [0,0,0];
    tot = tot+1;
end

x_range = -abs(rec_angle_limit):resolution:abs(rec_angle_limit);
y_range = -abs(rec_angle_limit):resolution:abs(rec_angle_limit);
z_range = -LZ_Limit/2:resolution:abs(LZ_Limit);

gvec_list = zeros(length(x_range)*length(y_range)*length(z_range),3);

disp(length(x_range));
for i = 1:length(x_range)
    disp(i);
    for j = 1:length(y_range)
        for k = 1:length(z_range)
        	tot = tot + 1;
            curadd_plane = x_range(i).*gvec_1 + y_range(j).*gvec_2;
            curadd_HOLZ =  - z_range(k).*kvec_norm;
            if (sqrt(sum(curadd_plane.^2)) <= rec_angle_limit)
                gvec_list(tot,:) = curadd_plane+curadd_HOLZ;
            end
        end
    end
end

tot = 0;
foundzero = 0;
for i = 1:length(gvec_list)
    curgvec = gvec_list(i,:);
    if (sqrt(sum(curgvec.^2)) == 0 && foundzero == 0)
        foundzero = 1;
        tot = tot + 1;
    elseif (sqrt(sum(curgvec.^2)) == 0 && foundzero == 1)
    else
        tot = tot + 1;
    end
end

truegvec_list = zeros(tot,3);

for i = 1:tot
    curgvec = gvec_list(i,:);
    truegvec_list(i,:) = curgvec;
end
gvec_list = truegvec_list;

hvec_list = zeros(size(gvec_list));

for i = 1:length(gvec_list)
    curgvec = gvec_list(i,:);
    hvec_list(i,:) = (g2h_matr*curgvec')';
end

end