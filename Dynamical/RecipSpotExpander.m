function [gvec_list,hvec_list,tot] = RecipSpotExpander(gvec_list,hvec_list,g2h_matr, h2g_matr, g_lim)

gvec_arr = zeros(length(gvec_list) + 1,3);
hvec_arr = zeros(length(gvec_list) + 1,3);

tot = 1;
for i = 1:length(gvec_list)
    curgvec = gvec_list{i};
    curhvec = hvec_list{i};
    if sqrt(sum(curgvec.^2)) <= g_lim
        gvec_arr(tot + 1,:) = curgvec;
        hvec_arr(tot + 1,:) = curhvec;
        tot = tot + 1;
    end
end

numadd = 13;
newgvec_list = zeros(tot*numadd,3);
newhvec_list = zeros(tot*numadd,3);

dist = 0.08; % 1/nm
%dist in each dir should be 0.08 1/nm to have 10% strain variation within
%the inner ring of spots, calculated based upon latpar

for i = 1:tot
    curgvec = gvec_arr(i,:);
    curhvec = hvec_arr(i,:);
    newgvec_list((i-1)*numadd+1,:) = curgvec;
    newhvec_list((i-1)*numadd+1,:) = curhvec;
    
    gvec_add_1 = curgvec + [dist, 0, 0];
    gvec_add_2 = curgvec + [2*dist, 0, 0];
    gvec_add_3 = curgvec + [-dist, 0, 0];
    gvec_add_4 = curgvec + [-2*dist, 0, 0];

    gvec_add_5 = curgvec + [0, dist, 0];
    gvec_add_6 = curgvec + [0, 2*dist, 0];
    gvec_add_7 = curgvec + [0, -dist, 0];
    gvec_add_8 = curgvec + [0, -2*dist, 0];

    gvec_add_9 = curgvec + [0, 0, dist];
    gvec_add_10 = curgvec + [0, 0, 2*dist];
    gvec_add_11 = curgvec + [0, 0, -dist];
    gvec_add_12 = curgvec + [0, 0, -2*dist];
    
    newgvec_list((i-1)*numadd+2,:) = gvec_add_1;
    newgvec_list((i-1)*numadd+3,:) = gvec_add_2;
    newgvec_list((i-1)*numadd+4,:) = gvec_add_3;
    newgvec_list((i-1)*numadd+5,:) = gvec_add_4;
    newgvec_list((i-1)*numadd+6,:) = gvec_add_5;
    newgvec_list((i-1)*numadd+7,:) = gvec_add_6;
    newgvec_list((i-1)*numadd+8,:) = gvec_add_7;
    newgvec_list((i-1)*numadd+9,:) = gvec_add_8;
    newgvec_list((i-1)*numadd+10,:) = gvec_add_9;
    newgvec_list((i-1)*numadd+11,:) = gvec_add_10;
    newgvec_list((i-1)*numadd+12,:) = gvec_add_11;
    newgvec_list((i-1)*numadd+13,:) = gvec_add_12;
    
    newhvec_list((i-1)*numadd+2,:) = (g2h_matr*gvec_add_1')';
    newhvec_list((i-1)*numadd+3,:) = (g2h_matr*gvec_add_2')';
    newhvec_list((i-1)*numadd+4,:) = (g2h_matr*gvec_add_3')';
    newhvec_list((i-1)*numadd+5,:) = (g2h_matr*gvec_add_4')';
    newhvec_list((i-1)*numadd+6,:) = (g2h_matr*gvec_add_5')';
    newhvec_list((i-1)*numadd+7,:) = (g2h_matr*gvec_add_6')';
    newhvec_list((i-1)*numadd+8,:) = (g2h_matr*gvec_add_7')';
    newhvec_list((i-1)*numadd+9,:) = (g2h_matr*gvec_add_8')';
    newhvec_list((i-1)*numadd+10,:) = (g2h_matr*gvec_add_9')';
    newhvec_list((i-1)*numadd+11,:) = (g2h_matr*gvec_add_10')';
    newhvec_list((i-1)*numadd+12,:) = (g2h_matr*gvec_add_11')';
    newhvec_list((i-1)*numadd+13,:) = (g2h_matr*gvec_add_12')';
end

gvec_list = newgvec_list;
hvec_list = newhvec_list;

end