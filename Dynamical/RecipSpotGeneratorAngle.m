function [gvec_list_h,gvec_list_k,gvec_list_l,hvec_list_h,hvec_list_k,hvec_list_l,count] = RecipSpotGeneratorAngle(angled,count,gvec_list_h,gvec_list_k,gvec_list_l,hvec_list_h,hvec_list_k,hvec_list_l, g_unit, g_unit_2, g_mags, g_mags_2, g_basis_1, g_basis_2)

lambda = 2.5079e-3*1e-9; % m, 200 keV

angle = angled;
curgvec = [gvec_list_h(end), gvec_list_k(end), gvec_list_l(end)];
angles = @(x,y,z) asind(lambda*1e9/2*sqrt(x.*x+y.*y+z.*z));
curangle = angles(curgvec(1), curgvec(2), curgvec(3));
curcount = count;

if curangle > angle
elseif curangle < angle
    for i = 1:curcount
        curgvec = [gvec_list_h(i), gvec_list_k(i), gvec_list_l(i)];
        curhvec = [hvec_list_h(i), hvec_list_k(i), hvec_list_l(i)];
        
        newvec1 = curgvec + g_unit.*g_mags;
        newhvec1 = curhvec + g_basis_1;
        newvec2 = curgvec + g_unit_2.*g_mags_2;
        newhvec2 = curhvec + g_basis_2;
        newvec3 = curgvec - g_unit.*g_mags;
        newhvec3 = curhvec - g_basis_1;
        newvec4 = curgvec - g_unit_2.*g_mags_2;
        newhvec4 = curhvec - g_basis_2;

        find1 = 0;
        find2 = 0;
        find3 = 0;
        find4 = 0;
        
        for j = 1:length(hvec_list_h)
            if isequal([hvec_list_h(j), hvec_list_k(j), hvec_list_l(j)],newhvec1) || isequal(newhvec1,[0,0,0])
                find1 = 1;
            end
        end
        
        if find1 == 1 || angles(newvec1(1), newvec1(2), newvec1(3)) > angle
        else
            count = count + 1;
            gvec_list_h = [gvec_list_h newvec1(1)];
            gvec_list_k = [gvec_list_k newvec1(2)];
            gvec_list_l = [gvec_list_l newvec1(3)];
            hvec_list_h = [hvec_list_h newhvec1(1)];
            hvec_list_k = [hvec_list_k newhvec1(2)];
            hvec_list_l = [hvec_list_l newhvec1(3)];
            
            [gvec_list_h,gvec_list_k,gvec_list_l,hvec_list_h,hvec_list_k,hvec_list_l,count] = RecipSpotGeneratorAngle(angled,count,gvec_list_h,gvec_list_k,gvec_list_l,hvec_list_h,hvec_list_k,hvec_list_l, g_unit, g_unit_2, g_mags, g_mags_2, g_basis_1, g_basis_2);
        end
        
        for j = 1:length(hvec_list_h)
            if isequal([hvec_list_h(j), hvec_list_k(j), hvec_list_l(j)],newhvec2) || isequal(newhvec2,[0,0,0])
                find2 = 1;
            end
        end
        
        if find2 == 1 || angles(newvec2(1), newvec2(2), newvec2(3)) > angle
        else
            count = count + 1;
            gvec_list_h = [gvec_list_h newvec2(1)];
            gvec_list_k = [gvec_list_k newvec2(2)];
            gvec_list_l = [gvec_list_l newvec2(3)];
            hvec_list_h = [hvec_list_h newhvec2(1)];
            hvec_list_k = [hvec_list_k newhvec2(2)];
            hvec_list_l = [hvec_list_l newhvec2(3)];
            
            [gvec_list_h,gvec_list_k,gvec_list_l,hvec_list_h,hvec_list_k,hvec_list_l,count] = RecipSpotGeneratorAngle(angled,count,gvec_list_h,gvec_list_k,gvec_list_l,hvec_list_h,hvec_list_k,hvec_list_l, g_unit, g_unit_2, g_mags, g_mags_2, g_basis_1, g_basis_2);
        end
        
        for j = 1:length(hvec_list_h)
            if isequal([hvec_list_h(j), hvec_list_k(j), hvec_list_l(j)],newhvec3) || isequal(newhvec3,[0,0,0])
                find3 = 1;
            end
        end
        
        if find3 == 1 || angles(newvec3(1), newvec3(2), newvec3(3)) > angle
        else
            count = count + 1;
            gvec_list_h = [gvec_list_h newvec3(1)];
            gvec_list_k = [gvec_list_k newvec3(2)];
            gvec_list_l = [gvec_list_l newvec3(3)];
            hvec_list_h = [hvec_list_h newhvec3(1)];
            hvec_list_k = [hvec_list_k newhvec3(2)];
            hvec_list_l = [hvec_list_l newhvec3(3)];
            
            [gvec_list_h,gvec_list_k,gvec_list_l,hvec_list_h,hvec_list_k,hvec_list_l,count] = RecipSpotGeneratorAngle(angled,count,gvec_list_h,gvec_list_k,gvec_list_l,hvec_list_h,hvec_list_k,hvec_list_l, g_unit, g_unit_2, g_mags, g_mags_2, g_basis_1, g_basis_2);
        end
        
        for j = 1:length(hvec_list_h)
            if isequal([hvec_list_h(j), hvec_list_k(j), hvec_list_l(j)],newhvec4) || isequal(newhvec4,[0,0,0])
                find4 = 1;
            end
        end
        
        if find4 == 1 || angles(newvec4(1), newvec4(2), newvec4(3)) > angle
        else
            count = count + 1;
            gvec_list_h = [gvec_list_h newvec4(1)];
            gvec_list_k = [gvec_list_k newvec4(2)];
            gvec_list_l = [gvec_list_l newvec4(3)];
            hvec_list_h = [hvec_list_h newhvec4(1)];
            hvec_list_k = [hvec_list_k newhvec4(2)];
            hvec_list_l = [hvec_list_l newhvec4(3)];
            
            [gvec_list_h,gvec_list_k,gvec_list_l,hvec_list_h,hvec_list_k,hvec_list_l,count] = RecipSpotGeneratorAngle(angled,count,gvec_list_h,gvec_list_k,gvec_list_l,hvec_list_h,hvec_list_k,hvec_list_l, g_unit, g_unit_2, g_mags, g_mags_2, g_basis_1, g_basis_2);
        end
        
    end
    
else
    curangle = angle;
end

end