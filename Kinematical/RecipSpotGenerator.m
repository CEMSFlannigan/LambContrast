function [gvec_list,hvec_list] = RecipSpotGenerator(count,gvec_list,lim,hvec_list)

curcount = count;
if curcount >= lim
elseif curcount < lim
    for i = 1:curcount
        curgvec = gvec_list{i};
        curhvec = hvec_list{i};
        
        newvec1 = curgvec + gvec_list{1};
        newhvec1 = curhvec + hvec_list{1};
        newvec2 = curgvec + gvec_list{2};
        newhvec2 = curhvec + hvec_list{2};
        newvec3 = curgvec + gvec_list{3};
        newhvec3 = curhvec + hvec_list{3};
        newvec4 = curgvec + gvec_list{4};
        newhvec4 = curhvec + hvec_list{4};

        find1 = 0;
        find2 = 0;
        find3 = 0;
        find4 = 0;
        
        for j = 1:length(hvec_list)
            if isequal(hvec_list{j},newhvec1) || isequal(newhvec1,[0,0,0])
                find1 = 1;
            end
        end
        
        if find1 == 1 || count > lim
        else
            count = count + 1;
            gvec_list{count} = newvec1;
            hvec_list{count} = newhvec1;
        end
        
        for j = 1:length(hvec_list)
            if isequal(hvec_list{j},newhvec2) || isequal(newhvec2,[0,0,0])
                find2 = 1;
            end
        end
        
        if find2 == 1 || count > lim
        else
            count = count + 1;
            gvec_list{count} = newvec2;
            hvec_list{count} = newhvec2;
        end
        
        for j = 1:length(hvec_list)
            if isequal(hvec_list{j},newhvec3) || isequal(newhvec3,[0,0,0])
                find3 = 1;
            end
        end
        
        if find3 == 1 || count > lim
        else
            count = count + 1;
            gvec_list{count} = newvec3;
            hvec_list{count} = newhvec3;
        end
        
        for j = 1:length(hvec_list)
            if isequal(hvec_list{j},newhvec4) || isequal(newhvec4,[0,0,0])
                find4 = 1;
            end
        end
        
        if find4 == 1 || count > lim
        else
            count = count + 1;
            gvec_list{count} = newvec4;
            hvec_list{count} = newhvec4;
        end
        
    end
    
    [gvec_list,hvec_list] = RecipSpotGenerator(count,gvec_list,lim,hvec_list);
else
    count = lim;
end

end