function [idxs, closest_dif_y, closest_dif_x] = findClosest(matx,maty,x,y)

difx = matx - ones(size(matx)).*x;

idxs = [1 1];
closest_dif_y = 1e100;
closest_dif_x = 1e100;
closest_y = 1e100;

for i = 1:size(matx,2)
    pot_indx = find(diff(sign(difx(:,i))))+1;
    pot_indx = pot_indx(end);
    if abs(matx(pot_indx,i) - x) > abs(matx(pot_indx-1,i)-x)
        pot_indx = pot_indx-1;
    end
    
    for j = 1:size(maty,2)
        if abs(maty(pot_indx,j) - y) < abs(closest_dif_y)
            closest_dif_y = maty(pot_indx,j) - y;
            closest_dif_x = difx(pot_indx,j);
            closest_y = maty(pot_indx,j);
            idxs = [pot_indx j];
        end
    end
    
end

end