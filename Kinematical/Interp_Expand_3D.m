Alph_Tilts_Interp = linspace(min(Alph_Tilts),max(Alph_Tilts),length(Alph_Tilts)*10);
Bet_Tilts_Interp = linspace(min(Bet_Tilts),max(Bet_Tilts),length(Bet_Tilts)*10);
Max_Osc_Interp = zeros(length(Max_Oscillations)*10,length(Max_Oscillations)*10);

for i = 1:length(Alph_Tilts_Interp)
    cur_Alph = Alph_Tilts_Interp(i);
    for j = 1:length(Bet_Tilts_Interp)
        cur_Bet = Bet_Tilts_Interp(j);
        
        Closest_Distance_1 = 1e10;
        Closest_Distance_2 = 1e100;
        Closest_Distance_3 = 1e1000;
        Closest_Distance_4 = 1e10000;
        
        idx1 = 0;
        idx2 = 0;
        idx3 = 0;
        idx4 = 0;
        
        for k = 1:length(Alph_Tilts)
            curdist = sqrt((Alph_Tilts(k) - cur_Alph)^2+(Bet_Tilts(k) - cur_Bet)^2);
            if curdist < Closest_Distance_1
                Closest_Distance_4 = Closest_Distance_3;
                Closest_Distance_3 = Closest_Distance_2;
                Closest_Distance_2 = Closest_Distance_1;
                Closest_Distance_1 = curdist;
                
                idx4 = idx3;
                idx3 = idx2;
                idx2 = idx1;
                idx1 = k;
            elseif curdist < Closest_Distance_2 && curdist > Closest_Distance_1
                Closest_Distance_4 = Closest_Distance_3;
                Closest_Distance_3 = Closest_Distance_2;
                Closest_Distance_2 = curdist;
                
                idx4 = idx3;
                idx3 = idx2;
                idx2 = k;
            elseif curdist < Closest_Distance_3 && curdist > Closest_Distance_2
                Closest_Distance_4 = Closest_Distance_3;
                Closest_Distance_3 = curdist;
                
                idx4 = idx3;
                idx3 = k;
            elseif curdist < Closest_Distance_4 && curdist > Closest_Distance_3
                Closest_Distance_4 = curdist;
                
                idx4 = k;
            end
        end
        plane_points = [Alph_Tilts(idx1), Bet_Tilts(idx1), Max_Oscillations(idx1); Alph_Tilts(idx2), Bet_Tilts(idx2), Max_Oscillations(idx2); Alph_Tilts(idx3), Bet_Tilts(idx3), Max_Oscillations(idx3); Alph_Tilts(idx4), Bet_Tilts(idx4), Max_Oscillations(idx4)]
        plane_constants = [ones(4,1), plane_points(:,1:2)] \ plane_points(:,3) % z = ax + by + c, [c, a, b]
        Max_Osc_Interp(i,j) = plane_constants(2)*cur_Alph + plane_constants(3)*cur_Bet + plane_constants(1);
    end
end