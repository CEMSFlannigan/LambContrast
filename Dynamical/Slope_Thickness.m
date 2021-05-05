%top_inv_line = @(x) top_inv_slope*(x-curx) + cury;
            %bot_inv_line = @(x) bot_inv_slope*(x-curx) + cury;
            
            %top_rr_top_in_x = fzero(@(x) top_inv_line(x) - y_top(x), curx);
            %top_rr_top_in_y = y_top(top_rr_top_in_x);
            %top_rr_bot_in_x = fzero(@(x) top_inv_line(x) - y_bot(x), curx);
            %top_rr_bot_in_y = y_bot(top_rr_bot_in_x);
            %top_rr_thickness = sqrt((top_rr_top_in_x - top_rr_bot_in_x)^2 + (top_rr_top_in_y - top_rr_bot_in_y)^2);
            
            %bot_rr_top_in_x = fzero(@(x) bot_inv_line(x) - y_top(x), curx);
            %bot_rr_top_in_y = y_top(bot_rr_top_in_x);
            %bot_rr_bot_in_x = fzero(@(x) bot_inv_line(x) - y_bot(x), curx);
            %bot_rr_bot_in_y = y_bot(bot_rr_bot_in_x);
            %bot_rr_thickness = sqrt((bot_rr_top_in_x - bot_rr_bot_in_x)^2 + (bot_rr_top_in_y - bot_rr_bot_in_y)^2);
            