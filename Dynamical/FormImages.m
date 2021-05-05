% max_bf_contrast = max(max(bf_contrast(:,1:33)));
% min_bf_contrast = min(min(bf_contrast(:,1:33)));
% 
% ex_image = cell(1,length(timevals));
% example_depth = 5000; % nm
% 
% for k = 1:33%length(timevals)
%     cur_img = zeros(example_depth,round((length(xvals))*x_spacing)+3);%);
%     append_img = ones(size(cur_img(1,4:end)));
%     
%      for i = 1:length(xvals)
%          append_img(1,round((i-1)*x_spacing)+1:(round((i)*x_spacing))+1) = bf_contrast(i,k);
%      end
%     
%     for l = 1:example_depth
%         cur_img(l,1) = min_bf_contrast;
%         cur_img(l,2) = max_bf_contrast;
%         cur_img(l,3:end) = append_img(1,:);
%     end
%     cur_img(cur_img < 0) = 0;
%     cur_img(cur_img > 1) = 1;
%     ex_image{k} = cur_img;
% end

max_bf_contrast = max(max(bf_contrast));%max([max(max(bf_contrast_5A_2B(:,1:33))) max(max(bf_contrast_n13p02A_1p55B(:,1:33))) max(max(bf_contrast_0A_0B(:,1:33)))]);
min_bf_contrast = min(min(bf_contrast));%min([min(min(bf_contrast_5A_2B(:,1:33))) min(min(bf_contrast_n13p02A_1p55B(:,1:33))) min(min(bf_contrast_0A_0B(:,1:33)))]);

ex_image = cell(1,length(timevals));
example_depth = 5000; % nm

for k = 1:length(timevals)
    cur_img = zeros(example_depth,round((length(xvals))*x_spacing)+3);%);
    append_img = ones(size(cur_img(1,4:end)));
    
     for i = 1:length(xvals)
         append_img(1,round((i-1)*x_spacing)+1:(round((i)*x_spacing))+1) = bf_contrast(i,k);
     end
    
    for l = 1:example_depth
        cur_img(l,1) = min_bf_contrast;
        cur_img(l,2) = max_bf_contrast;
        cur_img(l,3:end) = append_img(1,:);
    end
    cur_img(cur_img < 0) = 0;
    cur_img(cur_img > 1) = 1;
    ex_image{k} = cur_img;
end