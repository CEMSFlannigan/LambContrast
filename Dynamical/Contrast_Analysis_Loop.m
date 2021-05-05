beta_add = -1:0.05:1; %0.01
alpha_add = -5:0.25:5; %0.01

latpar = 0.5658; % nm
x_spacing = 50*latpar; % nm
xvals = 0:x_spacing:1000; % nm

bf_contrast_array = ones(length(xvals),length(beta_add)*length(alpha_add));

for i = 1:length(beta_add)
    for j = 1:length(alpha_add)
        input_array = {beta_add(i), alpha_add(j)};
        
        [bf_contrast,strain_matr,thickness_vector] = Contrast_Analysis_Parallel_Method(input_array);
        
        bf_contrast_array(:,(i-1)*length(beta_add)+j) = bf_contrast;
    end
end

assignin('base','bf_contrast',bf_contrast);

save('tiltvalues.mat');