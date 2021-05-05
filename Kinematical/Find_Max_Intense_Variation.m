curmax = 0;
foundmax = 0;

for k = 1:37%length(timevals)
    for k2 = 1:37%length(timevals)
        bf_diff = abs(bf_contrast(:,k) - bf_contrast(:,k2));
        curmax = max(bf_diff);
        
        if curmax > foundmax
            foundmax = curmax
        end
    end
end
