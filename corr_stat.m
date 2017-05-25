function correlation  = corr_stat (slice_combo)

slice_1 = slice_combo(:, 1);
slice_2 = slice_combo(:, 2);
stats = 0;

 for i  = 1:length(slice_1)
     if slice_1(i) == 1 && slice_2(i) == 1
         stats = stats + 1;
     end
 end
 correlation = stats;