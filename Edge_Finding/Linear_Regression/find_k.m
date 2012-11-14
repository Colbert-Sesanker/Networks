% This function finds hill function's k's from data

function k = find_k(target, source1, source2, spline_step)
target_diff = diff(target);
d = map(@sign, target_diff);
true_min = (find(d==1)(1) >= ceil(50/spline_step)); % true if first min is a a true min and not a reponse to initial condition 
d2 = diff(d);
if d(1) == 1 || (d(1) == -1 && true_min),    
  min2 = find(d2==2)(1) + 1     % index of second minimum in target  
  k = get_k(target, min2, source1, source2);
else
  min2 = find(d2==2)(2) + 1; 
  k = get_k(target, min2, source1, source2);
end
endfunction

function k = get_k(target,  min2, source1, source2)
m = diff(target(min2:end));
m = map(@sign, m);
max2 = find(m==-1)(1) + min2 - 1
k_center_idx = floor(length(target(min2:max2))/2) + min2
k = zeros(2,5);
for n = 0 : 4
  k(1, n+1) = source1(k_center_idx - 4 + 2*n);
  k(2, n+1) = source2(k_center_idx - 4 + 2*n);
end
endfunction
