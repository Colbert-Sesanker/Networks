function p=cor(a,b)
  c = corrcoef(a,b);
  p = c(2,1);
end
