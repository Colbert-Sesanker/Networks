% Returns true if matrix has any NaN entries
function  isTrue = hasNaN(matrix)
isTrue = any(isnan(matrix(:)));
end
