function result = compute_maxloc(matrix)
    [rows, cols] = size(matrix);
    result = matrix;
    for r = 1:rows
        for c = 1:cols
            neighbors = [
                matrix(max(r-1, 1), c), ... 
                matrix(min(r+1, rows), c), ... 
                matrix(r, max(c-1, 1)), ... 
                matrix(r, min(c+1, cols)) ...
                matrix(r,c)
            ];
            result(r, c) = max(neighbors);
        end
    end
end
