function result = compute_maxloc3D(matrix)
    [rows, cols, deps] = size(matrix);
    result = matrix;
    for r = 1:rows
        for c = 1:cols
            for d = 1:deps
                neighbors = [
                    matrix(max(r-1, 1), c, d), ... 
                    matrix(min(r+1, rows), c, d), ... 
                    matrix(r, max(c-1, 1), d), ... 
                    matrix(r, min(c+1, cols), d), ...
                    matrix(r, c, max(d-1, 1)), ...
                    matrix(r, c, min(d+1, deps))
                ];
                result(r, c, d) = max(neighbors);
        end
    end
end
