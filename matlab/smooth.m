function A2 = smooth(A, factor)
    [m, n] = size(A);
    A2 = zeros(m, n);
    
    for i = 2:m-1
        for j = 2:n-1
            A2(i, j) = A(i, j) + 1.0/4.1/factor * (A(i-1, j) - 2 * A(i, j) + A(i+1, j) + A(i, j-1) - 2 * A(i, j) + A(i, j+1));
        end
    end
    
    % Copy boundary values
    A2(1, :) = A(1, :);
    A2(end, :) = A(end, :);
    A2(:, 1) = A(:, 1);
    A2(:, end) = A(:, end);
end
