function A2 = smooth(A, factor)
    [m, n, h] = size(A);
    A2 = zeros(m, n, h);
    
    for i = 2:m-1
        for j = 2:n-1
            for k = 2:h-1
                A2(i,j,k) = A(i,j,k) + 1.0/6.1/factor*(...
                    ( ( A(i+1,j,k)-A(i,j,k) ) - ( A(i,j,k)-A(i-1,j,k) ) ) + ...
                    ( ( A(i,j+1,k)-A(i,j,k) ) - ( A(i,j,k)-A(i,j-1,k) ) ) + ...
                    ( ( A(i,j,k+1)-A(i,j,k) ) - ( A(i,j,k)-A(i,j,k-1) ) ) );
            end
        end
    end
    
    % Copy boundary values
    A2(1,:,:) = A(1,:,:);
    A2(end,:,:) = A(end,:,:);
    A2(:,1,:) = A(:,1,:);
    A2(:,end,:) = A(:,end,:);
    A2(:,:,1) = A(:,:,1);
    A2(:,:,end) = A(:,:,end);
end