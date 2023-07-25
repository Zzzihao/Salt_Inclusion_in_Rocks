function A = Z_mean(B)
    A = 0.5 * ( B(:,1:end-1) + B(:,2:end) );
end