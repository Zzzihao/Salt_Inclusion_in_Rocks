function A = XZ_mean3D(B)
    A = 0.25 * (  B(1:end-1,:,1:end-1) + B(1:end-1,:,2:end) + ...
        B(2:end,:,1:end-1) + B(2:end,:,2:end)  ); 
end