function psi_itp_2D(psi_itp, zval, X)
    x = X[1]; y = X[2]; 
    Nx = length(x); Ny = length(y); 
    psi_itp2D = im.*zeros(Nx, Ny);
    for i in 1:Nx
        for j in 1:Ny
                psi_itp2D[i, j] = psi_itp(x[i], y[j], zval);
        end
    end
    return psi_itp2D;   
end
