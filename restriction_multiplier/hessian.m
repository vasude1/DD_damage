function hess = hessian(x,B,Mass_sigma,nelements,C,area,C_sig_lag)
    C_sig_lag_ar = (area.*C_sig_lag')';
    nnodes = nelements+1;
    hess = zeros(2*nnodes+nelements,2*nnodes+nelements);
    hess(1:nnodes,1:nnodes) = C*mtimes(B',mtimes( mtimes(Mass_sigma,diag(area)),B));
    hess(nnodes+1:nnodes+nelements,nnodes+1:nnodes+nelements) = 1.0/C*  mtimes(Mass_sigma,diag(area)); 
    hess(nnodes+1:nnodes+nelements,nnodes+nelements+1:2*nnodes+nelements) = C_sig_lag_ar'; 
    hess(nnodes+nelements+1:2*nnodes+nelements,nnodes+1:nnodes+nelements) = hess(nnodes+1:nnodes+nelements,nnodes+nelements+1:2*nnodes+nelements)'; 
end
