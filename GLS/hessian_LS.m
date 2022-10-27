function hess = hessian_LS(x,B,Mass_sigma,nelements,C,area)
    weight = 1E-1;
    nnodes_sigma = nelements;
    nnodes_disp = nelements+1;
    hess = zeros(nnodes_disp+nnodes_sigma,nnodes_disp+nnodes_sigma);
    hess(1:nnodes_disp,1:nnodes_disp) = C*mtimes(B',mtimes( mtimes(Mass_sigma,diag(area)),B));
    hess(nnodes_disp+1:nnodes_disp+nnodes_sigma,nnodes_disp+1:nnodes_disp+nnodes_sigma) = 1.0/C*  mtimes(Mass_sigma,diag(area))+...
        weight*mtimes(B,B'); 
end
