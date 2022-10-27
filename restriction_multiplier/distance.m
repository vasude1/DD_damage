function [distance_,grad] = distance(mechanical,B,material,Mass_sigma,nelements,C,area,C_sig_lag)
    C_sig_lag_ar = (area.*C_sig_lag')';
    nnodes = nelements+1;
    strain_mech = mtimes(B,mechanical(1:nnodes,1));
    stress_mech = mechanical(nnodes+1:nnodes+nelements)';
    strain_mat = material(1:nelements)';
    stress_mat = material(nelements+1:2*nelements)';
    deltastress = stress_mech'-stress_mat';
    deltastrain = strain_mech-strain_mat';
    dist1 = 1.0/2.0/C*mtimes(deltastress',mtimes(Mass_sigma,area.*deltastress));
    dist2 = 1.0/2.0*C* mtimes(deltastrain',mtimes(Mass_sigma,area.*deltastrain));
    dist3 = mtimes(mechanical(nnodes+nelements+1:2*nnodes+nelements,1)', mtimes(C_sig_lag_ar, stress_mech' ));
    distance_= dist1+dist2+dist3;

    if nargout > 1 % gradient required
        grad = zeros(2*nnodes+nelements,1);
        grad(1:nnodes,1) = C*mtimes(B',mtimes(Mass_sigma,area.*deltastrain));
        grad(nnodes+1:nnodes+nelements,1) = 1.0/C*mtimes(Mass_sigma,area.*deltastress)+...
            mtimes(C_sig_lag_ar', mechanical(nnodes+nelements+1:2*nnodes+nelements));
        grad(nnodes+nelements+1:2*nnodes+nelements) = mtimes(C_sig_lag_ar, stress_mech');
%         if nargout > 2 % Hessian required
%             hess = zeros(2*nelements+1,2*nelements+1);
%             hess(1:nnodes_disp,1:nnodes_disp) = C*mtimes(B',mtimes(Mass_sigma,B));
%             hess(nnodes_disp+1:end,nnodes_disp+1:end) = 1.0/C* Mass_sigma;
%         end
    end
end
