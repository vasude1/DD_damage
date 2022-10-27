function [distance_,grad] = distance_LS(mechanical,B,material,Mass_sigma,nelements,C,area,bfv)

    strain_mech = mtimes(B,mechanical(1:nelements+1,1));
    stress_mech = mechanical(nelements+2:end)';
    strain_mat = material(1:nelements)';
    stress_mat = material(nelements+1:2*nelements)';
    deltastress = stress_mech'-stress_mat';
    deltastrain = strain_mech-strain_mat';
    dist1 = 1.0/2.0/C*mtimes(deltastress',mtimes(Mass_sigma,area.*deltastress));
    dist2 = 1.0/2.0*C* mtimes(deltastrain',mtimes(Mass_sigma,area.*deltastrain));
    weight = 1E-1;
    bbt = mtimes(B,B');
    dist3 = weight/2.0*mtimes(stress_mech, mtimes(bbt,stress_mech'))+mtimes(bfv',mtimes(B',stress_mech'));
    distance_= dist1+dist2+dist3;   

    if nargout > 1 % gradient required
        grad = zeros(2*nelements+1,1);
        grad(1:nelements+1,1) = C*mtimes(B',mtimes(Mass_sigma,area.*deltastrain));
        grad(nelements+2:end,1) = 1.0/C*mtimes(Mass_sigma,area.*deltastress);
        grad(nelements+2:end,1) = grad(nelements+2:end,1)+...
            weight*mtimes(bbt,stress_mech')+mtimes(B,bfv);
%         if nargout > 2 % Hessian required
%             hess = zeros(2*nelements+1,2*nelements+1);
%             hess(1:nnodes_disp,1:nnodes_disp) = C*mtimes(B',mtimes(Mass_sigma,B));
%             hess(nnodes_disp+1:end,nnodes_disp+1:end) = 1.0/C* Mass_sigma; 
%         end
    end
end