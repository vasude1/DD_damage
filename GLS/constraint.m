function [LHS,rhs, LHS_lip,rhs_lip] =  constraint(B,bfv,le,nelements,ud,area)
    nnodes=nelements+1;
      
    % Equilibrium as equality constraint and first and lasr rows as
    % dirichlet
    LHS = zeros(nnodes-2,nnodes+nelements);
    rhs = zeros(nnodes-2,1);
    area_B = area.*B;
    BT_short = area_B';    
    LHS(1:nnodes-2,nnodes+1:nnodes+nelements) = le*BT_short(2:end-1,:); %le*
%     LHS(nnodes+1,1) = 1.0;
%     LHS(nnodes+2,nnodes) = 1.0;

    rhs(1:nnodes-2,1) = bfv(2:end-1,1);
    % rhs(nnodes+1,1) = ud(1,1); rhs(nnodes+2,1) = ud(2,1);
    
    % Inequality constraint
    LHS_lip = zeros(2*nnodes-4,nnodes+nelements);
    rhs_lip = zeros(2*nnodes-4,1);
    inter_mat = mtimes(B',B);
    LHS_lip(1:nnodes-2,1:nnodes) = inter_mat(2:end-1,:);
    LHS_lip(nnodes-1:2*nnodes-4,1:nnodes) = -inter_mat(2:end-1,:);
    rhs_lip(:,1) = 15.0;
    
end