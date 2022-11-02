function [LHS,rhs, LHS_lip,rhs_lip] =  constraint(B,bfv,le,nelements,ud,area)
    nnodes=nelements+1;
      
    LHS = zeros(2*(nnodes-2),2*nnodes+nelements);
    rhs = zeros(2*(nnodes-2),1);
    area_B = area.*B;
    C_int = le*area_B(1:end-1,1:end-1);
    P = mtimes(C_int,B);
    
    LHS(1:(nnodes-2),1:nnodes) = P; %le*
    LHS((nnodes-2)+1:2*(nnodes-2),nnodes+nelements+1:2*nnodes+nelements) = P; %le*

    % rhs(1:nnodes-2,1) = bfv(2:end-1,1);
    rhs(1:(nnodes-2),1) = -le*0.0;
    rhs((nnodes-2)+1:end,1) = 0.0;
    % rhs(nnodes+1,1) = ud(1,1); rhs(nnodes+2,1) = ud(2,1);
    
    % Inequality constraint
    LHS_lip = zeros(4*(nnodes-2),2*nnodes+nelements);
    rhs_lip = zeros(4*(nnodes-2),1);
    inter_mat = mtimes(B',B);
    
    LHS_lip(1:nnodes-2,1:nnodes) = inter_mat(2:end-1,:);
    LHS_lip(nnodes-1:2*(nnodes-2),1:nnodes) = -inter_mat(2:end-1,:);
    
    LHS_lip(2*(nnodes-2)+1:3*(nnodes-2),nnodes+nelements+1:2*nnodes+nelements) = inter_mat(2:end-1,:);
    LHS_lip(3*(nnodes-2)+1:4*(nnodes-2),nnodes+nelements+1:2*nnodes+nelements) = -inter_mat(2:end-1,:);
    
    rhs_lip(:,1) = 50;
end