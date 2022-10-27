function [solu,fval,exitflag,lambda_] = call_minimizer(B,bfv,le,nelements,ud,mechanical,material,Mass_sigma,C,area,C_sig_lag)
    nnodes=nelements+1;
    
    init_guess = zeros(2*nnodes+nelements,1);
    init_guess(1:nnodes+nelements,1) = mechanical(:,1);
    
    [LHS,rhs, LHS_lip,rhs_lip] =  constraint(B,bfv,le,nelements,ud,area);
    
    distancewithgrad_ = @(x)  distance(x,B,material,Mass_sigma,nelements,C,area,C_sig_lag);
    hess = @(x,lambda) hessian(x,B,Mass_sigma,nelements,C,area,C_sig_lag);

    lb = -inf*ones(2*nnodes+nelements,1); ub = inf*ones(2*nnodes+nelements,1);
    lb(1,1) = ud(1,1); ub(1,1) = ud(1,1);
    lb(nnodes,1) = ud(2,1); ub(nnodes,1) = ud(2,1);
    lb(nnodes+nelements+1,1) = 0.0; ub(nnodes+nelements+1,1) = 0.0;
    lb(2*nnodes+nelements,1) = 0.0; ub(2*nnodes+nelements,1) = 0.0;
    
     options = optimoptions('fmincon','Display','none','Algorithm','sqp','ConstraintTolerance',1E-9,'MaxIterations',1E5,...
        'StepTolerance',1E-9,'SpecifyObjectiveGradient',true,'HessianFcn',hess,'MaxFunctionEvaluations',1E5);
        
    [solu,fval,exitflag,~,lambda_,~,~] = fmincon(distancewithgrad_,init_guess,LHS_lip,...
        rhs_lip,LHS,rhs,lb,ub,[],options);
     
    w = warning('query','last');
    id = w.identifier;
    warning('off',id);
end





