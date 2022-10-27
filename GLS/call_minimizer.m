function [solution,fval,exitflag,lambda_] = call_minimizer(B,bfv,le,nelements,ud,mechanical,material,Mass_sigma,C,area)
    
    [LHS,rhs, LHS_lip,rhs_lip] =  constraint(B,bfv,le,nelements,ud,area);
    
    distancewithgrad_ = @(x)  distance_LS(x,B,material,Mass_sigma,nelements,C,area,bfv);
    hess = @(x,lambda) hessian_LS(x,B,Mass_sigma,nelements,C,area);
    
    options = optimoptions('fmincon','Display','none','Algorithm','sqp','ConstraintTolerance',1E-9,'MaxIterations',1E5,...
        'StepTolerance',1E-9,'SpecifyObjectiveGradient',true,'HessianFcn',hess,'MaxFunctionEvaluations',1E5);

    lb = -inf*ones(2*nelements+1,1); ub = inf*ones(2*nelements+1,1);
    lb(1,1) = ud(1,1); ub(1,1) = ud(1,1);
    lb(nelements+1,1) = ud(2,1); ub(nelements+1,1) = ud(2,1);
    
    [solution,fval,exitflag,~,lambda_,~,~] = fmincon(distancewithgrad_,mechanical,LHS_lip,...
        rhs_lip,[],[],lb,ub,[],options);
        w = warning('query','last');
        id = w.identifier;
        warning('off',id);
end





