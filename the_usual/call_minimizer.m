function [solu,fval,exitflag,lambda_] = call_minimizer(B,bfv,le,nelements,ud,mechanical,material,Mass_sigma,C,area)
    nnodes=nelements+1;
    [LHS,rhs, LHS_lip,rhs_lip] =  constraint(B,bfv,le,nelements,ud,area);
    
    distancewithgrad_ = @(x)  distance(x,B,material,Mass_sigma,nelements,C,area);
    hess = @(x,lambda) hessian(x,B,Mass_sigma,nelements,C,area);

    lb = -inf*ones(2*nelements+1,1); ub = inf*ones(2*nelements+1,1);
    lb(1,1) = ud(1,1); ub(1,1) = ud(1,1);
    lb(nelements+1,1) = ud(2,1); ub(nelements+1,1) = ud(2,1);
    
     options = optimoptions('fmincon','Display','none','Algorithm','sqp','ConstraintTolerance',1E-9,'MaxIterations',1E5,...
        'StepTolerance',1E-9,'SpecifyObjectiveGradient',true,'HessianFcn',hess,'MaxFunctionEvaluations',1E5);
        
    [solu,fval,exitflag,~,lambda_,~,~] = fmincon(distancewithgrad_,mechanical,LHS_lip,...
        rhs_lip,LHS,rhs,lb,ub,[],options);
    lag = zeros(nnodes-2,1);
    solution = zeros(nnodes+nelements,1);
    
    while(norm(solu-solution)>1E-8)
        % lag - lambda_.ineqlin(1:nnodes-2,1) + lambda_.ineqlin(nnodes-1:2*nnodes-4,1)
            % Corrector step
         solu(:) = solution(:);
        lag = lambda_.ineqlin(1:nnodes-2,1) - lambda_.ineqlin(nnodes-1:2*nnodes-4,1);
        bfv_lag1 = zeros(nnodes-1,1); bfv_lag1(1:end-1,1) = lag;
        bfv_lag2 = zeros(nnodes-1,1); bfv_lag2(2:end,1) = lag;
        bfv_lag = 1.0/le*(bfv_lag1-bfv_lag2);
%         area_B = area.*B;
%         BT_short = area_B'; 
        BT_bfv_lag = mtimes(B',area.*bfv_lag);

        options = optimoptions('fmincon','Display','none','Algorithm','sqp','ConstraintTolerance',1E-9,'MaxIterations',1E5,...
            'StepTolerance',1E-9,'SpecifyObjectiveGradient',true,'HessianFcn',hess,'MaxFunctionEvaluations',1E5);

        [solution,fval,exitflag,~,lambda_,~,~] = fmincon(distancewithgrad_,solu,LHS_lip,...
            rhs_lip,LHS,rhs+BT_bfv_lag(2:end-1),lb,ub,[],options);
        % solu(:) = solution(:);
        w = warning('query','last');
        id = w.identifier;
        warning('off',id);
        % disp('Here');
    end
    
    % disp('Outta here');
end





