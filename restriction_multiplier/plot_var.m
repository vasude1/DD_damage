iter_ = step-1;
var=2;
x = linspace(le/2,1.0-(le/2.0),nelements);
% 
yyaxis left;
plot(x,elem_mech(:,iter_,var),'r-*');
hold on
plot(x,elem_mat(:,iter_,var),'b-*');

% 
% element = floor(nelements/2)+1;
% plot(elem_mech(element,1:iter_,1),elem_mech(element,1:iter_,2),'r-*')
% hold on
% plot(elem_mat(element,1:iter_,1),elem_mat(element,1:iter_,2),'b-*')


% yyaxis right;
% lag = elem_lagrange(1:nnodes-2,iter_,1) - elem_lagrange(nnodes-1:2*nnodes-4,iter_,1);
% bfv_lag1 = zeros(nnodes-1,1); bfv_lag1(1:end-1,1) = lag;
% bfv_lag2 = zeros(nnodes-1,1); bfv_lag2(2:end,1) = lag;
% bfv_lag = 1.0/le*(bfv_lag1-bfv_lag2);
% area_B = area.*B_u;
% BT_short = area_B'; 
% BT_bfv_lag = mtimes(BT_short,bfv_lag);
% % 
% plot(x,bfv_lag,'g-*');
% 
% disp(exit_flags(iter_,1));


% nod_coords = linspace(0,length,nnodes);
% plot(nod_coords,lag_mult_equ(:,iter_,1));