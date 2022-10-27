% Data driven case
% Linear Elements

delta_0 = 1.0;
delta_f = 50.0;
f=0.0; % Body force
delta_t = 1E-2; % Time increment for dynamics

beta=0.45;  % Numerical integration
gamma = 0.8;    % Numerical integration

% Discretization
nnodes = 21;
nsteps=2500;
length = 1.0;
nelements = nnodes-1; ne = nelements;
nele=linspace(1,nelements,nelements);
nodal_coords = linspace(0,length,nnodes);
le=length/nelements;
area = ones(nelements,1);
area(floor(nelements/2),1) = 0.97;

%% Setup arrays for FE
C=1E2;
displacement = zeros(nnodes,1);
velocity = zeros(nnodes,1);
acceleration = zeros(nnodes,1);

prev_disp = zeros(nnodes,1);
prev_velo = zeros(nnodes,1);
prev_acc = zeros(nnodes,1);

% Mechanical states
sigma = zeros(nelements,1);
epsilon = zeros(nelements,1);
prev_sigma = zeros(nelements,1);
prev_epsilon = zeros(nelements,1);

% Material states
sigma_t = ones(nelements,1);
epsilon_t = zeros(nelements,1);
energy = zeros(nelements,1);
damage = zeros(nelements,1);
dissipation = zeros(nelements,1);
C_ = 100.0*ones(nelements,1);
%redo_var = np.zeros((nelements,1))

% Material states at n-1
prev_sigma_t = zeros(nelements,1);
prev_epsilon_t = zeros(nelements,1);
prev_energy = zeros(nelements,1);
prev_damage = zeros(nelements,1);
prev_dissipation = zeros(nelements,1);
node_temp = zeros(nelements,1);

% Temporarily store the mateiral states to look for convergence
sigma_temp = zeros(nelements,1);
epsilon_temp = zeros(nelements,1);

% Store element mech and mat states at all times
elem_mech = zeros(nelements,nsteps,2);
elem_mat = zeros(nelements,nsteps,2);
elem_damage = zeros(nelements,nsteps,2);
elem_dissipation = zeros(nelements,nsteps,1);
elem_lagrange = zeros(2*nnodes-4,nsteps,1);
node_beta = zeros(nelements,nsteps,1);
elem_dist2 = zeros(nelements,nsteps,1);
prev_node = 200*ones(nelements,nsteps,1);
exit_flags = zeros(nsteps,1);

delem_mech = zeros(nelements,nsteps,1); % Why is this used?
delem_mat = zeros(nelements,nsteps,1);  % why is this used?
node_displacement = zeros(nnodes,nsteps,2); % Duplicate array?

% FE parameters
dist2 = zeros(nelements,1); % Stores the distance function at elements
eta = zeros(nnodes,1); % Lagrange multiplier
B_u = zeros(nelements,nnodes);  % Strain displacement
K = zeros(nnodes,nnodes);   %B^T B - Stiffness
RHS1 = zeros(nnodes,1);
LHS2 = zeros(nelements-1,nelements-1);
RHS2 = zeros(nelements-1,1);

for i = 1: nnodes-1
    B_u(i,i) = -1.0/le;
    B_u(i,i+1) = 1.0/le;
end

K(1:nnodes,1:nnodes) = B_u' * B_u;
rho=1.0;
Mass = zeros(nnodes,nnodes);
for i = 1:nnodes-1
    Mass(i:i+1,i:i+1) = Mass(i:i+1,i:i+1)+ rho*le/6.0*[[2.0 1.0];[1.0 2.0]];
end
Mass_sigma = le*eye(nelements);

% Body forces
bfv = f* ones(nnodes,1);
bfv = Mass*bfv;
Mass = zeros(nnodes,nnodes);

% Boundary conditions
u_d = zeros(nsteps,1);
interval = length*linspace(0.8, 0.6,300);
u_d(1:100,1)=length*linspace(0.0, 0.995,100 );
u_d(101:end,1) = length*linspace(0.995, 1.5,nsteps-100 );


%% Setup the solver
converge = false;
redo_step = false;

material = zeros(2*nelements,1);
mechanical_st = zeros(2*nelements+1,1);
ud = zeros(2,1);
A_large = zeros(2*nnodes,nnodes+nelements);
rhs_large = zeros(2*nnodes,1);

A1_large = zeros(nnodes,nnodes+nelements);
rhs1_large = zeros(nnodes,1);

iter=0;
count=0;
hess_ =  zeros(2*nelements,2*nelements);

%% Begin the solver

for step = 1:nsteps
    fprintf('Step = %d \n',step);
    fprintf('u_d = %f \n',u_d(step,1));
    while(~converge)
        epsilon_t(:,:) = epsilon_temp(:,:);
        sigma_t(:,:) = sigma_temp(:,:);

        u_pred = prev_disp + delta_t * prev_velo+(0.5-beta)*prev_acc * delta_t*delta_t;
        v_pred = prev_velo+(1-gamma)*prev_acc * delta_t;

        material(1:nelements,1) = epsilon_t(:,1);
        material(nelements+1:2*nelements,1) = sigma_t(:,1);

        ud(1,1) = 0.0;
        ud(2,1) = u_d(step,1);

        mechanical_st(1:nnodes,1) = displacement(:,1);
        mechanical_st(nnodes+1:nnodes+nelements,1) = sigma(:,1);

        [solution,fval,exitflag,lambda_] = call_minimizer(B_u,bfv,le,nelements,ud,mechanical_st,material,Mass_sigma,C,area);

        displacement = solution(1:nnodes);
        sigma = solution(nnodes+1:nnodes+nelements);
        epsilon = mtimes(B_u,displacement);

        acceleration = 1.0/beta/(delta_t*delta_t)*(displacement - u_pred);
        velocity = v_pred + gamma*delta_t*acceleration;

        % Find the closest material point
        for i = 1:nelements
            FE_point = [epsilon(i,1),sigma(i,1),0.5*sigma(i,1)*epsilon(i,1)];
            prev_state= [prev_epsilon_t(i,1),prev_sigma_t(i,1),prev_energy(i,1),prev_damage(i,1),prev_dissipation(i,1)];
            dam_mod = C*(1-prev_damage(i,1));
            % sig_1=0; eps_1=0; dam_1=0;
            % sig_2=0; eps_2=0; dam_2=0;

%             if(FE_point[0]<=delta_0):
%                 eps_1 = (C*C*FE_point[0] + dam_mod*FE_point[1])/(C*C+dam_mod*dam_mod)
%                 sig_1 = dam_mod*eps_1
%                 dam_1 = prev_damage[i,0]
%
%             if(FE_point[0]>=delta_0):
%                 eps_1 = (C*C*FE_point[0] + dam_mod*FE_point[1])/(C*C+dam_mod*dam_mod)
%                 sig_1 = (dam_mod*eps_1>prev_sigma_t[i,0])*1E6+(dam_mod*eps_1<prev_sigma_t[i,0])*dam_mod*eps_1
%                 dam_1 = prev_damage[i,0]

            eps_2=1E9; sig_2=1E9; dam_2=0.0;
            if(prev_damage(i,1)>0.0 || FE_point(1)>=delta_0)
                m=C*delta_0/(delta_0-delta_f); c=-delta_f*m;
                eps_2 = (C*C*FE_point(1)+m*FE_point(2)-m*c)/(m*m+C*C);
                if(eps_2<delta_0)
                    eps_2=delta_0;
                end
                sig_2 = m*eps_2+c;

                bl=((FE_point(1)>delta_0)*1.0)*((FE_point(1)<delta_f)*1.0);

                dam_2 = bl*delta_f*(FE_point(1)-delta_0)/FE_point(1)/(delta_f-delta_0)+(FE_point(1)>delta_f)*1.0;
                if(dam_2<elem_damage(i,step-1,1))
                    eps_2=elem_mat(i,step-1,1);sig_2=elem_mat(i,step-1,2); dam_2 = elem_damage(i,step-1,1);
                end         
            end

            eps_1 = (C*C*FE_point(1) + dam_mod*FE_point(2))/(C*C+dam_mod*dam_mod);
            sig_1 = dam_mod*eps_1;
            dam_1 = prev_damage(i,1);

            if(sig_1>sig_2)
                sig_1=1E9;
            end

            dist1 = C/2.0*(FE_point(1)-eps_1)*(FE_point(1)-eps_1);
            dist1 = dist1 + 1.0/(2.0*C)*(FE_point(2)-sig_1)*(FE_point(2)-sig_1);

            dist2 = C/2.0*(FE_point(1)-eps_2)*(FE_point(1)-eps_2);
            dist2 = dist2 + 1.0/(2.0*C)*(FE_point(2)-sig_2)*(FE_point(2)-sig_2);

            epsilon_temp(i,1) = (dist1>dist2)*eps_2+(dist1<=dist2)*eps_1;
            sigma_temp(i,1) = (dist1>dist2)*sig_2+(dist1<=dist2)*sig_1;
            damage(i,1) = (dist1>dist2)*dam_2+(dist1<=dist2)*dam_1;

%             data = generate_data(prev_state,1E2,500,1.0)
%             arg,temp = closest_node(FE_point,prev_state,data,C,le) #(1-damage[i,0])*
%             arg,temp = closest_node_graph(FE_point,prev_node[i,step-1,0],data,G,C,le)
%             epsilon_temp[i,0] = data[arg,0]
%             sigma_temp[i,0] = data[arg,1]
%             energy[i,0] = data[arg,2]
%             damage[i,0] = data[arg,3]
%             dissipation[i,0] = data[arg,4]
%             dist2[i,0]+=temp
%             node_temp[i,0] = arg
%         print(np.linalg.norm(sigma_t[:]-sigma_temp[:]))
        end
        if(norm(sigma_t-sigma_temp)<1E-4)

            prev_epsilon(:) = epsilon(:);
            prev_sigma(:) = sigma(:);

            prev_epsilon_t(:) = epsilon_temp(:);
            prev_sigma_t(:) = sigma_temp(:);

            prev_node(:,step,1) = node_temp(:,1);

            prev_energy(:) = energy(:);
            prev_damage(:) = damage(:);
            prev_dissipation(:) = dissipation(:);

            prev_disp(:) = displacement(:);
            prev_velo(:) = velocity(:);
            prev_acc(:) = acceleration(:);

            converge=true;
        end
    end

    elem_mech(:,step,1) = epsilon(:,1);
    elem_mech(:,step,2) = sigma(:,1);
    elem_mat(:,step,1) = epsilon_temp(:,1);
    elem_mat(:,step,2) = sigma_temp(:,1);
    elem_damage(:,step,1) = damage(:,1);
    elem_dissipation(:,step,1) = dissipation(:,1);
    node_displacement(:,step,1) = displacement(:,1);
    elem_lagrange(:,step,1) = lambda_.ineqlin;
    exit_flags(step,1) = exitflag;
    converge=false;
end
