function solution = model_solving_cont(par,S,phi_in,Mask,alpha,beta)
% solve the convex continuous problem for each initial point with given
% alpha and beta
format long;
digits(100)
% parameters
% eps = 1e-5;

% shrink the size
% scale = 1e-6;

K= par.K;
V = par.V;
E = par.E;
C = par.C;
VE_mat = par.VE_mat;
K2C = par.K2C;
c_prod = par.c_prod;
c_trans= par.c_trans;
I_0 = par.I_0;
d_mat = par.d_mat;
conv_mat = par.conv_mat;
f_prod = par.f_prod;
f_trans = par.f_trans;
rho = par.rho;
tau = S.tau;
prob_mask_V = Mask.prob_mask_V;
prob_mask_E = Mask.prob_mask_E;

ind_raw = par.ind_raw;
ind_no_prod = par.ind_no_prod;

% a = par.a;
% b = par.b;
% c_0 = par.c_0*scale;
% e_0 = par.e_0*scale;
% dsrp_nodes = par.dsrp_nodes;
dsrp_prob = par.dsrp_prob;



N_K = length(K);
N_V = length(V);
N_E = length(E);
N_C = length(C);
N_W = length(dsrp_prob);

% varimainables
p = sdpvar(N_V*N_W,N_K,'full');
Delta = sdpvar(N_V*N_W,N_K,'full');
I = sdpvar(N_V*N_W,N_K,'full');
x = sdpvar(N_V*N_W,N_K,'full');
y = sdpvar(N_E*N_W,N_K,'full');


mu = sdpvar(tau,1,'full');
nu_1 = sdpvar(1,1,'full');
nu_3 = sdpvar(tau,1,'full');
% nu_4 = binvar(tau,1,'full');
xi = sdpvar(1,1);
nu_2_p = sdpvar(N_V*N_W,N_K,'full');
nu_2_y = sdpvar(N_E*N_W,N_K,'full');
nu_2_x = sdpvar(N_V*N_W,N_K,'full');
nu_2_I = sdpvar(N_V*N_W,N_K,'full');
nu_2_Delta = sdpvar(N_V*N_W,N_K,'full');
% nu_2_alpha = sdpvar(N_E,1,'full');
% nu_2_beta = sdpvar(N_V,N_K,'full');

one_K = ones(N_K,1);
one_W = ones(N_W,1);
one_WK = zeros(N_W,N_V*N_W);
for i_w = 1:N_W
    one_WK(i_w,(i_w-1)*N_V+1:i_w*N_V)=1;
end

% mask_nu_p = repmat(eye(N_V*N_W),tau,1);
% mask_nu_y = repmat(eye(N_E*N_W),tau,1);
% 
constraints = [(VE_mat==-1)*y-(VE_mat==1)*y+p*conv_mat'-p+x-I_0+I == 0,...
    p*ones(N_K,1) <= c_prod,...
    Delta  <= d_mat,...
    Delta  >= d_mat-x,...
    Delta >= 0,...
    x >= 0,...
    y >= 0,...
    p >= 0,...
    I >= 0,...
    alpha <=1,...
    alpha >=0,...
    beta <=1,...
    beta >=0,...
%     zeta*ones(N_C,1) <= 1,...
    beta(:,ind_raw>0) == 0,...
    beta(ind_no_prod,:) == 0,...
    phi_in' * mu <= xi,...
    ones(1,tau)*mu == 1,...
    sum(S.y.*repmat(reshape(mu,1,1,tau),N_E*N_W,N_K),3) == y,...
    sum(S.p.*repmat(reshape(mu,1,1,tau),N_V*N_W,N_K),3) == p,...
    sum(S.x.*repmat(reshape(mu,1,1,tau),N_V*N_W,N_K),3) == x,...
    sum(S.I.*repmat(reshape(mu,1,1,tau),N_V*N_W,N_K),3) == I,...
    sum(S.Delta.*repmat(reshape(mu,1,1,tau),N_V*N_W,N_K),3) == Delta,...
    mu.*nu_3==0,...
%     nu_3 <= M_1*nu_4,...
%     mu <= M_2*(1-nu_4),...
    mu >= 0,...
    nu_3 >= 0,...
    nu_2_p(c_prod==0,:) == 0,...
];

for i_w = 1:N_W
    constraints = [constraints,...
%         p <= beta .* (c_prod * ones(N_K,1)'),...
%         y*ones(N_K,1) <= c_trans.*alpha,...
        p((i_w-1)*N_V+1:i_w*N_V,:) <= beta .* (c_prod((i_w-1)*N_V+1:i_w*N_V)*ones(1,N_K)),...
        y((i_w-1)*N_E+1:i_w*N_E,:)*one_K <= alpha .* c_trans((i_w-1)*N_E+1:i_w*N_E),...
%         y((i_w-1)*N_E+1:i_w*N_E,:)*one_K <=  c_trans((i_w-1)*N_E+1:i_w*N_E),...
        ];
end

% 
% for i_c = 1:N_C
% %     num_c = sum(K2C==i_c);
%     constraints = [constraints,...
%         zeta(:,i_c) >= sum(beta(:,K2C==i_c),2)/N_V
%         ];
% end

for t = 1:tau
    constraints = [constraints,...
%                 phi_in(t) - nu_1 - sum(nu_2_p.*S.p(:,:,t),"all")- sum(nu_2_y.*S.y(:,:,t),"all")+ nu_3(t) == 0,...
        phi_in(t) - nu_1 - sum(nu_2_p.*S.p(:,:,t)+nu_2_x.*S.x(:,:,t)...
        +nu_2_I.*S.I(:,:,t)+nu_2_Delta.*S.Delta(:,:,t),"all") - sum(nu_2_y.*S.y(:,:,t),"all") + nu_3(t) == 0
        ];


end


obj= prob_mask_V'* (rho.*Delta)*one_K + f_trans'*alpha + ones(1,N_V)*(f_prod.*beta)*one_K + xi;
% obj = E_penalty +f_trans_cost + f_trans_cost +xi;

% obj = (ones(1,N_W*N_V)*p1 + ones(1,N_W*N_E)*y1)*one_K;
% obj =  prob_mask_V'* (rho.*Delta)*one_K + xi;



options = sdpsettings('verbose',0,'solver','gurobi');
% options.gurobi.IntFeasTol = 1e-5;
% options.gurobi.TimeLimit = 1e3;
% options = sdpsettings('verbose',0,'solver','gurobi');

% provide initial guesses
% if init_guess.p~=0
%     assign(p,init_guess.p);
%     assign(y,init_guess.y);
%     assign(alpha,init_guess.alpha);
%     assign(beta,init_guess.beta);
% end

sol = optimize(constraints,obj,options);

% Analyze error flags
if sol.problem == 0
 % Extract and display value
%  solution.obj = value(obj-(ones(1,N_W*N_V)*p1 + ones(1,N_W*N_E)*y1)*one_K);
 solution.obj = value(obj);
 solution.p = value(p);
 solution.Delta = value(Delta);
 solution.x = value(x);
 solution.y = value(y);
 solution.alpha = value(alpha);
 solution.beta = value(beta);
%  solution.zeta = value(zeta);
 
 solution.I = value(I);
 solution.mu = value(mu);
 solution.nu_1 = value(nu_1);
 solution.nu_2_p = value(nu_2_p);
 solution.nu_2_y = value(nu_2_y);
 solution.nu_2_x = value(nu_2_x);
 solution.nu_2_I = value(nu_2_I);
 solution.nu_2_Delta = value(nu_2_Delta);
%  solution.nu_2_alpha = value(nu_2_alpha);
%  solution.nu_2_beta = value(nu_2_beta);
 solution.nu_3 = value(nu_3);
%  solution.nu_4 = value(nu_4);
 solution.xi = value(xi);
%  solution.p1 = value(p1);
%  solution.p2 = value(p2);
%  solution.y1 = value(y1);
%  solution.y2 = value(y2);

else
 display('Hmm, something went wrong!');
 sol.info
 yalmiperror(sol.problem);
end



end