 function solution = concave_minimization(par_w,option)
% w = 1;      % weight for time

format long;
digits(100)
% scale = 1e-2;

scale = 1;

K = par_w.K;
V = par_w.V;
E = par_w.E;
C = par_w.C;

c_trans = par_w.c_trans;
c_prod = par_w.c_prod;
d_mat = par_w.d_mat;
a = par_w.a;
b = par_w.b;
c_0 = par_w.c_0*scale;
e_0 = par_w.e_0*scale;
I_0 = par_w.I_0;


dsrp_nodes = par_w.dsrp_nodes;
dsrp_prob = par_w.dsrp_prob;

N_K = length(K);
N_V = length(V);
N_E = length(E);
N_W = length(dsrp_prob);
N_C = length(C);
c_trans_ori = c_trans(1:N_E);
c_prod_ori = c_prod(1:N_V); 


% calculate some masks for calculation
prob_mask_V = zeros(N_V*N_W,1);    % prob masks
prob_mask_E = zeros(N_E*N_W,1);    % prob masks

for i_w = 1:N_W
    prob_mask_V((i_w-1)*N_V+1:i_w*N_V) = dsrp_prob(i_w);
    prob_mask_E((i_w-1)*N_E+1:i_w*N_E) = dsrp_prob(i_w);
end

Mask.prob_mask_V = prob_mask_V;
Mask.prob_mask_E = prob_mask_E;


one_K = ones(N_K,1);
one_W = ones(N_W,1);
one_V = ones(N_V,1);
one_E = ones(N_E,1);

% samples
N_p = N_V*N_W*N_K;
N_y = N_E*N_W*N_K;
N_x = N_V*N_W*N_K;
N_I = N_V*N_W*N_K;
N_alpha = N_E;
N_beta = N_V*N_K;
N_Delta = N_V*N_W*N_K;
% N_zeta = N_V*N_C;

num_phi = N_p+N_y+N_x+N_I+N_alpha+N_beta+N_Delta;

tau_0 = num_phi;
UB = 1e6;
LB = -1e6;
eps = 0.01;

% initialize S
S.tau = tau_0+1;        % the last one contains purely zeros
p_3d = zeros(N_V*N_W,N_K,S.tau);
y_3d = zeros(N_E*N_W,N_K,S.tau);
x_3d = zeros(N_V*N_W,N_K,S.tau);
I_3d = zeros(N_V*N_W,N_K,S.tau);
alpha_3d = zeros(N_E,1,S.tau);
beta_3d = zeros(N_V,N_K,S.tau);
Delta_3d = zeros(N_V*N_W,N_K,S.tau);
% zeta_3d = zeros(N_V,N_C,S.tau);

count_layer = 0;
% t = 1;
% initialize variables p,y,I,beta
for kk = 1:N_K
    if kk == 1  % set up once
%         p_slice = reshape(diag(num_phi*rep_W(c_prod_ori,N_W)),N_V*N_W,1,N_V*N_W); % create a 3d slice
        p_slice = reshape(diag(num_phi*c_prod),N_V*N_W,1,N_V*N_W); % create a 3d slice
        p_3d(:,kk,1:N_V*N_W) = p_slice;
        count_layer = count_layer+N_p;
%         y_slice = reshape(diag(num_phi*rep_W(c_trans_ori,N_W)),N_E*N_W,1,N_E*N_W); % create a 3d slice
        y_start = count_layer+1;
        y_slice = reshape(diag(num_phi*c_trans),N_E*N_W,1,N_E*N_W); % create a 3d slice
        y_3d(:,kk,count_layer+1:count_layer+N_E*N_W) = y_slice;
        count_layer = count_layer+N_y;

        I_start =  count_layer+1;
        I_slice = reshape(num_phi*(diag(sum(I_0,'all')+c_prod)),N_V*N_W,1,N_V*N_W); % create a 3d slice
        I_3d(:,kk,count_layer+1:count_layer+N_V*N_W) = I_slice;
        count_layer = count_layer+N_I;

        beta_start =  count_layer+1;
        beta_slice = reshape(eye(N_V)*num_phi,N_V,1,N_V); % create a 3d slice
        beta_3d(:,kk,count_layer+1:count_layer+N_V) = beta_slice;
        count_layer = count_layer+N_beta;

        x_start = count_layer+1;

    else
        p_3d(:,kk,(kk-1)*N_V*N_W+1:kk*N_V*N_W) = p_slice;
        y_3d(:,kk,y_start+(kk-1)*N_E*N_W:y_start+kk*N_E*N_W-1) = y_slice;
        I_3d(:,kk,I_start+(kk-1)*N_V*N_W:I_start+kk*N_V*N_W-1) = I_slice;
        beta_3d(:,kk,beta_start+(kk-1)*N_V:beta_start+kk*N_V-1) = beta_slice;
    end
end

Delta_start = x_start+N_x;
% initilize x,Delta
for jj = 1:N_V*N_W
    for kk = 1:N_K
        x_3d(jj,kk,x_start+(jj-1)*N_K+(kk-1)) = num_phi*d_mat(jj,kk);
        Delta_3d(jj,kk,Delta_start+(jj-1)*N_K+(kk-1)) = num_phi*d_mat(jj,kk);
    end
end    

alpha_start  = Delta_start+N_Delta;
alpha_3d(:,:,alpha_start:alpha_start+N_alpha-1) = reshape(num_phi*eye(N_E),N_E,1,N_E);



S.p = p_3d;
S.y = y_3d;
S.x = x_3d;
S.I = I_3d;
S.Delta = Delta_3d;
S.alpha = alpha_3d;
S.beta = beta_3d;

% compute initial phi_in (concave part)
eps_s = 1;
phi_in = zeros(S.tau,1);
for t = 1:N_p+N_y   % concave part only nonzero if p/y nonzero 
    p_temp = p_3d(:,:,t);
    y_temp = y_3d(:,:,t);
    phi_in(t) = sum(prob_mask_E'*(c_0.*(y_temp+eps_s).^a-c_0.*eps_s.^a),'all')+sum(prob_mask_V'*(e_0.*(p_temp+eps_s).^b-e_0.*eps_s.^b),'all');
end
% save('temp_save');
% load('temp_save','phi_in');

% compute bigM parameters
K_phi = sum(c_0.*a.*eps_s.^(a-1),'all') + sum(e_0.*b.*eps_s.^(b-1),'all');

M_2 = 1; 
p_max = max(S.p,[],'all');
y_max = max(S.y,[],'all');

M_1 = K_phi * sqrt(N_p*p_max^2 +N_y*y_max.^2)...
    +sum(c_0 .* (c_trans*one_K').^a ,'all') +sum(e_0 .* (c_prod*one_K').^b,'all');

% M_1 = K_phi * sqrt(N_p*p_max^2 +N_y*y_max.^2)...
%     +sum(c_0 .* (rep_W(c_trans_ori,N_W)*one_K').^a ,'all') +sum(e_0 .* (rep_W(c_prod_ori,N_W)*one_K').^b,'all');

% M_1 = K_phi * sqrt(sum(c_trans.^2,'all') +sum(c_prod.^2,'all'))...
%     +sum(c_0 .* (c_trans*one_K').^a ,'all') +sum(e_0 .* (c_prod*one_K').^b,'all');
UB_list = [];
LB_list = [];
obj_list = [];

oldSol.obj  = UB;
LB_old =LB;
max_iter = 10000;
tstart0 = tic;
% init_guess.p = 0;
% init_guess.y = 0;
% init_guess.alpha = 0;
% init_guess.beta = 0;
for count = 1:max_iter
    count
    UB-LB
    % solve MIP with current S samples
    if option == 1  % exclusive planning
        sol = exclusive_model(par_w,S,phi_in,Mask,M_1,M_2,scale);
    elseif option == 0  % plan with given beta (no w)
        sol = nw_solving(par_w,S,phi_in,Mask,M_1,M_2,scale);
    else
        tstart = tic;
        sol = model_solving_w(par_w,S,phi_in,Mask,M_1,M_2,scale);
        telapsed = toc(tstart);
    end
    % append
    p_temp = sol.p;
    y_temp = sol.y;


    % compute UB, LB
    f_temp = sol.obj-sol.xi;
    phi_hat_temp = phi_in'*sol.mu;
%     phi_temp = sum(prob_mask_E'*(c_0.*y1_temp.^a-c_0.*t_1.^(a-1).*y1_temp),'all')+sum(prob_mask_V'*(e_0.*p1_temp.^b-e_0.*t_3.^(b-1).*p1_temp),'all');
    phi_temp = sum(prob_mask_E'*(c_0.*(y_temp+eps_s).^a-c_0.*eps_s.^a),'all')+sum(prob_mask_V'*(e_0.*(p_temp+eps_s).^b-e_0.*eps_s.^b),'all');
    LB = f_temp + phi_hat_temp;
    UB = f_temp + phi_temp;
    UB_list = [UB_list;UB];
    LB_list = [LB_list;LB];
    obj_list = [obj_list;sol.obj];
    % compute phi_in
    phi_in = [phi_in;phi_temp];
    S.tau = S.tau+1;
    S.p = cat(3,S.p,p_temp);
    S.y = cat(3,S.y,y_temp);
    S.x = cat(3,S.x,sol.x);
    S.I = cat(3,S.I,sol.I);
    S.Delta = cat(3,S.Delta,sol.Delta);
    S.alpha = cat(3,S.alpha,sol.alpha);
    S.beta = cat(3,S.beta,sol.beta);
    if count ==1
        sol_min = sol;
        sol_min.obj = UB;
        init_guess.p = sol_min.p;
        init_guess.y = sol_min.y;
        init_guess.alpha = sol_min.alpha;
        init_guess.beta = sol_min.beta;
    else
        if UB<=sol_min.obj
            sol_min = sol;
            sol_min.obj = UB;
            init_guess.p = sol_min.p;
            init_guess.y = sol_min.y;
            init_guess.alpha = sol_min.alpha;
            init_guess.beta = sol_min.beta;
        end
    end
    

    if UB-LB < 2 | count>=30 %(mod(count,10)==0)
       break;    

    end
%     sol.obj - oldSol.obj
    
%     if (abs(sol.obj - oldSol.obj) <0.1)&&(count >=5)
%         break;
%     end
    LB_old = LB;

    oldSol.obj = sol.obj;
    oldSol.p = sol.p;
    oldSol.y = sol.y;
%     oldSol.p1 = sol.p1;
%     oldSol.y1 = sol.y1;
%     oldSol.p2 = sol.p2;
%     oldSol.y2 = sol.y2;
    oldSol.I = sol.I;
    oldSol.Delta = sol.Delta;
    oldSol.x = sol.x;
    oldSol.alpha = sol.alpha;
    oldSol.beta = sol.beta;
%     oldSol.zeta = sol.zeta;
    oldSol.mu = sol.mu;
    oldSol.nu_1 = sol.nu_1;
    oldSol.nu_2_p = sol.nu_2_p;
    oldSol.nu_2_y = sol.nu_2_y;
    oldSol.nu_3 = sol.nu_3;
    oldSol.nu_4 = sol.nu_4;
    oldSol.xi = sol.xi;

end
% sol.obj = UB;   % substitute with real objective value 
telapsed0 = toc(tstart0);
solution = sol_min;

end

% 
% function rep_mat = rep_W(mat,N_W)
% rep_mat = zeros(size(mat,1)*N_W, size(mat,2));
% for i_w=1:N_W
%     rep_mat((i_w-1)*size(mat,1)+1:i_w*size(mat,1),:) = mat;
% end
% end