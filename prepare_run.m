function prepare_run(in_filename,out_filename,scale_f,severity,demand_level,scale_rho)
%% parameter
% format long;
% digits(100)


% special cases
% scale_f = 10;
% severity = 0; % (0,1]
% demand_level = 0.9; % (0,1]

% generate VE_mat data
% N_farm =  5;
% N_mid = 15;
% N_retail = 5;
% N_V = N_farm+N_mid+N_retail;
VE_mat_ori = readmatrix(in_filename,'Sheet','VE_mat');
VE_mat = VE_mat_ori(1:end-2,1:end-1);

N_V = size(VE_mat,1);
N_E = size(VE_mat,2);


% capacity for production
c_prod = VE_mat_ori(1:end-2,end);

% capacity for transportation
c_trans = VE_mat_ori(end-1,1:end-1);
c_trans = c_trans';
f_trans = VE_mat_ori(end,1:end-1);  % fixed setup cost for transportation
f_trans = f_trans';

% fixed setup cost for production
f_prod_ori = readmatrix(in_filename,'Sheet','f_prod');
f_prod = f_prod_ori(1:end-1,:);

% class of products
K2C = f_prod_ori(end,:);

% initial inventory
I_0 = readmatrix(in_filename,'Sheet','I_0');

% conversion rate
conv_mat = readmatrix(in_filename,'Sheet','conv_mat');

% demand
d_mat = readmatrix(in_filename,'Sheet','d_mat');

% unsatisfied penalty
rho = readmatrix(in_filename,'Sheet','rho');

% transportation costs

c_0 = readmatrix(in_filename,'Sheet','c_0');
a = readmatrix(in_filename,'Sheet','a');


e_0 = readmatrix(in_filename,'Sheet','e_0');
b = readmatrix(in_filename,'Sheet','b');

N_K = size(conv_mat,1);
N_C = max(K2C);

K = 1:N_K;
V = 1:N_V;
E = 1:N_E;
C = 1:N_C;  % classes of products




% f_prod = zeros(N_V,N_K);
b = ones(N_V,N_K)*0.9;
% e_0 = zeros(N_V,N_K);

% f_trans = zeros(N_E,1);
a = ones(N_E,N_K)*0.9;
% c_0 = zeros(N_E,N_K);
f_prod = f_prod*scale_f;
% f_trans = f_trans*10;

par.K = K;
par.V = V;
par.E = E;
par.C = C;
par.VE_mat = VE_mat;
par.K2C = K2C;
par.c_prod = c_prod;
par.c_trans = c_trans;
par.I_0 = I_0;
par.d_mat = d_mat*demand_level;
par.conv_mat = conv_mat;
par.f_prod = f_prod;
par.f_trans = f_trans;
par.rho = rho*scale_rho;


% %test


% par.t_1 = t_1;
% par.t_2 = t_2;
% par.t_3 = t_3;
% par.t_4 = t_4;
par.a = a;
par.b = b;
% par.c_0 = c_0;
% par.e_0 = e_0;


ind_raw = double(sum(conv_mat,1)==0)';   % identify raw materials (not made from anything)
% raw_provider = double(sum(VE_mat==1,2)==0)'; % identify raw material provider (not a destination from any link)
ind_no_prod = find(double(par.c_prod==0));

par.ind_raw = ind_raw;
par.ind_no_prod = ind_no_prod;


%% disruption
disruption = readmatrix(in_filename,'Sheet','disruption');
if size(disruption,1)==1
    dsrp_nodes = [zeros(1,N_V)];   % add "no-disruption" case
    dsrp_prob = 1;
else
    dsrp_nodes = [zeros(1,N_V);disruption(2:end,1:end-1)];   % add "no-disruption" case
    dsrp_prob = disruption(2:end,end);
    dsrp_prob = [1-sum(dsrp_prob);dsrp_prob];
end
N_W = length(dsrp_prob);
par_w.dsrp_nodes = dsrp_nodes;
par_w.dsrp_prob = dsrp_prob;


% original parameters (unchanged with disruptions)
par_w.K = K;
par_w.V = V;
par_w.E = E;
par_w.C = C;
par_w.K2C = K2C;
par_w.conv_mat = conv_mat;
par_w.ind_raw = ind_raw;
par_w.ind_no_prod = ind_no_prod;
par_w.f_prod = f_prod;
par_w.f_trans = f_trans;

% reshape parameters (vertically repeat some parameters N_W times)
par_w.d_mat = demand_level*rep_W(d_mat,N_W);
par_w.c_trans = rep_W(c_trans,N_W);
par_w.rho = rep_W(rho,N_W)*scale_rho;
par_w.a = rep_W(a,N_W);
par_w.b = rep_W(b,N_W);
par_w.c_0 = rep_W(c_0,N_W);
par_w.e_0 = rep_W(e_0,N_W);
par_w.I_0 = rep_W(I_0,N_W);

% t_3_w = rep_W(t_3,N_W);
b_w = rep_W(b,N_W);
c_prod_w = zeros(N_V*N_W,1);
VE_mat_w = zeros(N_V*N_W,N_E*N_W);
for i_w=1:N_W
    %     mask_temp = (dsrp_nodes(i_w,1:N_V)==0)';
    %     c_prod_w((i_w-1)*N_V+1:i_w*N_V) = c_prod.*mask_temp ;    %

    mask_temp = (dsrp_nodes(i_w,1:N_V)==1)'*severity+(dsrp_nodes(i_w,1:N_V)==0)';

    c_prod_w((i_w-1)*N_V+1:i_w*N_V) = c_prod.*mask_temp ;    %


    VE_mat_w((i_w-1)*N_V+1:i_w*N_V,(i_w-1)*N_E+1:i_w*N_E) = VE_mat; % VE matrix repeated on diagonal
end
par_w.c_prod = c_prod_w;
par_w.VE_mat = VE_mat_w;
% par_w.t_3 = t_3_w;
% par_w.b = b_w;

%% check
VE_mat_check = 1;
for ee = 1:N_E
    if length(find(VE_mat(:,ee))) ~= 2
        VE_mat_check = 0;
        break;
    end
end
G = VE2graph(VE_mat);   % SCN graph

% save('SCN_info_simp.mat')

% save('SCN_info_food.mat')

%% solving
% filename_old = sprintf('results_simp/simp_m2e_c_0.9_fp_%d.mat',0);
% solution = load(filename_old,'solution_m');
% solution_old = solution.solution_m;
% par_w.beta = solution_old.beta;
% par_w.alpha = solution_old.alpha;
% par_w.zeta = solution_old.zeta;

% solution_e = concave_minimization(par_w,1);   % exclusive
solution_m = concave_minimization(par_w,2);   % multi
save(out_filename)


% brute force method
% solution_bf = concave_minimization_bf(par_w,2,out_filename);   % multi


% filename = sprintf('results_food/nw_c_0.9_fp_%d.mat',scale_f*100);


% save(filename)

%
% filename_1 = sprintf('results_food/nw_c_0.9_fp_%d.mat',scale_f*100);
% solution = load(filename_1,'solution_m');
% solution = solution.solution_m;
% par_w.beta = solution.beta;
% par_w.alpha = solution.alpha;
%
% solution_m_nw2w = concave_minimization(par_w,0);   % multi, nw2w

% filename_2 = sprintf('results_food/nw2w_c_0.9_fp_%d.mat',scale_f*100);
% save(filename_2)




% filename_1 = sprintf('results/simp_sol_conc_0.9_fp_%d_noW.mat',scale_f*100);
% filename_1 = sprintf('results_food/nw_c_0.9_fp_%d.mat',scale_f*100);
% solution = load(filename_1,'solution_m');
% solution = solution.solution_m;
% par_w.beta = solution.beta;
% par_w.alpha = solution.alpha;
%
% solution_m_nw2w = concave_minimization(par_w,0);   % multi
% filename = sprintf('results/nw2w_c_0.5_fp_%d.mat',scale_f*100);
%
%
%
% filename = sprintf('results_food/c_0.9_fp_%d.mat',scale_f*100);
% save(filename)
%


end

%% Helper function
function rep_mat = rep_W(mat,N_W)
rep_mat = zeros(size(mat,1)*N_W, size(mat,2));
for i_w=1:N_W
    rep_mat((i_w-1)*size(mat,1)+1:i_w*size(mat,1),:) = mat;
end
end


function G = VE2graph(VE_mat)
% convert VE matrix to graph = (s,t)
N_E = size(VE_mat,2);
s = zeros(1,N_E);
t = zeros(1,N_E);
for e = 1:N_E
    s(e) = find(VE_mat(:,e)==-1);
    t(e) = find(VE_mat(:,e)==1);
end
G = digraph(s,t);
end