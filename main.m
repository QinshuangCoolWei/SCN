%% 3 node bruteforce vs 
scale_f = 10;
severity = 0; % (0,1]
demand_level = 1; % (0,1]
scale_rho = 1;

in_filename = '3node_simp_case.xlsx';
out_filename = sprintf('results_3node/3node_w_bf_fp_%d.mat',scale_f);
prepare_run(in_filename,out_filename,scale_f,severity,demand_level,scale_rho);

%%

% filename = 'food_data_2.xlsx';
% filename = '4node_simp_case.xlsx';
in_filename = 'simp_case.xlsx';

scale_f = 10;
severity = 0; % (0,1]
demand_level = 1; % (0,1]

% scale_rho = 0.9; %[0,1]
rho_list = 0;
for k = rho_list
    scale_rho = k/10;
    out_filename = sprintf('results_simp/simp_fp_%d_rho_%.1f.mat',scale_f*100,scale_rho);
    prepare_run(in_filename,out_filename,scale_f,severity,demand_level,scale_rho);
end



%%
in_filename = 'food_data_2.xlsx';

% scale severity
scale_f = 10;
% severity = 0; % (0,1]
demand_level = 1; % (0,1]
scale_rho = 1; %[0,1]

sev_list = 0:10;
for s = sev_list
    severity = s/10;
    out_filename = sprintf('results_food/food_fp_%d_sev_%.1f.mat',scale_f*100,severity);
    prepare_run(in_filename,out_filename,scale_f,severity,demand_level,scale_rho);
end



% scale demand
scale_f = 10;
severity = 0; % (0,1]
% demand_level = 1; % (0,1]
scale_rho = 1; %[0,1]

dl_list = 0:10;
for dl = dl_list
    demand_level = dl/10;
    out_filename = sprintf('results_food/food_fp_%d_demand_%.1f.mat',scale_f*100,demand_level);
    prepare_run(in_filename,out_filename,scale_f,severity,demand_level,scale_rho);
end


% scale rho
scale_f = 10;
severity = 0; % (0,1]
demand_level = 1; % (0,1]
% scale_rho = 0.9; %[0,1]

rho_list = 0:10;
for k = rho_list
    scale_rho = k/10;
    out_filename = sprintf('results_food/food_fp_%d_rho_%.1f.mat',scale_f*100,scale_rho);
    prepare_run(in_filename,out_filename,scale_f,severity,demand_level,scale_rho);
end