%% network/general


%% severity diff

% with disruptions
obj_m_list = zeros(1,11);

scale_f = 10;
sev_list = 0:10;
for s = sev_list
    severity = s/10;
    filename = sprintf('results_simp/simp_m2e_fp_%d_sev_%.1f.mat',scale_f*100,severity);
    load(filename,'solution_m');


    obj_m_list(s+1) = solution_m.obj;
    

end


figure
set(gcf,'Units','Inches')
set(gcf,'Position',[4 4 4 3])
set(gca,'units','inches')
set(gcf, 'PaperUnits','inches');        
set(gcf, 'PaperSize', [4 3]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 4 3]);
plot(1-sev_list*0.1,obj_m_list,':k','Marker','*','Linewidth',2)


xlabel('Disruption Severity','Interpreter','latex','FontSize',12);
ylabel('$min$ $\mathcal J$','Interpreter','latex','FontSize',12)
% saveas(gcf,'simp_J_vs_severity','pdf')


%% obj vs demand
% with disruptions
obj_list = zeros(1,11);
dl_ori = 15;
scale_f = 10;
dl_list = 0:10;
for dl = dl_list
    demand_level = dl/10;
    filename = sprintf('results_simp/simp_fp_%d_demand_%.1f.mat',scale_f*100,demand_level);
    load(filename,'solution_m');
    obj_list(dl+1) = solution_m.obj;
   
end



figure
set(gcf,'Units','Inches')
set(gcf,'Position',[4 4 4 3])
set(gca,'units','inches')
set(gcf, 'PaperUnits','inches');        
set(gcf, 'PaperSize', [4 3]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 4 3]);
plot(dl_list*0.1*dl_ori,obj_list,':k','Marker','*','Linewidth',2)


xlabel('Demand','Interpreter','latex','FontSize',12);
ylabel('$min$ $\mathcal J$','Interpreter','latex','FontSize',12)
saveas(gcf,'simp_J_vs_demand','pdf')


%% obj vs rho
% with disruptions
rho_list1 = [0.01:0.01:0.09];
rho_list2 = [0.1:0.1:0.9];
rho_list = [0,rho_list1,rho_list2];

obj_list = zeros(1,length(rho_list));
rho_ori = 5000;

scale_f = 10;
filename = sprintf('results_simp/simp_fp_%d_rho_%.1f.mat',scale_f*100,0);
load(filename,'solution_m');
obj_list(1) = solution_m.obj;
for k = 1:length(rho_list1)
    scale_rho = rho_list1(k);
    filename = sprintf('results_simp/simp_fp_%d_rho_%.2f.mat',scale_f*100,scale_rho);
    load(filename,'solution_m');
    obj_list(k+1) = solution_m.obj;
end
for k = 1:length(rho_list2)
    scale_rho = rho_list2(k);
    filename = sprintf('results_simp/simp_fp_%d_rho_%.1f.mat',scale_f*100,scale_rho);
    load(filename,'solution_m');
    obj_list(k+1+length(rho_list1)) = solution_m.obj;
end

figure
set(gcf,'Units','Inches')
set(gcf,'Position',[4 4 4 3])
set(gca,'units','inches')
set(gcf, 'PaperUnits','inches');        
set(gcf, 'PaperSize', [4 3]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 4 3]);
plot(rho_list*rho_ori,obj_list,':k','Marker','*','Linewidth',2)


xlabel('$\rho$','Interpreter','latex','FontSize',12);
ylabel('$min$ $\mathcal J$','Interpreter','latex','FontSize',12)
saveas(gcf,'simp_J_vs_rho','pdf')

