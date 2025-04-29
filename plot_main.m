%% network/general
load('SCN_info_simp.mat');
frame_width=3;
frame_height=2.5;

x_pos = [0,0,0,2,2,4,6];
y_pos = [4,2,0,3,1,2,2];


node_color = [0.4660 0.6740 0.1880];
edge_color = [0 0.4470 0.7410];
node_color_raw = '#8B8681'; 
node_color_retail = '#4169E1'; 



text_offset = 0.55*ones(1,N_V);
text_offset(5) = -0.55;


item_color = ["#FFF200",...  % 4 -> yellow
    "#7E2F8E",...    % 5 -> purple
    "#FFA500"];  % 6 -> orange




% graph products
figure
set(gcf,'Units','Inches')
set(gcf,'Position',[4 4 frame_width frame_height])
set(gca,'units','inches')
set(gcf, 'PaperUnits','inches');        
set(gcf, 'PaperSize', [frame_width frame_height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 frame_width frame_height]);

item_x = [0,0,0,2,2,4];
item_y = [4,2,0,3,1,2];

% convert conversion matrix to graph = (s,t)
% N_KE = size(conv_mat,1);
% s = zeros(1,N_KE);
% t = zeros(1,N_KE);
% for item = 1:N_KE
%     s(e) = find(conv_mat(:,item)>0);
%     t(e) = item;
% end
[s_item,t_item] = find(conv_mat);
G_item = digraph(s_item,t_item);

prod_N = plot(G_item, 'XData', item_x, 'YData', item_y, 'LineWidth',4,'EdgeColor',edge_color,'EdgeFontSize',12,'EdgeLabelColor',edge_color,  ...
    'NodeFontSize',12,'MarkerSize',24,'NodeColor',node_color_raw);
for item = 4:6
    highlight(prod_N,item,'NodeColor',item_color(item-3))
end


% Add new labels that are to the upper, right of the nodes
node_labels = cell(1,N_K);
for k = 1:N_K
    string =  sprintf('%d',k);
    node_labels{k} = string;
end
prod_N.NodeLabel =node_labels;

N_EK = length(s_item);
edge_labels = cell(1,N_EK);
for e = 1:N_EK
    string =  sprintf('%.1f',conv_mat(s_item(e),t_item(e)));
    edge_labels{e} = string;
end
prod_N.EdgeLabel =edge_labels;

text(prod_N.XData, prod_N.YData ,prod_N.NodeLabel, ...
    'VerticalAlignment','middle',...
    'HorizontalAlignment', 'center',...
    'FontSize', 12)
% Remove old labels
prod_N.NodeLabel = {}; 
axis off

% saveas(gcf,'prod_N_simp','pdf')


%% graph SCN 
figure
set(gcf,'Units','Inches')
set(gcf,'Position',[4 4 frame_width frame_height])
set(gca,'units','inches')
set(gcf, 'PaperUnits','inches');        
set(gcf, 'PaperSize', [frame_width frame_height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 frame_width frame_height]);



h = plot(G, 'XData', x_pos, 'YData', y_pos, 'LineWidth',4,'EdgeColor',edge_color,'EdgeFontSize',12,'EdgeLabelColor',edge_color,  ...
    'NodeFontSize',12,'MarkerSize',24,'NodeColor',node_color);

highlight(h,N_V,'NodeColor',node_color_retail)
for v = 1:3
    highlight(h,v,'NodeColor',node_color_raw)
end


% Add new labels that are to the upper, right of the nodes

node_labels = cell(1,N_V);
for v = 1:N_V
    string =  sprintf('%d,%d',par.V(v),par.c_prod(v));
    node_labels{v} = string;
end
h.NodeLabel =node_labels;

% edge_labels = cell(1,N_E);
% for e = 1:N_E
%     string =  sprintf('%d',par.E(e));
%     edge_labels{e} = string;
% end
% h.EdgeLabel =edge_labels;

text(h.XData, h.YData ,h.NodeLabel, ...
    'VerticalAlignment','middle',...
    'HorizontalAlignment', 'center',...
    'FontSize', 12)
% Remove old labels
h.NodeLabel = {}; 
axis off

% saveas(gcf,'SCN_simp','pdf')


%% obj diff

% with disruptions
obj_m_list = zeros(1,11);
obj_e_list = zeros(1,11);
sol_m_list = cell(11);
sol_e_list = cell(11);

k_list = 0:10;
for k = k_list
    filename = sprintf('results/simp_sol_conc_0.9_fp_%d.mat',k*100);
    load(filename,'solution_m');
    load(filename,'solution_e');
    sol_m_list{k+1} = solution_m;
    sol_e_list{k+1} = solution_e;

    obj_m_list(k+1) = solution_m.obj;
    obj_e_list(k+1) = solution_e.obj;
    

end


% without disruptions
obj_m_list_nw = zeros(1,11);
obj_e_list_nw = zeros(1,11);
sol_m_list_nw = cell(11);
sol_e_list_nw = cell(11);

k_list = 0:10;
for k = k_list
    filename_nw = sprintf('results/simp_sol_conc_0.9_fp_%d_noW.mat',k*100);
    load(filename_nw,'solution_m');
    load(filename_nw,'solution_e');
    sol_m_list_nw{k+1} = solution_m;
    sol_e_list_nw{k+1} = solution_e;

    obj_m_list_nw(k+1) = solution_m.obj;
    obj_e_list_nw(k+1) = solution_e.obj;
    

end


% nw2w
obj_m_list_nw2w = zeros(1,11);
sol_m_list_nw2w = cell(11);

sol_m_list_nw2w{1} = sol_m_list{1};
obj_m_list_nw2w(1) = obj_m_list(1);
k_list = 0:10;
for k = k_list(2:end)
    filename_nw2w = sprintf('results/nw2w_conc_0.9_fp_%d.mat',k*100);
    load(filename_nw2w,'solution_m_nw2w');
    sol_m_list_nw2w{k+1} = solution_m_nw2w;
    obj_m_list_nw2w(k+1) = solution_m_nw2w.obj;
end


% figure
% set(gcf,'Units','Inches')
% set(gcf,'Position',[4 4 4 3])
% set(gca,'units','inches')
% set(gcf, 'PaperUnits','inches');        
% set(gcf, 'PaperSize', [4 3]);
% set(gcf, 'PaperPositionMode', 'manual');
% set(gcf, 'PaperPosition', [0 0 4 3]);
% 
% % Minimal Expected Cost vs Fixed Production Cost
% plot(k_list*100,(obj_e_list_nw-obj_m_list_nw)./obj_e_list_nw,'-k','Marker','*','Linewidth',2)
% hold on 
% plot(k_list*100,(obj_e_list-obj_m_list)./obj_e_list,'-b','Marker','*','Linewidth',2)
% hold off
% % axis([min(k_list),max(k_list),min(obj_arr)-10,max(obj_arr)+100]);
% % title('Minimal Expected Cost vs Fixed Production Cost','FontSize',12);
% xlabel('$f_p$','Interpreter','latex','FontSize',12);
% ylabel('(J_E-J_M)/J_E')
% legend('w/o disruption', 'with disruption')
% saveas(gcf,'J_diff_vs_fp','pdf')
% 
%% obj 

figure
set(gcf,'Units','Inches')
set(gcf,'Position',[4 4 4 3])
set(gca,'units','inches')
set(gcf, 'PaperUnits','inches');        
set(gcf, 'PaperSize', [4 3]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 4 3]);
plot(k_list*100,obj_m_list_nw,':k','Marker','*','Linewidth',2)
hold on 
plot(k_list*100,obj_e_list_nw,':b','Marker','*','Linewidth',2)
hold on 
plot(k_list*100,obj_m_list,'-k','Marker','*','Linewidth',2)
hold on 
plot(k_list*100,obj_e_list,'-b','Marker','*','Linewidth',2)
hold off
leg=legend('multi-production w/o disruption','exclusive production w/o disruption','multi-production with disruption','exclusive production with disruption');
set(leg, 'Position', [0.4, 0.45, .25, .2])

xlabel('$f_p$','Interpreter','latex','FontSize',12);
ylabel('$min$ $\mathcal J$','Interpreter','latex','FontSize',12)
saveas(gcf,'simp_J_vs_fp','pdf')


% % 
% % linear without disruptions
% obj_m_list_lin_nw = zeros(1,11);
% obj_e_list_lin_nw = zeros(1,11);
% sol_m_list_lin_nw = cell(11);
% sol_e_list_lin_nw = cell(11);
% 
% k_list = 0:10;
% for k = k_list
%     filename_lin_nw = sprintf('results/simp_sol_conc_1_fp_%d_noW.mat',k*100);
%     load(filename_lin_nw,'solution_m');
%     load(filename_lin_nw,'solution_e');
%     sol_m_list_lin_nw{k+1} = solution_m;
%     sol_e_list_lin_nw{k+1} = solution_e;
% 
%     obj_m_list_lin_nw(k+1) = solution_m.obj;
%     obj_e_list_lin_nw(k+1) = solution_e.obj;
% end

% figure
% plot(k_list,obj_m_list_nw)
% hold on
% plot(k_list,obj_m_list_lin_nw)
% hold on
% plot(k_list,obj_m_list)
% hold off
% legend('nw_m','nw_m_l','w_m')
%% nw2w
figure
set(gcf,'Units','Inches')
set(gcf,'Position',[4 4 4 3])
set(gca,'units','inches')
set(gcf, 'PaperUnits','inches');        
set(gcf, 'PaperSize', [4 3]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 4 3]);


% Minimal Expected Cost vs Fixed Production Cost
plot(k_list*100,obj_m_list,'-k','Marker','*','Linewidth',2)
hold on 
plot(k_list*100,obj_m_list_nw2w,':b','Marker','*','Linewidth',2)
hold off
% axis([min(k_list),max(k_list),min(obj_arr)-10,max(obj_arr)+100]);
% title('Minimal Expected Cost vs Fixed Production Cost','FontSize',12);
xlabel('$f_p$','Interpreter','latex','FontSize',12);
ylabel('$min$ $\mathcal J$','Interpreter','latex','FontSize',12)
lgd = legend('risk-aware', 'risk-unaware');
lgd.Position = [0.35, 0.75, 0.15, 0.15]; % [left, bottom, width, height]
saveas(gcf,'simp_w_vs_nw','pdf')




%% fancy plots prepare

kk = 6;
% no disruption
p_m_nw = sol_m_list_nw{kk}.p;
p_e_nw = sol_e_list_nw{kk}.p;
beta_m_nw = sol_m_list_nw{kk}.beta;
beta_e_nw = sol_e_list_nw{kk}.beta;

% with disruption
p_m_w = sol_m_list{kk}.p;
p_e_w = sol_e_list{kk}.p;
p_m_nw2w = sol_m_list_nw2w{kk}.p;

beta_m_w = sol_m_list{kk}.beta;
beta_e_w = sol_e_list{kk}.beta;

% strong concave
filename = sprintf('results/simp_sol_conc_0.5_fp_%d.mat',1*100);
sol_m_temp = load(filename,'solution_m');
sol_m_temp = sol_m_temp.solution_m;
sol_e_temp = load(filename,'solution_e');


%% importance of capacity planning (multi vs exclusive) beta

% plot f = 100, with w, beta plots
temp_file = beta_e_w;
figure
set(gcf, 'Units', 'Inches');
set(gcf,'Position',[4 4 frame_width frame_height])

SCN = plot(G, 'XData', x_pos, 'YData', y_pos, 'LineWidth',4,'EdgeColor',edge_color,'EdgeFontSize',12,'EdgeLabelColor',edge_color,  ...
    'NodeFontSize',12,'MarkerSize',20,'NodeColor',node_color);

highlight(SCN,N_V,'NodeColor',node_color_retail)
for v = 1:3
    highlight(SCN,v,'NodeColor',node_color_raw)
end


for v = 4:6
    highlight(SCN,v,'MarkerSize',12)
end


node_labels = cell(1,N_V);
for v = 1:N_V
    string =  sprintf('%d',par.V(v));
    node_labels{v} = string;
end
SCN.NodeLabel =node_labels;

text(SCN.XData, SCN.YData+text_offset ,SCN.NodeLabel, ...
    'VerticalAlignment','middle',...
    'HorizontalAlignment', 'center',...
    'FontSize', 12)
% Remove old labels
SCN.NodeLabel = {}; 
axis off

ax_pos = get(gca, 'Position'); % Get axis position in figure coordinates
ax_xlim = get(gca, 'XLim');    % X-axis limits
ax_ylim = get(gca, 'YLim');    % Y-axis limits
hold on
for v = 4:6
    pie_r = 0.13;
    % Convert data coordinates to normalized figure coordinates
    x_norm = (SCN.XData(v) - ax_xlim(1)) / (ax_xlim(2) - ax_xlim(1)) * ax_pos(3) + ax_pos(1);
    y_norm = (SCN.YData(v) - ax_ylim(1)) / (ax_ylim(2) - ax_ylim(1)) * ax_pos(4) + ax_pos(2);
    
    % Create axes for pie at calculated normalized coordinates
    ax1 = axes('Position', [x_norm - pie_r/2, y_norm - pie_r/2, pie_r, pie_r]);
    pie_list = double(temp_file(v,4:end)>0); % enabled choices
    pie_num = size(pie_list);   % # enabled choices
    pie_labels = repmat({''},1,3);
%     pie_1 = pie([1,1,1],pie_labels);
    pie_1 = pie(pie_list,pie_labels);
    for pp = 1:length(pie_1)
        if isa(pie_1(pp), 'matlab.graphics.primitive.Patch')
        % Set the face color
            pie_1(pp).FaceColor = item_color((pp+1)/2);
        end       
    end
    hold on
end
hold on;
% plot legend
legend_labels = {'item 4', 'item 5', 'item 6'};
lgd_size = 20;
for i = 1:length(legend_labels)
    h(i) = plot(NaN, NaN, 's', 'MarkerFaceColor', item_color(i), 'MarkerEdgeColor', 'k', 'MarkerSize', lgd_size);
end
hold off;

% Add the legend to the main graph
lgd = legend(h, legend_labels);
lgd.Position = [0.65, 0.75, 0.1, 0.1]; % [left, bottom, width, height]
lgd.FontSize = 12;
hold off
set(gcf,'Position',[4 4 frame_width frame_height])
set(gca,'units','inches')
set(gcf, 'PaperUnits','inches');        
set(gcf, 'PaperSize', [frame_width frame_height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 frame_width frame_height]);

saveas(gcf,'beta_e_w_9_k_6','pdf')

%% importance of capacity planning (multi vs exclusive) p

% plot f = 100, with w, p plots
w_scenario = 4;
ind_list = (w_scenario-1)*N_V+1:w_scenario*N_V;
temp_file = p_m_w(ind_list,:);
figure
set(gcf, 'Units', 'Inches');
set(gcf,'Position',[4 4 frame_width frame_height])

SCN = plot(G, 'XData', x_pos, 'YData', y_pos, 'LineWidth',4,'EdgeColor',edge_color,'EdgeFontSize',12,'EdgeLabelColor',edge_color,  ...
    'NodeFontSize',12,'MarkerSize',20,'NodeColor',node_color);


highlight(SCN,N_V,'NodeColor',node_color_retail)
for v = 1:3
    highlight(SCN,v,'NodeColor',node_color_raw)
end

for v = 4:6
    highlight(SCN,v,'MarkerSize',12)
end

node_labels = cell(1,N_V);
for v = 1:N_V
    string =  sprintf('%d',par.V(v));
    node_labels{v} = string;
end
SCN.NodeLabel =node_labels;


text(SCN.XData, SCN.YData+text_offset ,SCN.NodeLabel, ...
    'VerticalAlignment','middle',...
    'HorizontalAlignment', 'center',...
    'FontSize', 12)
% Remove old labels
SCN.NodeLabel = {}; 
axis off
ax_pos = get(gca, 'Position'); % Get axis position in figure coordinates
ax_xlim = get(gca, 'XLim');    % X-axis limits
ax_ylim = get(gca, 'YLim');    % Y-axis limits
hold on
for v = 4:6
    if dsrp_nodes(w_scenario,v)==1
        highlight(SCN,v,'NodeColor','#B22222','MarkerSize',20)
    end
     
    if sum(temp_file(v,4:end)>0)==0
        continue;
    end
    pie_r = 0.13;
    % Convert data coordinates to normalized figure coordinates
    x_norm = (SCN.XData(v) - ax_xlim(1)) / (ax_xlim(2) - ax_xlim(1)) * ax_pos(3) + ax_pos(1);
    y_norm = (SCN.YData(v) - ax_ylim(1)) / (ax_ylim(2) - ax_ylim(1)) * ax_pos(4) + ax_pos(2);
    
    % Create axes for pie at calculated normalized coordinates
    ax1 = axes('Position', [x_norm - pie_r/2, y_norm - pie_r/2, pie_r, pie_r]);
    pie_list = double(temp_file(v,4:end)>0); % enabled choices
    pie_num = size(pie_list);   % # enabled choices
    pie_labels = repmat({''},1,3);
%     pie_1 = pie([1,1,1],pie_labels);
    pie_1 = pie(pie_list,pie_labels);
    for pp = 1:length(pie_1)
        if isa(pie_1(pp), 'matlab.graphics.primitive.Patch')
        % Set the face color
            pie_1(pp).FaceColor = item_color((pp+1)/2);
        end       
    end
    hold on
end
hold on;
% plot legend
legend_labels = {'item 4', 'item 5', 'item 6'};
lgd_size = 20;
for i = 1:length(legend_labels)
    h(i) = plot(NaN, NaN, 's', 'MarkerFaceColor', item_color(i), 'MarkerEdgeColor', 'k', 'MarkerSize', lgd_size);
end
hold off;

% Add the legend to the main graph
lgd = legend(h, legend_labels);
lgd.Position = [0.65, 0.75, 0.1, 0.1]; % [left, bottom, width, height]
lgd.FontSize = 12;
hold off
set(gcf,'Position',[4 4 frame_width frame_height])
set(gca,'units','inches')
set(gcf, 'PaperUnits','inches');        
set(gcf, 'PaperSize', [frame_width frame_height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 frame_width frame_height]);


saveas(gcf,'p_m_w_9_k_6','pdf')



%% importance of capacity planning (w vs nw2w) p

% plot f = 100, with w, p plots
w_scenario = 4;
ind_list = (w_scenario-1)*N_V+1:w_scenario*N_V;
temp_file = p_m_nw2w(ind_list,:);
figure
set(gcf, 'Units', 'Inches');
set(gcf,'Position',[4 4 frame_width frame_height])

SCN = plot(G, 'XData', x_pos, 'YData', y_pos, 'LineWidth',4,'EdgeColor',edge_color,'EdgeFontSize',12,'EdgeLabelColor',edge_color,  ...
    'NodeFontSize',12,'MarkerSize',20,'NodeColor',node_color);


highlight(SCN,N_V,'NodeColor',node_color_retail)
for v = 1:3
    highlight(SCN,v,'NodeColor',node_color_raw)
end

for v = 4:6
    highlight(SCN,v,'MarkerSize',12)
end

node_labels = cell(1,N_V);
for v = 1:N_V
    string =  sprintf('%d',par.V(v));
    node_labels{v} = string;
end
SCN.NodeLabel =node_labels;


text(SCN.XData, SCN.YData+text_offset ,SCN.NodeLabel, ...
    'VerticalAlignment','middle',...
    'HorizontalAlignment', 'center',...
    'FontSize', 12)
% Remove old labels
SCN.NodeLabel = {}; 
axis off
ax_pos = get(gca, 'Position'); % Get axis position in figure coordinates
ax_xlim = get(gca, 'XLim');    % X-axis limits
ax_ylim = get(gca, 'YLim');    % Y-axis limits
hold on
for v = 4:6
    if dsrp_nodes(w_scenario,v)==1
        highlight(SCN,v,'NodeColor','#B22222','MarkerSize',20)
    end
     
    if sum(temp_file(v,4:end)>0)==0
        continue;
    end
    pie_r = 0.13;
    % Convert data coordinates to normalized figure coordinates
    x_norm = (SCN.XData(v) - ax_xlim(1)) / (ax_xlim(2) - ax_xlim(1)) * ax_pos(3) + ax_pos(1);
    y_norm = (SCN.YData(v) - ax_ylim(1)) / (ax_ylim(2) - ax_ylim(1)) * ax_pos(4) + ax_pos(2);
    
    % Create axes for pie at calculated normalized coordinates
    ax1 = axes('Position', [x_norm - pie_r/2, y_norm - pie_r/2, pie_r, pie_r]);
    pie_list = double(temp_file(v,4:end)>0); % enabled choices
    pie_num = size(pie_list);   % # enabled choices
    pie_labels = repmat({''},1,3);
%     pie_1 = pie([1,1,1],pie_labels);
    pie_1 = pie(pie_list,pie_labels);
    for pp = 1:length(pie_1)
        if isa(pie_1(pp), 'matlab.graphics.primitive.Patch')
        % Set the face color
            pie_1(pp).FaceColor = item_color((pp+1)/2);
        end       
    end
    hold on
end
hold on;
% plot legend
legend_labels = {'item 4', 'item 5', 'item 6'};
lgd_size = 20;
for i = 1:length(legend_labels)
    h(i) = plot(NaN, NaN, 's', 'MarkerFaceColor', item_color(i), 'MarkerEdgeColor', 'k', 'MarkerSize', lgd_size);
end
hold off;

% Add the legend to the main graph
lgd = legend(h, legend_labels);
lgd.Position = [0.65, 0.75, 0.1, 0.1]; % [left, bottom, width, height]
lgd.FontSize = 12;
hold off
set(gcf,'Position',[4 4 frame_width frame_height])
set(gca,'units','inches')
set(gcf, 'PaperUnits','inches');        
set(gcf, 'PaperSize', [frame_width frame_height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 frame_width frame_height]);

saveas(gcf,'p_nw2w_9_k_6','pdf')


%% importance of capacity planning (exclusive) p

% plot f = 100, with w, p plots
w_scenario = 4;
ind_list = (w_scenario-1)*N_V+1:w_scenario*N_V;
temp_file = p_e_w(ind_list,:);
figure
set(gcf, 'Units', 'Inches');
set(gcf,'Position',[4 4 frame_width frame_height])

SCN = plot(G, 'XData', x_pos, 'YData', y_pos, 'LineWidth',4,'EdgeColor',edge_color,'EdgeFontSize',12,'EdgeLabelColor',edge_color,  ...
    'NodeFontSize',12,'MarkerSize',20,'NodeColor',node_color);


highlight(SCN,N_V,'NodeColor',node_color_retail)
for v = 1:3
    highlight(SCN,v,'NodeColor',node_color_raw)
end


node_labels = cell(1,N_V);
for v = 1:N_V
    string =  sprintf('%d',par.V(v));
    node_labels{v} = string;
end
SCN.NodeLabel =node_labels;


text(SCN.XData, SCN.YData+text_offset ,SCN.NodeLabel, ...
    'VerticalAlignment','middle',...
    'HorizontalAlignment', 'center',...
    'FontSize', 12)
% Remove old labels
SCN.NodeLabel = {}; 
axis off
ax_pos = get(gca, 'Position'); % Get axis position in figure coordinates
ax_xlim = get(gca, 'XLim');    % X-axis limits
ax_ylim = get(gca, 'YLim');    % Y-axis limits
hold on
for v = 4:6
    if dsrp_nodes(w_scenario,v)==1
        highlight(SCN,v,'NodeColor','#B22222','MarkerSize',20)
    end
     
    if sum(temp_file(v,4:end)>0)==0
        continue;
    end
    pie_r = 0.13;
    % Convert data coordinates to normalized figure coordinates
    x_norm = (SCN.XData(v) - ax_xlim(1)) / (ax_xlim(2) - ax_xlim(1)) * ax_pos(3) + ax_pos(1);
    y_norm = (SCN.YData(v) - ax_ylim(1)) / (ax_ylim(2) - ax_ylim(1)) * ax_pos(4) + ax_pos(2);
    
    % Create axes for pie at calculated normalized coordinates
    ax1 = axes('Position', [x_norm - pie_r/2, y_norm - pie_r/2, pie_r, pie_r]);
    pie_list = double(temp_file(v,4:end)>0); % enabled choices
    pie_num = size(pie_list);   % # enabled choices
    pie_labels = repmat({''},1,3);
%     pie_1 = pie([1,1,1],pie_labels);
    pie_1 = pie(pie_list,pie_labels);
    for pp = 1:length(pie_1)
        if isa(pie_1(pp), 'matlab.graphics.primitive.Patch')
        % Set the face color
            pie_1(pp).FaceColor = item_color((pp+1)/2);
        end       
    end
    hold on
end
hold on;
% plot legend
legend_labels = {'item 4', 'item 5', 'item 6'};
lgd_size = 20;
for i = 1:length(legend_labels)
    h(i) = plot(NaN, NaN, 's', 'MarkerFaceColor', item_color(i), 'MarkerEdgeColor', 'k', 'MarkerSize', lgd_size);
end
hold off;

% Add the legend to the main graph
lgd = legend(h, legend_labels);
lgd.Position = [0.65, 0.75, 0.1, 0.1]; % [left, bottom, width, height]
lgd.FontSize = 12;
hold off
set(gcf,'Position',[4 4 frame_width frame_height])
set(gca,'units','inches')
set(gcf, 'PaperUnits','inches');        
set(gcf, 'PaperSize', [frame_width frame_height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 frame_width frame_height]);


saveas(gcf,'p_e_w_9_k_6','pdf')



%% When f_p is large, multi-item choice will be bad (exclusive better) 
% simp with fixed beta,alpha for multi-item, compare exclusive and multi

figure
set(gcf,'Units','Inches')
set(gcf,'Position',[4 4 4 3])
set(gca,'units','inches')
set(gcf, 'PaperUnits','inches');        
set(gcf, 'PaperSize', [4 3]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 4 3]);


% with disruptions
k_list = 0:3:39;
len = length(k_list);
obj_m_list = zeros(1,len);
obj_e_list = zeros(1,len);
sol_m_list = cell(len);
sol_e_list = cell(len);

for k = 1:len
    filename = sprintf('results_simp/simp_m2e_c_0.9_fp_%d.mat',(k-1)*3*100);
    load(filename,'solution_m');
    load(filename,'solution_e');
    sol_m_list{k} = solution_m;
    sol_e_list{k} = solution_e;

    obj_m_list(k) = solution_m.obj;
    obj_e_list(k) = solution_e.obj;
    
end


% Minimal Expected Cost vs Fixed Production Cost
plot(k_list*100,obj_e_list,'-k','Marker','*','Linewidth',2)
hold on 
plot(k_list*100,obj_m_list,':b','Marker','*','Linewidth',2)
hold off
axis([min(k_list*100),max(k_list*100)+100,min(obj_m_list)-2000,max(obj_m_list)+1000]);
% title('Minimal Expected Cost vs Fixed Production Cost','FontSize',12);
xlabel('$f_p$','Interpreter','latex','FontSize',12);
ylabel('$min$ $\mathcal J$','Interpreter','latex','FontSize',12)
lgd = legend('exclusive production', 'multi-item production');
lgd.Position = [0.5, 0.2, 0.15, 0.15]; % [left, bottom, width, height]
saveas(gcf,'simp_m_vs_e_fix','pdf')
