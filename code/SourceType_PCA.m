clc; clear all; close all;
rng('default')

% import and structure data
Folder = cd;
Folder = fullfile(Folder, '..');
warning('off', 'MATLAB:table:ModifiedAndSavedVarnames')
filename = fullfile(Folder, '/data/SI_T2.xlsx');
data = readtable(filename,'Format','auto');
data.SourceType = categorical(data.SourceType);

PFAS_names = {'PFHxA','PFHpA', 'PFOA',  'PFNA', 'PFDA','PFBS', 'PFHxS', 'PFOS'};
PFAS_inf = {'PFHxA_inf','PFHpA_inf', 'PFOA_inf',  'PFNA_inf', 'PFDA_inf','PFBS_inf', 'PFHxS_inf', 'PFOS_inf'};
PFAS_eff = {'PFHxA_eff','PFHpA_eff', 'PFOA_eff',  'PFNA_eff', 'PFDA_eff','PFBS_eff', 'PFHxS_eff', 'PFOS_eff'};


for i = 1:8
    %convert influent values to double and log concentrations, remove ND
    data_og = data{:,PFAS_inf(i)};
    data(:,PFAS_inf(i)) = [];
    data{:,PFAS_inf(i)} = log10(cell_str_2_num(data_og));

    %convert effluent values to double and log concentrations, remove ND
    data_og = table2array(data(:,PFAS_eff(i)));
    data(:,PFAS_eff(i)) = [];
    data{:,PFAS_eff(i)} = log10(cell_str_2_num(data_og));

end

% remove data with no effluent observations
data(find(all(isnan(data{:,PFAS_eff}),2)),:) = [];
pca_data = data{:,PFAS_eff};

data.SourceType = categorical(data.SourceType);
data.SourceType = reordercats(data.SourceType, {'Municipal', 'Mixed', 'NA', ...
    'High Ind.'});

data.Month = month(datetime(data.Month, 'InputFormat','MMMM'));
data.Month(isnan(data.Month)) = 6;
data.Year = data.Year + (data.Month-1)/12 + 15/365;

var_label= PFAS_names;
n_v = length(var_label)*3;

%% Correlation 
[rho,pval]= corr(data{:,PFAS_eff},'rows','complete'); %,'Type','Spearman' );
rho_H = heatmap(rho); 
rho_H.YDisplayLabels = PFAS_names;
rho_H.XDisplayLabels =PFAS_names;
title('Correlations between PFAS in Effluent Data')
set(gcf,'color','w');


%% Perform PCA

b = 2;  % number of pc

% calculate z-scores
pca_data = pca_data  - nanmean(pca_data,1);
pca_data = pca_data ./ nanstd(pca_data,1);

% perform PCA
[coeff,score,pcvar,mu,v] = ppca(pca_data,b);


%% Plotting

% source
f_source = figure();
f_source.Position = [1441 696 620 355];

h2 = biplot(coeff(:,1:b),'Scores',score(:,1:b),'VarLabels',var_label);
g = data.SourceType; g= fillmissing(g,'constant','unknown');
groups = categories(removecats(g));
n_groups = length(groups);
colors = [ 0 0 1; 0 1 0; 0 1 0; 1 0 0];
markersymbol= {'o','*', '|', '^'};
for i = 1:n_groups
    indx = find(g == groups(i));
    h2_g(i) = h2(indx(1)+n_v);
    for j = 1:length(indx)
        h2(indx(j)+n_v).Marker = markersymbol(i);
        h2(indx(j)+n_v).MarkerEdgeColor = colors(i,:);
        h2(indx(j)+n_v).MarkerSize = 7;
        h2(indx(j)+n_v).MarkerFaceColor = colors(i,:);
        
    end
end
h2(22).HorizontalAlignment = 'right';
h2(23).HorizontalAlignment = 'right';
h2(23).VerticalAlignment = "top";

xlim([-0.6 0.61])
format short g
ax= gca;
ax.FontSize = 14;
xlabel("PC 1, variance = "+ round(pcvar(1),2, 'significant'), 'FontSize',16)
ylabel("PC 2, variance = "+ round(pcvar(2),2, 'significant'), 'FontSize',16)
legend(h2_g,groups)
set(gcf,'color','w')
set(gca,'color',[0.9 0.9 0.9])
title('PFAS Biplot and Source Type', 'FontSize', 14)

Folder = cd;
Folder = fullfile(Folder, '..');
saveas(gcf,fullfile(Folder, '/figures and results/SourceType.png'));


% % measurement year
% indx = find(isnan(data.Year));
% y_min = min(data.Year);
% y_max = max(data.Year);
% 
% figure()
% h2 = biplot(coeff(:,1:b),'Scores',score(:,1:b),'VarLabels',var_label);
% 
% for i = 1: size(score,1)
%     f = (data{i,'Year'} - y_min)/(y_max - y_min);
%     if any(indx(:)==i)
%         h2(i+n_v).MarkerSize = 1;
%     else
%     color = f.* [1 1 1];
%     h2(i+n_v).MarkerEdgeColor = color;
%     h2(i+n_v).MarkerSize = 15;
%     end
% end
% 
% set(gca,'color',[1 0.7 1])  % background color
% set(gcf,'color','w')
% colormap("gray")
% cb=colorbar;
% cb.Limits = [0 1];
% cb.TickLabels = round(linspace(y_min,y_max,11),1);
% title("PFAS Biplot and Sample Year")
% 
% 
