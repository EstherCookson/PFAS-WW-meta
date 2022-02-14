% Global patterns and Temporal Trends of PFAS in Wastewater:  
% A Meta-analysis 
%
%~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

clc; clear all; close all; 

%% Preprocess data 

%import and structure data
Folder = cd;
Folder = fullfile(Folder, '..');
warning('off', 'MATLAB:table:ModifiedAndSavedVarnames')
filename = fullfile(Folder, '/data/SI_T1.xlsx');
data = readtable(filename,'Format','auto');

% remove industrial data and data without known measurement year
data.SourceType = categorical(data.SourceType);
ind_indx = data.SourceType == 'High Ind.';
data(ind_indx,:)=[];
data(isnan(data.Year),:)=[];

data.Author = categorical(data.Author);
data.Continent = categorical(data.Continent);
data.Country = categorical(data.Country);

PFAS_names = {'PFHxA','PFHpA', 'PFOA',  'PFNA', 'PFDA','PFBS', 'PFHxS', 'PFOS'};
PFAS_inf = {'PFHxA_inf','PFHpA_inf', 'PFOA_inf',  'PFNA_inf', 'PFDA_inf','PFBS_inf', 'PFHxS_inf', 'PFOS_inf'};
PFAS_eff = {'PFHxA_eff','PFHpA_eff', 'PFOA_eff',  'PFNA_eff', 'PFDA_eff','PFBS_eff', 'PFHxS_eff', 'PFOS_eff'};

% convert concentrations to log concentrations (set ND=0.5*LOD)
for i = 1:8
    data_og = data{:,PFAS_inf(i)};
    data(:,PFAS_inf(i)) = [];
    data{:,PFAS_inf(i)} = log10(cell_str_2_num(data_og));

    data_og = data{:,PFAS_eff(i)};
    data(:,PFAS_eff(i)) = [];
    data{:,PFAS_eff(i)} = log10(cell_str_2_num(data_og));
end

% define sample date by month and year
data.Month = month(datetime(data.Month, 'InputFormat','MMMM'));
data.Month(isnan(data.Month)) = 6;
data.Year = data.Year + (data.Month-1)/12 + 15/365;
mean_yr_all = nanmean(data.Year);
data.CenteredYear = data.Year- mean_yr_all;

% sort countries by continent
G_summary = groupsummary(data, ["Continent","Country"]);
all_countries = G_summary.Country;
n_C = length(categories(data.Country));

% add GDP data
Folder = cd;
Folder = fullfile(Folder, '..');
warning('off', 'MATLAB:table:ModifiedAndSavedVarnames')
filename = fullfile(Folder, '/data/GDP_per_capita.xlsx');
GDP = readtable(filename,'Format','auto', 'Range','A4:E40' );

GDP.CountryName = categorical(GDP.CountryName);
Logical = data.Country == GDP.CountryName';
data.GDP = nan(size(data,1),1);
for i = 1: size(GDP,1)
    data.GDP(Logical(:,i)) = GDP.x2019(i);
end
C_hat2019 = struct;  
GDP2019 = struct; 
Centered2019 = (2019- mean_yr_all);

% Initialize table for export of regression results
T_lmeResults = table();
A = repmat(PFAS_names,n_C,1);
T_lmeResults.PFAS = A(:); 
FixedEffectStats = {'Intercept', 'Slope', 'SE_Slope'};
T_lmeResults{:, FixedEffectStats} = nan(n_C*8,3);
A = repmat(G_summary.Continent,8,1); 
T_lmeResults.Continent = A(:);
A = repmat(all_countries,8,1);
T_lmeResults.Country = A(:); 
CountryEffectStats = {'n_p', 'n', 'b_prime', 'p_Int', 'SE_Int', 'Lower','Upper'};
T_lmeResults{:, CountryEffectStats} = nan(n_C*8,7);


%% Temporal Regression Analysis

f = figure();
set(gcf,'color','w')

for i = 1:8
    
    % organize data specific to the i'th PFAS
    T = [data(:, PFAS_eff(i)) data(:, {'Country', 'Year', 'CenteredYear', 'Continent','Author', 'GDP'})];
    T((ismissing(data(:,PFAS_eff(i)))) | isinf(data{:,PFAS_eff(i)}) ,:) =[];
    g = unique(T.Country);
    g2 = unique(T.Continent);
    Tg_summary = groupsummary(T, ["Continent","Country"]);
    Tg_summary.PFAS(:) = PFAS_names(i);

    %%%%%% LM %%%%%%%%%%%
    lm = fitlm(T.CenteredYear, T{:,PFAS_eff(i)});

    lm_c_CI = lm.coefCI;
    lm_slope(i) = lm.Coefficients.Estimate(2);
    lm_slopeCI(i,:) = lm_c_CI(2,:);
    lm_intCI(i,:) = lm_c_CI(1,:);

    %%%%%% LME with Random Intercept %%%%%%%%%%%
    formula = append(string(PFAS_eff(i)), "~CenteredYear+1 +(1|Country)");
    lme = fitlme(T,formula);
    
    % fixed effect estimates
    beta(i,:) = fixedEffects(lme);
    param_CI= lme.coefCI;
    lme_slopeCI(i,:) = param_CI(2,:);
    [~,~,FE]= lme.fixedEffects;

    % country effect estimates (b')
    [~,~,STATS] = randomEffects(lme); 
    STATS.Level = nominal(STATS.Level);
   
    % plot observations and regression predictions
    plotTemporalRegression(i, PFAS_eff, T, G_summary, lm, beta, PFAS_names);

    %%%%%%%%% Export results to Table %%%%%
    % Add fixed effects to export table
    T_lmeResults.Intercept((1:n_C) + ((i-1)*n_C))= repmat(beta(i,1), n_C,1);
    T_lmeResults.Slope((1:n_C) + ((i-1)*n_C)) = repmat(beta(i,2), n_C,1);
    T_lmeResults.p_Slope((1:n_C) + ((i-1)*n_C)) = repmat(FE.pValue(2), n_C,1);
    T_lmeResults.SE_Slope((1:n_C) + ((i-1)*n_C)) =  repmat(FE.SE(2), n_C,1);
    T_lmeResults.Low_Slope((1:n_C) + ((i-1)*n_C)) = repmat(param_CI(2,1), n_C,1);
    T_lmeResults.Upper_Slope((1:n_C) + ((i-1)*n_C)) = repmat(param_CI(2,2), n_C,1);
    
    % Add country specific Effect / Descritpions to table
    [~,idx_C] = ismember(Tg_summary.Country,all_countries,'rows');
    T_lmeResults.n((idx_C)+ ((i-1)*n_C)) = Tg_summary.GroupCount;
    Tg_summary = groupsummary(T, ["Continent","Country","Author"]);
    n_p = sum(Tg_summary.Country == all_countries');
    T_lmeResults.n_p((1:n_C) + ((i-1)*n_C))= n_p';
    for j = 1: n_C
        if ~isempty(STATS.Estimate(STATS.Level == all_countries(j)))
            T_lmeResults.b_prime(j+((i-1)*n_C)) = STATS.Estimate(STATS.Level == all_countries(j));
            T_lmeResults.Lower(j+((i-1)*n_C)) = STATS.Lower(STATS.Level == all_countries(j));
            T_lmeResults.Upper(j+((i-1)*n_C)) = STATS.Upper(STATS.Level == all_countries(j));
            T_lmeResults.p_Int(j+((i-1)*n_C)) = STATS.pValue(STATS.Level == all_countries(j));
            T_lmeResults.SE_Int(j+((i-1)*n_C)) = STATS.SEPred(STATS.Level == all_countries(j));
        end
    end

    % lme predicted observations
    C_hat = beta(i,1) + beta(i,2) * T.CenteredYear + ((T.Country==g')* STATS.Estimate);

    % adjust observations to 2019 predicted concentration
    C_obs2019.(PFAS_names{i}) = T{:,PFAS_eff(i)} - ((T.CenteredYear - Centered2019)* beta(i,2)); 
    C_obs2019.Country{i} = T.Country;
    GDP2019.(PFAS_names{i}) = T.GDP; 

 clear STATS; clear STATS2; clear T; clear x; clear y;
end

% percent change of concentrations each year
percent_change = (10.^(beta(:,2))-1)*100; 

Folder = cd;
Folder = fullfile(Folder, '..');
saveas(gcf,fullfile(Folder, '/figures and results/TemporalRegression.png'));
%% Export LME Results Table 

T_lmeResults.Total_Estimate = T_lmeResults.b_prime + T_lmeResults.Intercept;
Folder = cd;
Folder = fullfile(Folder, '..');
writetable(T_lmeResults,fullfile(Folder, '/figures and results/LME_stat_results.xlsx'));


%% PFAS and GDP per capita 
% perform regression on observeded log(C_eff) adjusted to 2019 (C_obs2019) with respect
% to 2019 national GDP per capita [US$] (GDP2019)

T_lmeResults.Country = categorical(T_lmeResults.Country);
T_lmeResults.PFAS = categorical(T_lmeResults.PFAS);


GDPpercapRegression(GDP2019, C_obs2019, GDP, T_lmeResults, beta, Centered2019)



%% Plot b' : mean estimates for each country
plot_bprime(T_lmeResults)

