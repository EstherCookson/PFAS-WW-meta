
%import and structure data
Folder = cd;
Folder = fullfile(Folder, '..');
warning('off', 'MATLAB:table:ModifiedAndSavedVarnames')
filename = fullfile(Folder, '/data/SI_T2.xlsx');
data = readtable(filename,'Format','auto');

PFAS_names = {'PFHxA','PFHpA', 'PFOA',  'PFNA', 'PFDA','PFBS', 'PFHxS', 'PFOS'};
PFAS_inf = {'PFHxA_inf','PFHpA_inf', 'PFOA_inf',  'PFNA_inf', 'PFDA_inf','PFBS_inf', 'PFHxS_inf', 'PFOS_inf'};
PFAS_eff = {'PFHxA_eff','PFHpA_eff', 'PFOA_eff',  'PFNA_eff', 'PFDA_eff','PFBS_eff', 'PFHxS_eff', 'PFOS_eff'};


for i = 1:8
    %convert influent values to double and log concentrations, remove ND
    data_og = data{:,PFAS_inf(i)};
    data(:,PFAS_inf(i)) = [];
    data{:,PFAS_inf(i)} = log10(ND2nan(data_og));

    %convert effluent values to double and log concentrations, remove ND
    data_og = table2array(data(:,PFAS_eff(i)));
    data(:,PFAS_eff(i)) = [];
    data{:,PFAS_eff(i)} = log10(ND2nan(data_og));

end

PFAS_InfEff= {'PFHxA_inf','PFHpA_inf', 'PFOA_inf',  'PFNA_inf', 'PFDA_inf','PFBS_inf', 'PFHxS_inf', 'PFOS_inf', 'PFHxA_eff','PFHpA_eff', 'PFOA_eff',  'PFNA_eff', 'PFDA_eff','PFBS_eff', 'PFHxS_eff', 'PFOS_eff'};
PFAS= data(:,PFAS_InfEff);


%% Analysis
f = figure();
f.Position = [2499 130 957 517];
l = 0.005;
width = 0.16;
height = 0.3;
bottom = height*2;

for i =1:8 
    % index influent and effluent
    int1= find(string(PFAS.Properties.VariableNames) == PFAS_inf(i));
    int2= find(string(PFAS.Properties.VariableNames) == PFAS_eff(i));
    pfas_inf= PFAS{:,int1};
    pfas_eff = PFAS{:,int2};

    % remove influent and effluent when either is nan
    ind= ((isnan(pfas_inf))| (isnan(pfas_eff)));
    pfas_inf(ind)=[];  
    pfas_eff(ind)=[];

    % determine average ratio between effluent and influent concentrations
    % 10^a = C_eff/C_inf -> a= log(C_eff) - log(C_inf)
    mdl = fitlm(pfas_inf, pfas_eff-pfas_inf,'constant');
    a(i) = mean(pfas_eff-pfas_inf); % = mdl.Coefficients.Estimate
    [rho(i) p_corr(i)] = corr(pfas_inf, pfas_eff);
    sigma(i) = std(pfas_eff-pfas_inf);
    pfas_eff_pred = pfas_inf+(a(i));
    
    x = -3.5:0.5:3.5;
    y= a(i)+x;

    % plotting
    if i < 6
        left = (l*i) + (width*(i-1));
    end
    if i > 5 
        bottom = (0.1);
        left = (l*(i-5)) + (width*(i-6));
    end
    axes('Position',[0.1+left bottom width height])
    x2 = [x, fliplr(x)];

    inBetween = [x, fliplr(repmat(3.5,1,15))];
    f = fill(x2, inBetween, [1 1 0.9]);
    hold on;
    inBetween = [x, fliplr(repmat(-3.5,1,15))];
    f = fill(x2, inBetween, [0.95 0.96 1]);

    plot(x,x, 'k', 'linewidth', 1.5)
    plot(pfas_inf, pfas_eff, 'b*'); hold on;
    plot(x,y,'Color',[0.8500, 0.3250, 0.0980], 'linewidth', 2);
    hold off;
  
    if i == 3|i==7
        xlabel('log(\bfC_{inf}\rm)')
    end
    if i == 1|i==6
        ylabel('log(\bfC_{eff}\rm)')
    else
        set(gca, 'yTickLabel', [])
    end

    title( PFAS_names(i) )
    subtitle("\langleC_{eff}/C_{inf}\rangle=" + round(10^a(i),2)); %+ ",  r="+ round(rho(i),2))
    xlim([-3.5 3.5])
    ylim([-3.5 3.5])
    set(gca, 'fontsize', 12);
    set(gca, 'Color', 'w');
    grid on;

end

lgnd = legend("$\mathrm{C_{eff} > C_{inf}}$", "$\mathrm{C_{eff} < C_{inf}}$","$\mathrm{C_{eff} = C_{inf}}$", "Observations", "$\mathrm{C_{eff}} = 10^a \times \mathrm{C_{inf}}$");
set(lgnd,'Interpreter','latex');
set(lgnd,'FontSize',14);
lgnd.Position= [0.6588 0.2190 0.1400 0.1286];

Folder = cd;
Folder = fullfile(Folder, '..');
saveas(gcf,fullfile(Folder, '/figures and results/Inf_Eff.png'));

%%

function output_cell = ND2nan(input_cell)

for i = 1:length(input_cell)
        
        str_val = input_cell{i};
        
        if contains(str_val,'<')
            str_val2 = nan; % *rand(1,1)
            input_cell{i}=str_val2;

        else
            input_cell{i}=str2double(str_val);
        end 
    end 

	output_cell = cell2mat(input_cell);

end
