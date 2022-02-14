
function GDPpercapRegression(GDP2019,C_obs2019, GDP, T_lmeResults, beta, Centered2019)

f= figure();
f.Position = [49 21 1291 776];
set(gcf,'color','w')
n_r = 2;
n_c = 4;
l = 0.001;
width = 0.85/n_c;
height = 0.73/n_r;
bottom = height*1.6;
PFAS_names = {'PFHxA','PFHpA', 'PFOA',  'PFNA', 'PFDA','PFBS', 'PFHxS', 'PFOS'};

for i = 1:8
    
    % perform linear regression on PFAS  and GDP per capita 
    lm = fitlm(GDP2019.(PFAS_names{i})/10000, C_obs2019.(PFAS_names{i}));
    
    lmCoef_p = lm.Coefficients.pValue;
    lmcoef = lm.Coefficients.Estimate;
    slope(i) = lmcoef(2);
    coef_CI= lm.coefCI;
    slope_uncertainty(i) = coef_CI(2,2) - slope(i);
    [y ci_y]= predict(lm);

    [ci_y I] = sort(ci_y);
    gdp_sorted  = GDP2019.(PFAS_names{i})/10000 ;
    gdp_sorted = gdp_sorted(I);
    Tgdp = T_lmeResults(find(T_lmeResults.PFAS == PFAS_names(i)),:);

    % scale marker size by n_obs
    [row col] = find(Tgdp{:, 'Country'} == GDP.CountryName');
    ms = Tgdp.n(row);
    ms = ((ms - mean(ms, 'omitnan'))/ std(ms,'omitnan'))*3.3 +10;
    
    % figure spacing
    if i < 5
        left = (l*i) + ((width+0.02)*(i-1));
    end
    if i > 4 
        bottom = (0.13);
        left = (l*(i-(n_c))) + ((width+0.02)*(i-(n_c+1)));
    end
    axes('Position',[0.06+left bottom width height])
    
    hold on
    for j = 1: size(ms)
        if ~isnan(Tgdp.b_prime(row(j)) + beta(i,1)+(Centered2019*beta(i,2)));
            p1 =plot(GDP.x2019(col(j))/10000, Tgdp.b_prime(row(j)) + beta(i,1)+(Centered2019*beta(i,2)),'ko', 'MarkerSize', ms(j), 'MarkerFaceColor', [0.5 0.5 0.5], 'HandleVisibility','off');
            if j == 1  
                p1 =plot(GDP.x2019(col(j))/10000, Tgdp.b_prime(row(j)) + beta(i,1)+(Centered2019*beta(i,2)),'ko', 'MarkerSize', ms(j), 'MarkerFaceColor', [0.5 0.5 0.5]);
            end
        end
    end

    p2 = plot(GDP2019.(PFAS_names{i})/10000,C_obs2019.(PFAS_names{i}), '.', 'Color', [0.6350 0.0780 0.1840]); 
    p3 = plot(GDP2019.(PFAS_names{i})/10000,y, 'k','linewidth',2);
    plot(gdp_sorted(:,1),ci_y(:,1), 'r--'); 
    plot(gdp_sorted(:,2),ci_y(:,2), 'r--');

    txt= ['slope = '  num2str(round(slope(i),2))  '\pm' num2str(round(slope_uncertainty(i) ,2))];
    text (3.5, -1.8, txt , 'FontSize', 14);
    txt =  sprintf('p = %1.2e',  lmCoef_p(2));
    text (3.5, -2.2, txt, 'FontSize', 14);
    text(3.5, -2.6,"R^2 = " + round(lm.Rsquared.Ordinary,2), 'FontSize',14);
    xlim([-0.1 6.9])
    ylim([-2.9 3.4])
    title(PFAS_names(i))

    if i >4
        xlabel('GDP per Capita [10,000 US$]')
    end

    if i == 1|i==(n_c+1)
        ylabel("log(\bfC_{eff}\rm)")
    end

    set(gca, 'fontsize', 18);
    set(gca, 'Color', 'w');
    box on;

    end  
lgnd= legend([p1 p2 p3], {"b+b' adjusted to 2019", 'Observations adjusted to 2019', 'GDP per Capita regression'});
lgnd.Position= [0.1027 0.0140 0.7483 0.0296] ;   
lgnd.Orientation = 'horizontal';

Folder = cd;
Folder = fullfile(Folder, '..');
saveas(gcf,fullfile(Folder, '/figures and results/GDPperCap.png'));

end