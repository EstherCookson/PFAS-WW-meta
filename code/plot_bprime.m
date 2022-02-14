% Plots Temporal LME regression results 

% INPUTS
%   T = table of LME regression results and stats

function plot_bprime(T)

PFAS_names = {'PFHxA','PFHpA', 'PFOA',  'PFNA', 'PFDA','PFBS', 'PFHxS', 'PFOS'};

% configure country ordering
T.Country = cellstr(T.Country);
G_summary = groupsummary(T, ["Continent","Country"]);
T.Country = categorical(T.Country);
T.Country = reordercats(T.Country, G_summary.Country);


%%
f= figure();
f.Position= [1699 335 1201 426];
set(gcf,'color','w')
l = 0.01;
bottom = 0.12;
width = 0.1;
height = 0.7;

for i = 1:8
    t = T(T.PFAS == PFAS_names(i),:);
    left = (l*i) + (width*(i-1));
    axes('Position',[0.1+left bottom width height])
    plot(zeros(size(t,1)),t.Country, 'linewidth',1.5, 'Color', [0.5 0 0.5]);
    hold on

    for j = 1:length(t.b_prime)

        if t.b_prime(j) > 0
            f = t.b_prime(j)/max(T.b_prime);
            r = 0.5*(1+f);
            b = 0.5* (1-f);
            color = [r 0 b];
        
         else if t.b_prime(j)< 0
            f = t.b_prime(j)/ min(T.b_prime);
            r = 0.5* (1-f);
            b = 0.5*(1+f);
            color = [r 0 b];

         else if isnan(t.b_prime)
            color = 'k';
        end
        end 
        end
        errorbar(t.b_prime(j), t.Country(j), t.Upper(j) - t.b_prime(j), 'horizontal','o', 'MarkerFaceColor', color, 'Color',color, 'LineWidth',1.5,'MarkerSize',6.5);
        
        xlim([-2.2 2.2])
    end
    Intercept = round(t.Intercept(1),3);
    Slope = round(t.Slope(1),2);
    Slope_uncertainty = round(t.Upper_Slope(1) - t.Slope(1),2);
    title(PFAS_names(i))
    stitle = sprintf([ 'b = ',  num2str(Intercept), ', \nm = ' , num2str(Slope), char(177) , num2str(Slope_uncertainty)]);
    subtitle(stitle)
    hold off
 
    if i>1
        set(gca, 'yTickLabel', [])
    end
        xlabel("b'")

    set(gca, 'fontsize', 16);
    grid on
    box on;
    grid minor
    
end
Folder = cd;
Folder = fullfile(Folder, '..');
saveas(gcf,fullfile(Folder, '/figures and results/CountryMean.png'));

end
