function plot_population(population)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       US population dynamics simulation between 2001 and 2011           %
%      (C) Michael Pokojovy                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    population(:, 2:3) = population(:, 2:3)*1E-6;

    set(gcf, 'PaperUnits', 'centimeters');
    xSize = 24; ySize = 12;
    xLeft = (21 - xSize)/2; yTop = (30 - ySize)/2;
    set(gcf,'PaperPosition', [xLeft yTop xSize ySize]);
    set(gcf,'Position',[0 0 xSize*50 ySize*50]);
    
    subplot_tight(1, 2, 1, [0.08 0.065]);
    barh(population(:, 1), population(:, 2), 'FaceColor', [0.0 0.5 1.0], 'EdgeColor', [0.0 0.5 1.0]);
    axis([0 2.5 0 100]);
    view(180, - 90);
    set(gca, 'YAxisLocation', 'right');
    xlabel('Number of male individuals in millions', 'interpreter', 'latex', 'FontSize', 18);
    
    subplot_tight(1, 2, 2, [0.08 0.065]);
    barh(population(:, 1), population(:, 3), 'FaceColor', [255 20 147]/255, 'EdgeColor', [255,20,147]/255);
    axis([0 2.5 0 100]);
    xlabel('Number of female individuals in millions', 'interpreter', 'latex', 'FontSize', 18);
    ylabel('Age in years', 'interpreter', 'latex', 'FontSize', 18);
end