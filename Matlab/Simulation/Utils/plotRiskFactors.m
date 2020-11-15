function plotRiskFactors(E, titlestr)
    tickSize = 22;
    titleSize = 34;
    ylabelSize = 30;
    lineW = 3;
    legendSize = 26;
    plot(1:3650, E, 'LineWidth', lineW);

    a = gca;
    a.FontSize = tickSize;
    xticks(0:365:3650)
    xlim([0 3650])
    xticklabels({0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10})   



    set(a,'box','off','color','none')
    set(a,'TickDir','out')
    set(a,'XLim',[0 3650])
    b = axes('Position',get(a,'Position'),'box','on','xtick',[],'ytick',[]);
    axes(a)



    
    pause(0.1);
    hold off
    
    title(titlestr, 'FontSize', titleSize, 'FontName', 'Times New Roman', 'FontWeight','Normal');
    legend("Shift", "Twist", "Butterfly", "4th", "5th", "6th", 'color','none', 'FontSize', legendSize, 'Location','northeast', 'NumColumns', 2, 'Interpreter', 'latex')
end

