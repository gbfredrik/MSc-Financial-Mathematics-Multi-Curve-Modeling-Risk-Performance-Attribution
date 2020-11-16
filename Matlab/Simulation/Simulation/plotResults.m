function plotResults(currCurve, xValues, results, tradeDatesAll_OOS, i, titleStr, figureNum)

    tickSize = 22;
    titleSize = 34;
    ylabelSize = 30;

    if figureNum == 1
        lineW = 2;
    else
        lineW = 2;
    end

    xTicksStr = {num2str(datestr(tradeDatesAll_OOS(i))), num2str(datestr(tradeDatesAll_OOS(i) + 365)), num2str(datestr(tradeDatesAll_OOS(i) + 730)) ...
        num2str(datestr(tradeDatesAll_OOS(i) + 1095)), num2str(datestr(tradeDatesAll_OOS(i) + 1460)), num2str(datestr(tradeDatesAll_OOS(i) + 1825)), ...
        num2str(datestr(tradeDatesAll_OOS(i) + 2190)), num2str(datestr(tradeDatesAll_OOS(i) + 2555)), num2str(datestr(tradeDatesAll_OOS(i) + 2920)) ...
        num2str(datestr(tradeDatesAll_OOS(i) + 3285)), num2str(datestr(tradeDatesAll_OOS(i) + 3650))};
    
    figure(figureNum)
    plot(xValues,results);
    hold on
    plot(xValues, currCurve, 'LineWidth', lineW, 'Color', 'k');
    hold off
    
    a = gca;
    a.FontSize = tickSize;
    xticks(0:365:3650)
    xticklabels(xTicksStr)
    xtickangle(15)    
    xlim([0 3650])


    set(a,'box','off','color','none')
    set(a,'TickDir','out')
    set(a,'XLim',[0 3650])
    b = axes('Position',get(a,'Position'),'box','on','xtick',[],'ytick',[]);
    axes(a)



    
    pause(0.1);
    hold off
    
    title(titleStr, 'FontSize', titleSize, 'FontName', 'Times New Roman', 'FontWeight','Normal');
    ylabel('Forward Rate', 'FontName', 'Times New Roman', 'FontSize', ylabelSize)



end