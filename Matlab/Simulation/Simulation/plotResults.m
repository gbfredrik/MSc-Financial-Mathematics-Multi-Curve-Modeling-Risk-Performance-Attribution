function plotResults(currCurve, xValues, results, tradeDatesAll_OOS, i, titleStr, figureNum)

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
    title(titleStr, 'FontSize', 14);
    a = gca;
    set(a,'box','off','color','none')
    axes('Position',get(a,'Position'),'box','on','xtick',[],'ytick',[]);
    axes(a)
    %linkaxes([a b])
    ylabel('ForwardRate', 'FontSize', 14) 
    xticks(0:365:3650)
    set(gca,'TickDir','out')
    set(gcf,'color','w');
    set(gca,'XLim',[0 3650])
    xticklabels(xTicksStr);
    xtickangle(45)
    pause(0.1);
    hold off





end