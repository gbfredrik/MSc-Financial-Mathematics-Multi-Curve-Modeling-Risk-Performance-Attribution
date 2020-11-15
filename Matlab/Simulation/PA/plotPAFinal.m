function plotPAFinal(plotPAParams, tradeDatesAll_OOS)

    % fontsizes
    tickSize = 22;
    titleSize = 34;
    ylabelSize = 30;
    legendSize = 26;
    lineWidth = 2;
    
    %----------

    cNPV = plotPAParams{1}(2:end);
    totNPV =  plotPAParams{7};
    cSumRiskFactors = plotPAParams{2}(2:end);
    totSumRiskFactors = plotPAParams{8};
    cEps_A = plotPAParams{3}(2:end);
    totEps_A =  plotPAParams{9};
    cEps_I = plotPAParams{4}(2:end);
    totEps_I =  plotPAParams{10};
    cEps_P = plotPAParams{5}(2:end);
    totEps_P =  plotPAParams{11};
    cCarry = plotPAParams{6}(2:end);
    totCarry =  plotPAParams{12};
    resultDates = plotPAParams{13};
    
    cShift_1_f = plotPAParams{14}(2:end);
    cTwist_1_f = plotPAParams{15}(2:end);
    cButterfly_1_f = plotPAParams{16}(2:end);
    cFourth_sixth_1_f = plotPAParams{17}(2:end);
    cSum_second = plotPAParams{18}(2:end);
    
    cShift_1_pi = plotPAParams{19}(2:end);
    cTwist_1_pi = plotPAParams{20}(2:end);
    cButterfly_1_pi = plotPAParams{21}(2:end);
    cFourth_eight_1_pi = plotPAParams{22}(2:end);    
    
    totShift_1_f = plotPAParams{23};
    totTwist_1_f = plotPAParams{24};
    totButterfly_1_f = plotPAParams{25};
    totFourth_sixth_1_f = plotPAParams{26};
    totSum_second = plotPAParams{27};
    
    totShift_1_pi = plotPAParams{28};
    totTwist_1_pi = plotPAParams{29};
    totButterfly_1_pi = plotPAParams{30};
    totFourth_eight_1_pi = plotPAParams{31};   
    
    
    figure(1)
    x = 1:length(cNPV);
    % Carry
    plot(x, cCarry, 'Color', '#77AC30', 'LineWidth', lineWidth, 'HandleVisibility','off')
    hold on
    area(x, cCarry, 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'FaceColor', '#77AC30')
    
    
    % Tot risk factors   
    plot(x, cSumRiskFactors, 'Color', '#A2142F', 'LineWidth', lineWidth, 'HandleVisibility','off')
    area(x, cSumRiskFactors, 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'FaceColor', '#A2142F')       

    % eps_p
    plot(x, cEps_P, 'Color', '#7E2F8E', 'LineWidth', lineWidth, 'HandleVisibility','off')
    area(x, cEps_P, 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'FaceColor', '#7E2F8E')      
    
    % eps_a
    plot(x, cEps_A, 'Color', '#4DBEEE', 'LineWidth', lineWidth, 'HandleVisibility','off')
    area(x, cEps_A, 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'FaceColor', '#4DBEEE')  
    
    %eps_I
    plot(x, cEps_I, 'Color', '#D95319', 'LineWidth', 1.5, 'HandleVisibility','off')
    area(x, cEps_I, 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'FaceColor', '#D95319')   

    % NPV
    plot(x, cNPV, 'LineWidth', lineWidth+1, 'Color', [0 0 0])

    a = gca;
    a.FontSize = tickSize;
    set(a,'box','off','color','none')
    axes('Position',get(a,'Position'),'box','on','xtick',[],'ytick',[]);   
    axes(a)
    set(gca,'TickDir','out')
    set(gcf,'color','w');
    set(gca,'XLim',[0 length(totNPV)])    
    xticks([1 252 252*2 252*3 252*4 252*5 252*6 252*7])
    xticklabels({'2013', '2014', '2015', '2016', '2017', '2018', '2019', '2020'});  
    
    title('Accumulated performance attribution for all parameters', 'FontSize', titleSize, 'FontName', 'Times New Roman', 'FontWeight','Normal');
    legend('Carry', 'Risk Factors', '$\Delta \epsilon^P$', '$\Delta \epsilon^A$', '$\Delta \epsilon^I$', 'NPV', 'Interpreter', 'latex', 'FontSize', legendSize, 'Location','northwest', 'NumColumns', 1)
    ylabel('Accumulated attribution', 'FontName', 'Times New Roman', 'FontSize', ylabelSize)
  
    
    figure(2)
    % eps_p
    plot(x, cEps_P, 'Color', '#7E2F8E', 'LineWidth', lineWidth, 'HandleVisibility', 'off')
    hold on
    area(x, cEps_P, 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'FaceColor', '#7E2F8E')      
    
    % eps_a
    plot(x, cEps_A, 'Color', '#4DBEEE', 'LineWidth', lineWidth, 'HandleVisibility', 'off')
    area(x, cEps_A, 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'FaceColor', '#4DBEEE')  
    
    %eps_I
    plot(x, cEps_I, 'Color', '#D95319', 'LineWidth', lineWidth, 'HandleVisibility','off')
    area(x, cEps_I, 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'FaceColor', '#D95319')     
    
    a = gca;
    set(a,'box','off','color','none')
    axes('Position',get(a,'Position'),'box','on','xtick',[],'ytick',[]);   
    axes(a)
    a.FontSize = tickSize;      
    set(gca,'TickDir','out')
    set(gcf,'color','w');
    set(gca,'XLim',[0 length(totNPV)])    
    xticks([1 252 252*2 252*3 252*4 252*5 252*6 252*7])
    xticklabels({'2013', '2014', '2015', '2016', '2017', '2018', '2019', '2020'}); 
    title('Error terms', 'FontSize', titleSize, 'FontName', 'Times New Roman', 'FontWeight','Normal');
    legend('$\Delta \epsilon^P$', '$\Delta \epsilon^A$', '$\Delta \epsilon^I$', 'Interpreter', 'latex', 'FontSize', legendSize, 'Location','northwest')
    ylabel('Accumulated attribution', 'FontName', 'Times New Roman', 'FontSize', ylabelSize)
  
    
    figure(3)
    % shift_1_f
    plot(x, cShift_1_f, 'LineWidth', lineWidth, 'Color', [0 0 0])
    hold on 
    
    % twist_1_f
    plot(x, cTwist_1_f, 'LineWidth', lineWidth, 'Color', [0 0 1])
    
    % butterfly_1_f
    plot(x, cButterfly_1_f, 'LineWidth', lineWidth, 'Color', [0 1 1])  
    
    % fourth_sixth_1_f
    plot(x, cFourth_sixth_1_f, 'LineWidth', lineWidth, 'Color', [1 0 1])
      
    % shift_1_pi
    plot(x, cShift_1_pi, 'LineWidth', lineWidth, 'Color', [1 0 0]) 
    
    % twist_1_pi
    plot(x, cTwist_1_pi, 'LineWidth', lineWidth, 'Color', [0 1 0])
    
    % butterfly_1_pi
    plot(x, cButterfly_1_pi, 'LineWidth', lineWidth, 'Color', [.5 .5 .5])   
    
    % fourth_eight_1_pi
    plot(x, cFourth_eight_1_pi, 'LineWidth', lineWidth, 'Color', [.5 0 0])
    
    % sum_second
    plot(x, cSum_second, 'LineWidth', lineWidth, 'Color', [0 .5 0])     
        

    a = gca;
    set(a,'box','off','color','none')
    axes('Position',get(a,'Position'),'box','on','xtick',[],'ytick',[]);   
    axes(a)
    
    set(gca,'TickDir','out')
    set(gcf,'color','w');
    set(gca,'XLim',[0 length(totNPV)])    
    xticks([1 252 252*2 252*3 252*4 252*5 252*6 252*7])
    xticklabels({'2013', '2014', '2015', '2016', '2017', '2018', '2019', '2020'}); 
    a.FontSize = tickSize;
    
        
    title('Risk factors', 'FontSize', titleSize, 'FontName', 'Times New Roman', 'FontWeight','Normal');
    ylabel('Accumulated attribution', 'FontName', 'Times New Roman', 'FontSize', ylabelSize)
    legend('Shift$_{f}$', 'Twist$_{f}$', 'Butterfly$_{f}$', '4th-6th$_{f}$', 'Shift$_{\pi}$', 'Twist$_{\pi}$', 'Butterfly$_{\pi}$', '4th-8th$_{\pi}$', '$\Sigma$ 2:nd', 'FontSize', legendSize, 'Location','northwest', 'NumColumns', 3, 'Interpreter', 'latex')
        
    
end

