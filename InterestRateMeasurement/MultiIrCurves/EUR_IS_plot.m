close all
%% 1-day forward ON basis spread // T-Expected 1-day ON basis spread
figster = figure()
plot(datetime(tradeDatesAll, 'ConvertFrom', 'DateNum'), piAll(:,1))
datetick('x', 'yyyy')
xlabel('Trade date')
ylabel('Basis spread [%]')

exportgraphics(figster, 'Figures/EUR_IS_1day_Basis_Spread.png');



%% 

