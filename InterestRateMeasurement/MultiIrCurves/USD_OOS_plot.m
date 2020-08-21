% Plots for USD curves
close all;
fTauAll = fAll + piAll;

indices = 1:1;

fig = figure();
hold on

plot(fAll(indices,1:3655)', 'b')
plot(fTauAll(indices,1:3655)', 'g')
legend ('ON', '3M', 'Location', 'SouthEast')
xlabel('Days')
ylabel('Forward rate, %')

hold off

fprintf('Graphed trade dates: \n')
disp(datestr(tradeDatesAll(1, indices)))
exportgraphics(fig, 'Figures/USD_OS_1.png');
