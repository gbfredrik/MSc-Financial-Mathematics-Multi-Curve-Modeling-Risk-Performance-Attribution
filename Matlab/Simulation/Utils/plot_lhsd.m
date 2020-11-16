load('U_first')
U_first = U;
load('U_second')
load('ZZero')
load('rank_stat')
load('dfMZero')
load('dfCZero')
load('rho_zero')
tickSize = 12;
titleSize = 34;
ylabelSize = 30;
xlabelSize = 30;

figure(1)
subplot(1,2,1)
scatter(U_first(:,1),U_first(:,2), 70, 'filled');
a = gca;
a.FontSize = tickSize;
title("Original sample with Student's t-copula", 'FontSize', titleSize, 'FontName', 'Times New Roman', 'FontWeight','Normal');
xlabel('$U^1$', 'FontSize', xlabelSize, 'Interpreter', 'latex')
ylabel('$U^2$', 'FontSize', xlabelSize, 'Interpreter', 'latex')
axis([0 1 0 1])
pbaspect([1 1 1])
grid on

subplot(1,2,2)
scatter(U(:,1),U(:,2), 70, 'filled');
b = gca;
b.FontSize = tickSize;
title("LHSD sample", 'FontSize', titleSize, 'FontName', 'Times New Roman', 'FontWeight','Normal');
xlabel('$V^1$', 'FontSize', xlabelSize, 'Interpreter', 'latex')
ylabel('$V^2$', 'FontSize', xlabelSize, 'Interpreter', 'latex')
pbaspect([1 1 1])
grid on

%------
figure(2)
scatter(ZZero(:,1),ZZero(:,2), 70, 'filled');
c = gca;
c.FontSize = tickSize;
grid on
title("Inverse transformed sampled random variables", 'FontSize', titleSize, 'FontName', 'Times New Roman', 'FontWeight','Normal');
xlabel('$\epsilon_{0}^1$', 'FontSize', xlabelSize, 'Interpreter', 'latex')
ylabel('$\epsilon_{0}^2$', 'FontSize', xlabelSize, 'Interpreter', 'latex')
pbaspect([1 1 1])

corr(U_first)

corr(U)
