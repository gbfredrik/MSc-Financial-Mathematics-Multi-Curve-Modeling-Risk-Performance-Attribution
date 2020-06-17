measurementPath = 'X:\Examensarbete\Data\';
xSimulated = readmatrix(strcat(measurementPath,'KernelRnd1.csv'),'Delimiter',';');
xRealized = readmatrix(strcat(measurementPath,'KernelX.csv'),'Delimiter',';');

%%
f = kdeMulti(xSimulated, xRealized); % Maybe change input to kde
di = residual(f, f+0.0001);
[isBetter, probablility] = likelihoodRatioTest(di, 0.95);

