clear;
dataPath = 'Curves/';

% Set file name
targetFileNames = 'EUR_IS_10Yr';
% Options:
    % 'EUR_IS_10Yr'
    % 'EUR_OOS_10Yr'
    % 'USD_OOS_10Yr'
    % 'SEK_OOS_10Yr'
    
sourceFileName = 'EUR_IS_10YrCurves_Clean_Final2.mat';
% Options:
    % 'EUR_IS_10YrCurves.mat'
    % 'EUR_OOS_10YrCurves_Clean.mat'
    % 'USD_OOS_10YrCurves_Clean.mat'
    % TODO 'SEK_OOS_10YrCurves.mat'

load(strcat(dataPath, sourceFileName))

%% tradeDatesAll
dlmwrite(strcat(dataPath, targetFileNames, '_tradeDatesAll.csv'), size(tradeDatesAll, 2), 'delimiter',';');
dlmwrite(strcat(dataPath, targetFileNames, '_tradeDatesAll.csv'), tradeDatesAll, '-append', 'delimiter', ';', 'precision', '%.0f');

%% fDatesAll
dlmwrite(strcat(dataPath, targetFileNames, '_fDatesAll.csv'), size(fDatesAll), 'delimiter',';');
dlmwrite(strcat(dataPath, targetFileNames, '_fDatesAll.csv'), fDatesAll, '-append', 'delimiter', ';');

%% TAll
dlmwrite(strcat(dataPath, targetFileNames, '_TAll.csv'), size(TAll), 'delimiter',';');
dlmwrite(strcat(dataPath, targetFileNames, '_TAll.csv'), TAll, '-append', 'delimiter', ';', 'precision', 16);

%% fAll
dlmwrite(strcat(dataPath, targetFileNames, '_fAll.csv'), size(fAll), 'delimiter',';');
dlmwrite(strcat(dataPath, targetFileNames, '_fAll.csv'), fAll, '-append', 'delimiter', ';', 'precision', 16); 

%% piAll
dlmwrite(strcat(dataPath, targetFileNames, '_piAll.csv'), size(piAll), 'delimiter',';');
dlmwrite(strcat(dataPath, targetFileNames, '_piAll.csv'), piAll, '-append', 'delimiter', ';', 'precision', 16);

%% zAll
dlmwrite(strcat(dataPath, targetFileNames, '_zAll.csv'), size(zAll), 'delimiter',';');
dlmwrite(strcat(dataPath, targetFileNames, '_zAll.csv'), zAll, '-append', 'delimiter', ';', 'precision', 16);
