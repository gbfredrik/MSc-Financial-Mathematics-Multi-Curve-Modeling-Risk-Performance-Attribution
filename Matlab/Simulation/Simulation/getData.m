function [fAll_IS, piAll_IS, tradeDatesAll_IS, fAll_OOS, piAll_OOS, tradeDatesAll_OOS, ccy, fDatesAll] = getData(path_IS, path_OOS)

    load(path_IS);
    fAll_IS = fAll(:,1:3650);
    piAll_IS = piAll(:,1:3650);
    tradeDatesAll_IS = tradeDatesAll;
    clearvars -except fAll_IS piAll_IS tradeDatesAll_IS path_OOS
    %load('EUR_OOS_10YrCurves_Clean_Final.mat')
    load(path_OOS);
    %load('SEK_OOS_10YrCurves_Clean_Final.mat')
    fAll_OOS = fAll(:,1:3650);
    piAll_OOS = piAll(:,1:3650);
    fDatesAll = fDatesAll(:,1:3650);
    tradeDatesAll_OOS = tradeDatesAll;
    ccy = currency;
end
