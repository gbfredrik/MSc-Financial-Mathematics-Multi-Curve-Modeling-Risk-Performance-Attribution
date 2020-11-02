%load('EUR_IS_10YrCurves.mat')

%% klipp OOS vid 257 och infoga i IS
load('EUR_OOS_10YrCurves_Clean.mat')

%%
fAllTemp = fAll(1:256,:);
fDatesAllTemp = fDatesAll(1:256,:);
piAllTemp = piAll(1:256,:);
TAllTemp = TAll(1:256,:);
tradeDatesAllTemp = tradeDatesAll(1:256);
zAllTemp =  zAll(1:256,:);

clearvars -except fAllTemp fDatesAllTemp piAllTemp TAllTemp tradeDatesAllTemp zAllTemp
save('tempSaveOOS')
%%
load('EUR_OOS_10YrCurves_Clean.mat')
fAll = fAll(257:end-23,:);
fDatesAll = fDatesAll(257:end-23,:);
piAll = piAll(257:end-23,:);
TAll = TAll(257:end-23,:);
tradeDatesAll = tradeDatesAll(257:end-23);
zAll = zAll(257:end-23,:);

clearvars -except fAll fDatesAll piAll TAll tradeDatesAll zAll currency fileName measurementPath
save('EUR_OOS_10YrCurves_Clean_Final')

%%
load('EUR_IS_10YrCurves.mat')

load('tempSaveOOS.mat')
%%

fAll = [fAll; fAllTemp];
fDatesAll = [fDatesAll; fDatesAllTemp];
piAll = [piAll; piAllTemp];
TAll = [TAll; TAllTemp];
tradeDatesAll = [tradeDatesAll tradeDatesAllTemp];
zAll = [zAll; zAllTemp];

%%
clearvars -except fAll fDatesAll piAll TAll tradeDatesAll zAll currency fileName measurementPath
save('EUR_IS_10YrCurves_Clean_Final')
%%
load('EUR_OOS_10YrCurves_Clean_Final')
%%
load('SEK_OOS_10YrCurves')
fAll = fAll(251:end-23,:);
fDatesAll = fDatesAll(251:end-23,:);
piAll = piAll(251:end-23,:);
TAll = TAll(251:end-23,:);
tradeDatesAll = tradeDatesAll(251:end-23);
zAll = zAll(251:end-23,:);

clearvars -except fAll fDatesAll piAll TAll tradeDatesAll zAll currency fileName measurementPath
save('SEK_OOS_10YrCurves_Clean_Final')


%%
load('USD_OOS_10YrCurves_Clean')

fAll = fAll(252:end,:);
fDatesAll = fDatesAll(252:end,:);
piAll = piAll(252:end,:);
TAll = TAll(252:end,:);
tradeDatesAll = tradeDatesAll(252:end);
zAll = zAll(252:end,:);


clearvars -except fAll fDatesAll piAll TAll tradeDatesAll zAll currency fileName measurementPath
save('USD_OOS_10YrCurves_Clean_Final')

