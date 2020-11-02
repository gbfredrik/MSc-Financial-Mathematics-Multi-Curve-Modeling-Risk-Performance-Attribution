function [floatCashFlows, fixCashFlows, yield, fixingDatesCashFlows, RoP, IborDates, Ibor, Nom] = getPortfolioData(instruments, ccy)

    floatCashFlows = {};
    fixCashFlows = {};
    yield = [];
    fixingDatesCashFlows = {};
    RoP = {};

    loop = instruments;
    for i = 1:length(loop)
        currColumn = strcat(loop(i), ':', loop(i));
        floatCashFlows{i} = rmmissing(table2array(readtable('PortfolioData.xlsx', 'Sheet', strcat("Float Cashflows ", ccy), 'Range', currColumn)));
        floatCashFlows{i} = floatCashFlows{i}(2:end);
        fixCashFlows{i} = rmmissing(table2array(readtable('PortfolioData.xlsx', 'Sheet', strcat("Fix Cashflows ", ccy), 'Range', currColumn)));
        fixCashFlows{i} = fixCashFlows{i}(2:end);
        yield(i) = rmmissing(table2array(readtable('PortfolioData.xlsx', 'Sheet', strcat("Yield ", ccy), 'Range', currColumn))) / 100;
        fixingDatesCashFlows{i} = rmmissing(table2array(readtable('PortfolioData.xlsx', 'Sheet', strcat("Fixing dates ", ccy), 'Range', currColumn)));
        RoP(i) = rmmissing(table2array(readtable('PortfolioData.xlsx', 'Sheet', strcat("RoP ", ccy), 'Range', currColumn)));
    end

    IborDates = rmmissing(table2array(readtable('PortfolioData.xlsx', 'Sheet', strcat("IBOR ", ccy), 'Range', 'A:A')));
    Ibor = rmmissing(table2array(readtable('PortfolioData.xlsx', 'Sheet', strcat("IBOR ", ccy), 'Range', 'B:B'))) / 100;
    Nom = 1000 * ones(length(loop), 1);


end

