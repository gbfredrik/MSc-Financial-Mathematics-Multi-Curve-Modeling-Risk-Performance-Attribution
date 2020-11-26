function [floatCashFlows, fixCashFlows, yield, fixingDatesCashFlows, RoP, IborDates, Ibor, Nom] = getPortfolioData(instruments, ccy, paType)
    floatCashFlows = {};
    fixCashFlows = {};
    yield = [];
    fixingDatesCashFlows = {};
    RoP = {};

    loop = instruments;
    for i = 1:length(loop)
        currColumn = strcat(loop(i), ':', loop(i));

        if verLessThan('matlab', '9.8')
            floatCashFlows{i} = rmmissing(table2array(readtable('PortfolioData.xlsx', 'Sheet', strcat("Float Cashflows ", ccy), 'Range', currColumn)));
            floatCashFlows{i} = floatCashFlows{i}(2:end);
            fixCashFlows{i} = rmmissing(table2array(readtable('PortfolioData.xlsx', 'Sheet', strcat("Fix Cashflows ", ccy), 'Range', currColumn)));
            fixCashFlows{i} = fixCashFlows{i}(2:end);
            yield(i) = rmmissing(table2array(readtable('PortfolioData.xlsx', 'Sheet', strcat("Yield ", ccy), 'Range', currColumn))) / 100;
            fixingDatesCashFlows{i} = rmmissing(table2array(readtable('PortfolioData.xlsx', 'Sheet', strcat("Fixing dates ", ccy), 'Range', currColumn)));
            RoP(i) = rmmissing(table2array(readtable('PortfolioData.xlsx', 'Sheet', strcat("RoP ", ccy, " ", paType), 'Range', currColumn)));
        else
            floatCashFlows{i} = rmmissing(table2array(readtable('PortfolioData.xlsx', 'Sheet', strcat("Float Cashflows ", ccy), 'Range', currColumn)));
            floatCashFlows{i} = floatCashFlows{i}(3:end);
            fixCashFlows{i} = rmmissing(table2array(readtable('PortfolioData.xlsx', 'Sheet', strcat("Fix Cashflows ", ccy), 'Range', currColumn)));
            fixCashFlows{i} = fixCashFlows{i}(3:end);
            yield(i) = rmmissing(table2array(readtable('PortfolioData.xlsx', 'Sheet', strcat("Yield ", ccy), 'Range', currColumn, 'NumHeaderLines', 1))) / 100;
            fixingDatesCashFlows{i} = rmmissing(table2array(readtable('PortfolioData.xlsx', 'Sheet', strcat("Fixing dates ", ccy), 'Range', currColumn)));
            fixingDatesCashFlows{i} = fixingDatesCashFlows{i}(2:end);
            RoP(i) = rmmissing(table2array(readtable('PortfolioData.xlsx', 'Sheet', strcat("RoP ", ccy, " ", paType), 'Range', currColumn, 'format', 'auto')));
        end
    end
    
    IborDates = rmmissing(table2array(readtable('PortfolioData.xlsx', 'Sheet', strcat("IBOR ", ccy), 'Range', 'A:A')));
    Ibor = rmmissing(table2array(readtable('PortfolioData.xlsx', 'Sheet', strcat("IBOR ", ccy), 'Range', 'B:B'))) / 100;
    Nom = 1000 * ones(length(loop), 1);
end
