function [simParams] = initSimParams(N, E, DZero, DPi, fAll_IS, piAll_IS, fAll_OOS, piAll_OOS, tradeDatesAll_OOS)

    T = 1:3650;
    simParams = {}; 
    for i = 1:8
        for j = 1:2
            simParams{i}{j} = 0;
        end
    end
    simParams{9}{1} = DZero;
    simParams{9}{2} = DPi;
    simParams{10}{1} = E.Zero;
    simParams{10}{2} = E.Pi;
    simParams{11} = N;
    simParams{12} = T;
    simParams{13}{1} = 0;
    simParams{13}{2} = 0;
    
    % Initialize sim params
    [~, ~, simParams] = TermStructureSim(1, simParams, fAll_IS, piAll_IS, fAll_OOS, piAll_OOS, tradeDatesAll_OOS, true);
    
end

