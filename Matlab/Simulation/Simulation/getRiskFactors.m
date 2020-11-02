function [E, E_k, DZero, DPi] = getRiskFactors(kZero, kPi, fAll_IS, piAll_IS, type)
    
    % Ändra utifrån hur många år som ska vara med nedan
    DZero = fAll_IS(end-1511:end,:) - fAll_IS(end-1512:end-1,:);
    DPi = piAll_IS(end-1511:end,:) - piAll_IS(end-1512:end-1,:);

    if type == 1
        CZero = cov(DZero);
        CPi = cov(DPi);
        [V,D] = eigs(CZero, size(CZero, 1));
        [e,ind] = sort(diag(D),1, 'descend');
        E = {};
        E.Zero = V(:,ind);
        [V,D] = eigs(CPi, size(CZero, 1));
        [e,ind] = sort(diag(D),1, 'descend');
        E.Pi = V(:,ind);

        E_k.Zero = E.Zero(:,1:kZero);
        E_k.Pi = E.Pi(:,1:kPi);
    
    elseif type == 2 || type == 3
    
        load('EZero.mat')
        load('E_k_Zero_IRAM.mat')
        load('EPi.mat')
        load('E_k_Pi_IRAM.mat')

        E = {};
        E.Zero = table2array(CopyofBDCSVDeigvecEURIS10YrfAll);
        E.Pi = table2array(BDCSVDeigvecEURIS10YrpiAll);
        E_k = {};
        
        if type == 2
        
            % BDCSV is used for E_k
            E_k.Zero = E.Zero(:,1:6);
            E_k.Pi = E.Pi(:,1:8);

        elseif type == 3
            
            E_k.Zero = table2array(IRAMeigvecEURIS10YrfAll);
            E_k.Pi = table2array(IRAMeigvecEURIS10YrpiAll);       
            
            % Cleaning if IRAM is used
            E_k.Zero(:,3) = -E_k.Zero(:,3);
            E_k.Zero(:,5) = -E_k.Zero(:,5);
            E_k.Pi(:,1) = -E_k.Pi(:,1);
            E_k.Pi(:,2) = -E_k.Pi(:,2);
            E_k.Pi(:,4) = -E_k.Pi(:,4);
            E_k.Pi(:,6) = -E_k.Pi(:,6);
            E_k.Pi(:,7) = -E_k.Pi(:,7);

        end
    end 
end

