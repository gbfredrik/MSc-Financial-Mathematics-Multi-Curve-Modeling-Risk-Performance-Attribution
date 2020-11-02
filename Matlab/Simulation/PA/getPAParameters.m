function [paParams] = getPAParameters(paParams, A, f, pi)
    
    % Set prev yield curve to current
    paParams{16} = paParams{14};
    paParams{17} = paParams{15};  
    paParams{20} = A * paParams{16};
    paParams{21} = A * paParams{17}; 
    
    % Set new yield curves
    paParams{14} = f;
    paParams{15} = pi;
    paParams{18} = A * f;
    paParams{19} = A * pi;

end

