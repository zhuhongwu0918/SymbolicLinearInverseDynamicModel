function N = SwitchedObservabilityMatrix(C,A)
%% N = SwitchedObservabilityMatrix(C,A)
% Compute the switched Observability matrix associated with the pairs
% (C, A_i) where A_i are the elements of the cell array A.

    AT = A;
    for i = 1:length(A)
        AT{i} = A{i}';
    end
    
    %Duality!AT就是Aphi, C'是Ki_
    N = SwitchedControllabilityMatrix(AT,C')';
end
    