function R = SwitchedControllabilityMatrix(A,B)
%% R = SwitchedControllabilityMatrix(A,B)
% Compute the switched Controllability matrix associated with the pairs
% (A_i, B) where A_i are the elements of the cell array A.

    R = RangeBasis(B); % Range of R gives Reachable Subspace
    r = 0;
    % While subspace dimension is growing
    while rank(R) > r
        r = rank(R);
        R_new = R;
        for i = 1:length(A)%1：6
            R_new = [R_new, A{i}*R];%R是Vi_,A{i}是关节自由度的叉乘Phi×
        end
        R = RangeBasis(R_new);
    end
end
    