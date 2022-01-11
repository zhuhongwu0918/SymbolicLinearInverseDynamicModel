function PrintParameterSummary(model, N, M, V, C, a)
%% PrintParameterSummary(model, N, M, V, C, a)
% Prints a parameter summary from the outputs of the RPNA algorithm

    for i = 1:model.NB
        
        fprintf('===============================\n');
        fprintf('Minimal Parameters for Body %d\n',i);
        x=0;
        for k = 1:10
            if norm(M{i}(k,:)) > 0%M的第{1}个，如果这一行不是全为0，则属于最小参数集内，Minimal Parameters for Body 1
                fprintf('%s\n',a{k});
                x=x+1;
            end
        end
        if x == 0
            fprintf('None\n');
        end

        fprintf('\nIdentifiable Parmeters for Body %d\n',i);
        x=0;
        for k = 1:10
           if IsIdentifiable(model,N,i,k)%判断是否属于可独立辨识参数
               %如果这N这一行有值,则不可辨识
              fprintf('%s\n',a{k});
              x=x+1;
           end
        end
        if x == 0
            fprintf('None\n');
        end


        fprintf('\nUnidentifiable Parmeters for Body %d\n',i);
        x=0;
        for k = 1:10
           if IsUnidentifiable(model,N,i,k)%判断是否属于可辨识参数
              fprintf('%s\n',a{k});
              x=x+1;
           end
        end
        if x == 0
            fprintf('None\n');
        end

        fprintf('\nDim Null N(%d) = %d\n',i,10-size(N{i},1));
        fprintf('Dim VelocitySpan V(%d) = %d\n',i, rank(V{i}));  
        fprintf('Dim OuterProductSpan K(%d) = %d\n\n',i, rank(C{i}));  
        
    end