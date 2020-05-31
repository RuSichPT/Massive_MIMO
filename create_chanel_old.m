function [H, H1] = create_chanel(flag_my_chanel,prm)
switch flag_my_chanel
    case 'STATIC'
        rng(94); %167-сид(начальное число)
        for i= 1:prm.numTx
            for j = 1:prm.numRx           
                H_tmp(i,j) = (randn(1)+1i*randn(1))/sqrt(2);
            end
        end
        for ind = 1:prm.n_package
            H{ind} =  H_tmp(:,:);
        end    
       rng(94); %167-сид(начальное число)
        for i = 1:prm.numSTS
            for j = 1:prm.numSTS          
                H_tmp1(i,j) = (randn(1)+1i*randn(1))/sqrt(2);
            end
        end
        for ind = 1:prm.n_package
            H1{ind} =  H_tmp1(:,:);
        end
    case 'DINAMIC'
        rng(94);
        for ind = 1:prm.n_package        
            for i= 1:prm.numTx
                for j = 1:prm.numRx           
                     H{ind}(i,j) = (randn(1)+1i*randn(1))/sqrt(2);                 
                end
            end
        end
        rng(94);
        for ind = 1:prm.n_package        
            for i = 1:prm.numSTS
                for j = 1:prm.numSTS          
                    H1{ind}(i,j) = (randn(1)+1i*randn(1))/sqrt(2);
                end
            end
        end
    case 'ONE'        
        for ind = 1:prm.n_package
              H{ind} = eye(prm.numTx,prm.numRx);
              H1{ind} = H{ind}(1:prm.numSTS,1:prm.numSTS);
        end
end
end    



