function [Out_dataMod_ZF] = My_MIMO_Equalize_ZF(Y,H_estim,fast)
% Модель Y = X*H+ksi; 

% Y - принятые символы [поток,numTx]
% H_estim - оценка матрицы вида [numTx,numRx,поток]
if fast ==1
    inv_H_ZF = inv(H_estim*H_estim');
    Out_dataMod_ZF =  Y*H_estim'*inv_H_ZF;
else
    for i = 1:size(Y,1)    
        h_estim = squeeze(H_estim(:,:,i));
    %     tic
        inv_H_ZF = inv(h_estim*h_estim');
        Out_dataMod_ZF(i,:) =  Y(i,:)*h_estim'*inv_H_ZF;
    %     time_ZF = toc;
    end
end
% Выходные данные в символах IQ
end

