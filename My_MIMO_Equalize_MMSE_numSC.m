function [Out_dataMod_MMSE] = My_MIMO_Equalize_MMSE_numSC(Y,H_estim,sigma)
% Модель Y = X*H+ksi; 

% Эквалайзер для каждой поднесущей
% Y - принятые символы [msc,symb_ofdm,numTx]
% H_estim - оценка матрицы вида [msc,numTx,numRx]
% msc - кол-во поднесущих,symb_ofdm - кол-во символов ofdm
% sigma - СКО шума (диагональная)
for i = 1:size(Y,1)    
    h_estim = squeeze(H_estim(i,:,:));
    inv_H_MMSE = inv(h_estim*h_estim'+2*sigma(1:size(h_estim,1),1:size(h_estim,1))^2); %% НЕ ВСЕ СИГМЫ УЧАВСТУЮТ
    Out_dataMod_MMSE(i,:,:) =  squeeze(Y(i,:,:))*h_estim'*inv_H_MMSE;
end
% Выходные данные в символах IQ
end

