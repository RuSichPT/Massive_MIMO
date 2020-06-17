function [Out_Symbols_ML] = My_MIMO_Equalize_ML_numSC(Y,H_estim,M,numTx)
% Модель Y = X*H+ksi; модуляция - qam

% Эквалайзер для каждой поднесущей
% Y - принятые символы [msc,symb_ofdm,numTx]
% H_estim - оценка матрицы вида [msc,numTx,numRx]
% msc - кол-во поднесущих,symb_ofdm - кол-во символов ofdm
% M - порядок модуляции
% numTx - кол-во передающих антенн

Out_Symbols_ML = [];
bps = log2(M);
S = 0:2^(bps*numTx)-1; %Все возможные варианты x в десятичном виде 
allBits = de2bi(S, 'left-msb')'; %Все возможные варианты x в битах 
for i =1:numTx
    tmp = allBits(1+(i-1)*bps:bps+(i-1)*bps,:); % tmp - временная переменная
    allSymb(i,:) = bi2de(tmp','left-msb') ; %Все возможные варианты x в символах от 0 до M-1
end
allTxSig = qammod(allSymb,M); %Все возможные варианты x в символах IQ
allTxSig = allTxSig.';

for j = 1:size(Y,2)
    for i = 1:size(Y,1)    
        r = squeeze(Y(i,j,:)).';
        h_estim = squeeze(H_estim(i,:,:));
    %     tic
        [~, k] = min(sum(abs(repmat(r,[2^(bps*numTx),1]) - allTxSig*h_estim).^2,2));
        % Суммирование происходит по антеннам 
        estML = allSymb(:,k).';
    %     time_ML = toc;
        Out_Symbols_ML = [Out_Symbols_ML; de2bi(estML,bps,'left-msb').'];    
    end
end
% Выходные данные в символах от 0 до M-1
end

