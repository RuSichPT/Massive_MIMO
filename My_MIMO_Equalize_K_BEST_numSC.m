function [Out_Symbols_KB] = My_MIMO_Equalize_K_BEST_numSC(Y,num_k,H_estim,M,N,Inp_data)
% Модель Y = X*H+ksi; модуляция - qam
% y = x*H   y.' = H.'*x.' поэтому H_estim.' dataNoise.'

% Эквалайзер для каждой поднесущей
% Y - принятые символы [msc,symb_ofdm,numTx]
% H_estim - оценка матрицы вида [msc,numTx,numRx]
% msc - кол-во поднесущих,symb_ofdm - кол-во символов ofdm
% num_k - номер K-best
% N - numTx;
% M - модуляция
y = permute(Y,[3 2 1]);
H_estim = permute(H_estim,[2 3 1]);
% h_estim = squeeze(H_estim(1,:,:));
Out_Symbols_KB = [];
for ind = 1:size(H_estim,3)
    h_estim = squeeze(H_estim(:,:,ind));
    [Q,R(:,:,ind)] = qr(h_estim); %  Разложили на QR  y = QR*x
    y_tilda(:,:,ind) = Q'*squeeze(y(:,:,ind)); % Q'*Q = 1 получаем y = R*x
end
for indSymb = 1:size(y,2)
    for i = 1:size(y,3)
        r = squeeze(y_tilda(:,indSymb,i));
        R_est = squeeze(R(:,:,i));
        allSymb_help{1} = 0:M -1;  % Возможные варианты x для последней антены 
        allTxSig_help{1} = pskmod(allSymb_help{1},M);
        for j =1:N % цикл по антенам     
            for g = 1:num_k^(j-1)
                Symb_help = allSymb_help{g}; 
                TxSig_help = allTxSig_help{g};
                help_r = repmat(r(N-(j-1):N),[1,size(TxSig_help,2)]);
                for ind_k = 1:num_k % цикл по к            
                    if j==1 % первая итерация
                        [~, k] = min(abs(help_r - R_est(N,N)*TxSig_help).^2);                               
                    else % последующие
                        [~, k] = min(sum(abs(help_r - R_est(N-(j-1):N,N-(j-1):N)*TxSig_help).^2,1)); 
                    end              
                    est1{g,j}(:,ind_k) = Symb_help(:,k); %Оценка сигналов                    
                    TxSig_help(:,k) = 100;                         
                end
            end
            if j ~= N
                allSymb_help = [];
                ii = 1;
                for g = 1:num_k^(j-1)
                    for ind_k = 1:num_k
                        allSymb_help{ii} = [0:M-1;repmat(est1{g,j}(:,ind_k),[1,M])]; % Возможные варианты x 
                        allTxSig_help{ii} = qammod(allSymb_help{ind_k},M);
                        ii = ii+1;
                    end
                end        
            else
                allSymb_end = [];
                allTxSig_end = [];
                for g=1:num_k^(j-1)
                    allSymb_end = [allSymb_end est1{g,j}];                 
                end
                allTxSig_end = qammod(allSymb_end,M);
            end
        end
        help_r1 = repmat(r,[1,size(allTxSig_end,2)]);    
        [~, k] = min(sum(abs(help_r1 - R_est*allTxSig_end).^2,1));
        estKB = allSymb_end(:,k);
        Out_Symbols_KB = [Out_Symbols_KB; de2bi(estKB,log2(M),'left-msb').'];  
    end
end
% Выходные данные в символах от 0 до M-1
end
