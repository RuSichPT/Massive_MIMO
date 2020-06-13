clear;clc;close all;
% rng(67)
prm.numTx = 4; % Кол-во излучающих антен 
prm.numRx = 4; % Кол-во приемных антен
prm.numSTS = 2; % Кол-во потоков 2/4/8/16/32/64
prm.M = 16;% Порядок модуляции
prm.bps = log2(prm.M); % Коль-во бит на символ в секунду
symb = 10000;
prm.n = prm.bps*symb;% Длина бинарного потока
%%
c = 3e8;        % propagation speed
fc = 60e9;      % carrier frequency
lambda = c/fc;  % wavelength
prm.nRays = 10;
txcenter = [0;0;0];
rxcenter = [1500;500;0];
txarray = phased.ULA('NumElements',4,'ElementSpacing',lambda/2);
txmipos = getElementPosition(txarray)/lambda;
rxarray = phased.ULA('NumElements',4,'ElementSpacing',lambda/2);
rxmopos = getElementPosition(rxarray)/lambda;
chan = scatteringchanmtx(txmipos,rxmopos,prm.nRays);
%% Формируем данные
expFactor = prm.numTx/prm.numSTS ;
Inp_data = randi([0 1],prm.n,prm.numSTS); % Передаваемые данные
Inp_data = repmat(Inp_data,[1 expFactor]);
all_Inp_data = reshape(Inp_data,size(Inp_data,1)*prm.numRx,1);
for m = 1:prm.numTx 
    Inp_Mat = reshape(Inp_data(:,m),length(Inp_data(:,m))/prm.bps,prm.bps); %Группируем биты
    Inp_Sym(:,m) = bi2de(Inp_Mat);  % Входные данные в символах
end
%% Модулятор
Mod_data_inp = qammod(Inp_Sym,prm.M);% Модулятор QAM-M для полезной инф
Mod_data_inp_tx = Mod_data_inp;
% Mod_data_inp_tx = reshape(Mod_data_inp,length(Mod_data_inp)/(4),4);
%% Канал
wp = ones(1,prm.numTx);
wc = ones(prm.numRx,1);
[wt,wr] = diagbfweights(chan);
A = sqrt(sum(abs(Mod_data_inp).^2)/(size(Mod_data_inp,1)));
Mod_data_inp_pre_tx = Mod_data_inp_tx*wt;
H_STS = repmat(eye(prm.numSTS,prm.numSTS),[expFactor 1]);
snr_param = 0:40;
h_dig = wt*chan*wr
all_Out_data = [];
for m=1:numel(snr_param)
%     n = A*sqrt(db2pow(-snr_param(m))/2)*(randn(symb,prm.numRx)+1i*randn(symb,prm.numRx));
    Chanel_data = Mod_data_inp_pre_tx*chan*wr;%+n*wr;
    Chanel_data = awgn(Chanel_data,snr_param(m),'measured');
    Eq_data = Chanel_data*pinv(h_dig);
    Chanel_data_m = Mod_data_inp_tx*chan;
    Chanel_data_m = awgn(Chanel_data_m,snr_param(m),'measured');
    Eq_data_m = Chanel_data_m*pinv(chan);
    Chanel_data_sts = 1/sqrt(prm.numSTS)*Mod_data_inp_tx*chan;
    Chanel_data_sts = awgn(Chanel_data_sts,snr_param(m),'measured');
    Eq_data_sts = Chanel_data_sts*pinv(chan);
    Eq_data_sts1 = Eq_data_sts*H_STS;
    all_Out_data = [];all_Out_data_m = [];all_Out_data_sts=[];
    for i = 1:prm.numRx
%         Eq_data = Chanel_data(:,i)./h_dig(i,i);% нужно лишь для all_Out_data
        Out_Sym(:,i) = qamdemod(Eq_data(:,i),prm.M);
        Out_Mat = de2bi(Out_Sym(:,i));
        Out_data = Out_Mat(:);
        all_Out_data = [all_Out_data;Out_data];
        [~,ber(i,m)] = biterr(Inp_data(:,i),Out_data);
        Out_Sym_m  = qamdemod(Eq_data_m(:,i),prm.M);
        Out_Mat_m  = de2bi(Out_Sym_m);
        Out_data_m  = Out_Mat_m(:);
        all_Out_data_m = [all_Out_data_m;Out_data_m];
        [~,ber_m(i,m)] = biterr(Inp_data(:,i),Out_data_m);
    end
    for i = 1:prm.numSTS
        Out_Sym_sts  = qamdemod(Eq_data_sts1(:,i),prm.M);
        Out_Mat_sts  = de2bi(Out_Sym_sts );
        Out_data_sts  = Out_Mat_sts (:);
        all_Out_data_sts  = [all_Out_data_sts ;Out_data_sts ];
        [~,ber_sts(i,m)] = biterr(Inp_data(:,i),Out_data_sts);
    end
    
    [~,ber_mean(m)] = biterr(all_Inp_data,all_Out_data);
    [~,ber_mean_m(m)] = biterr(all_Inp_data,all_Out_data_m);
    [~,ber_mean_sts(m)] = biterr(all_Inp_data(1:prm.n*prm.numSTS),all_Out_data_sts); 
end

for i = 1:prm.numRx
semilogy(snr_param,ber(i,:));
hold on
end
semilogy(snr_param,ber_mean);
grid on
legend('1','2','3','4','mean')
figure()
for i = 1:prm.numRx
semilogy(snr_param,ber_m(i,:));
hold on
end
semilogy(snr_param,ber_mean_m);
grid on
legend('1','2','3','4','mean m')
figure()
for i = 1:prm.numSTS
semilogy(snr_param,ber_sts(i,:));
hold on
end
semilogy(snr_param,ber_mean_sts);
grid on
legend('1','2','mean sts')
