clear;clc;close all;
s = rng(94);
%% Управление
flag_DN = 0; % Enable/disable построение ДН
flag_spectr = 0;% Enable/disable спектр OFDM
flag_Steering = 0;       % Enable/disable steering Улучшает BER если prm.numR > prm.numSTS
flag_prekod = 0; % Enable/disable прекодирование  Не меняет
flag_cor = 0;
%% Каналы 
flag_my_chanel = 'STATIC'; %'STATIC'; 'DINAMIC'; 'ONE'
flag_chanel = 0; % Enable/disable MIMOCHANEL Не работает
%% Параметры
prm.numUsers = 1; % const
prm.M = 16;% Порядок модуляции
prm.bps = log2(prm.M); % Коль-во бит на символ в секунду
prm.numSTS = 2 ; %  4/8/16/32/64
prm.numTx = 2; % Кол-во излучающих антен
prm.numRx = 2; % Кол-во приемных антен
K_norm = prm.numSTS/prm.numRx; % Нормировка по энергии
prm.n_package = 10; % Кол-во посылок
SNR = 0:30;% dB
%% Параметры OFDM 
prm.Msc = 450; % Кол-во поднессущих
prm.N_FFT = 512; % Длина FFT для OFDM
prm.N_symb_pack = 2; % Кол-во символов OFDM в 1 пакете от 1 антены
prm.Nsymb = prm.N_symb_pack*prm.numSTS ;%  Кол-во символов OFDM от системы
prm.Ngi = 10; % длина защитных интервалов

prm.n = prm.bps*prm.Nsymb*prm.Msc;% Длина бинарного потока
prm.Nsymb_pilot = 1;
prm.n_pilot = prm.Msc*prm.numSTS;
%% Вспомогательные параметры
prm.numSTSVec = prm.numSTS;
prm.fc = 28e9;               % 28 GHz system несущая
prm.fs = 40e6; % Частота дискретизации
m = -prm.N_FFT/2:prm.N_FFT/2 -1;
Fanaliz = m*prm.fs/prm.N_FFT ; % Анализируемые частоты ДПФ
prm.cLight = physconst('LightSpeed');
prm.lambda = prm.cLight/prm.fc;
%% 
[isTxURA,prm.expFactorTx,isRxURA,prm.expFactorRx] = helperArrayInfo(prm,true);
%Возвращает 1 если можно использовать URA для Tx и Rx антенн 
%expFactor во сколько раз антенн больше чем потоков.
%% Координаты базовой/мобильной станций
prm.posTx = [0;0;0];         % BS/Transmit array position, [x;y;z], meters
maxRange = 1000;            % all MSs within 1000 meters of BS
prm.mobileRange = randi([1 maxRange],1,prm.numUsers); 
% Angles specified as [azimuth;elevation], az=[-180 180], el=[-90 90] elevation - угол места
prm.mobileAngle = [rand(1,prm.numUsers)*360-180; ... в градусах
                    rand(1,prm.numUsers)*180-90];
prm.steeringAngle = prm.mobileAngle;% Transmit steering angle (Близкие к mobileAngle)
[xRx,yRx,zRx] = sph2cart(deg2rad(prm.mobileAngle(1)),...
            deg2rad(prm.mobileAngle(2)),prm.mobileRange);
prm.posRx = [xRx;yRx;zRx];

[toRxRange,toRxAng] = rangeangle(prm.posTx,prm.posRx);
spLoss = fspl(toRxRange,prm.lambda);
%% Весовые коэф для прд/прм ФАР 
[wT, arrayTx] = Transmit_Beam_Steering(prm,isTxURA,flag_DN);
[wR, arrayRx] = Receive_Beam_Steering(prm,isRxURA,toRxAng,flag_DN);
%% Создание канальной матрицы
prm.nRays = 300;             % Number of rays 
[H,H1] = create_chanel_old(flag_my_chanel,prm);
% prm.posTxElem = getElementPosition(arrayTx)/prm.lambda;
% prm.posRxElem = getElementPosition(arrayRx)/prm.lambda;
% for ind = 1:prm.n_package
%       H{ind}= scatteringchanmtx(prm.posTxElem,prm.posRxElem,prm.nRays);
%       H1{ind} = H{ind}(1:prm.numSTS,1:prm.numSTS); 
% end
% for ind = 1:prm.n_package
N_zad = 1; % второй луч задержан на 2 отсчета 0.5e-7;
tau = 1/prm.fs;
pdB = 0;
chan = comm.MIMOChannel('MaximumDopplerShift',0, ...
    'SampleRate',                prm.fs,...
    'SpatialCorrelation',false, ...
    'NumTransmitAntennas',prm.numTx, ...
    'NumReceiveAntennas',prm.numRx,...
    'RandomStream','mt19937ar with seed', ...
    'Seed',73);
% chan = comm.MIMOChannel('MaximumDopplerShift',0, ...
%     'SampleRate',                prm.fs,...
%     'PathDelays',                tau,...
%     'AveragePathGains',          pdB,...
%     'SpatialCorrelation',false, ...
%     'NumTransmitAntennas',prm.numTx, ...
%     'NumReceiveAntennas',prm.numRx,...
%     'RandomStream','mt19937ar with seed', ...
%     'Seed',73);
chan1 = comm.MIMOChannel('MaximumDopplerShift',0, ...
    'SampleRate',                prm.fs,...
    'PathDelays',                tau,...
    'AveragePathGains',          pdB,...
    'SpatialCorrelation',false, ...
    'NumTransmitAntennas',prm.numSTS, ...
    'NumReceiveAntennas',prm.numSTS,...
    'RandomStream','mt19937ar with seed', ...
    'Seed',73);
% end
%%
for indSNR = 1:length(SNR)  
    All_In_data = [];      
    All_Out_data_ZF = [];
    All_Out_data_MMSE = [];
    All_Out_data_ML = [];
    All_Out_data_KB = [];
    All_Out_data_ZF1 = [];
    All_Out_data_MMSE1 = [];
    All_Out_data_ML1 = [];
    All_Out_data_KB1 = [];
    for ind = 1:prm.n_package
        %% Формируем данные
        In_data = randi([0 1],prm.n,1); % Передаваемые данные
        All_In_data = [All_In_data ; In_data];
        In_Matrix = reshape(In_data,length(In_data)/prm.bps,prm.bps); %Группируем биты 
        In_Symbols = bi2de(In_Matrix);  % Входные данные в символах

        Pilot = randi([0 1],prm.n_pilot,1);%  набор пилотов в битах
        %% Модулятор
        dataMod = pskmod(In_Symbols,prm.M);% Модулятор PSK
        Pilot_Mod = pskmod(Pilot,2);% Модулятор BPSK для пилотов       
        %% Зондирование канала
        % Модулятор OFDM 
        preambleSigSTS = OFDMmod_matrix_view(Pilot_Mod,prm.numSTS,prm.Msc,prm.Ngi,prm.N_FFT);
        preambleSigSTS = preambleSigSTS.';
        % Repeat over numTx
        preambleSig = zeros(size(preambleSigSTS,1),prm.numTx);
        for i = 1:prm.numSTS
            preambleSig(:,(i-1)*prm.expFactorTx+(1:prm.expFactorTx)) = ...
                repmat(preambleSigSTS(:,i),1,prm.expFactorTx);
        end
        % Прохождение канала
        dataChan_Pilot = preambleSig*H{ind};
        dataChan_Pilot1 = preambleSigSTS*H1{ind};
        % Собственный  шум         
        [dataNoise_Pilot ,~] = my_awgn(dataChan_Pilot,SNR(indSNR));
%         dataNoise_Pilot = awgn(dataChan_Pilot,SNR(indSNR),'measured');
       [dataNoise_Pilot1,~]= my_awgn(dataChan_Pilot1,SNR(indSNR));
%         dataNoise_Pilot1 = awgn(dataChan_Pilot1,SNR(indSNR),'measured');
        % Демодулятор OFDM
        for i = 1:prm.numRx
            tmp_dataMod_Pilot(:,i) = OFDMdemod(dataNoise_Pilot(:,i).',prm.Nsymb_pilot,prm.Msc,prm.Ngi,...
                            prm.N_FFT,1);
        end
        for i = 1:prm.numSTS
            tmp_dataMod_Pilot1(:,i) = OFDMdemod(dataNoise_Pilot1(:,i).',prm.Nsymb_pilot,prm.Msc,prm.Ngi,...
                            prm.N_FFT,1);
        end
        % Оценка канальной матрицы по пилотам 
        Pilot_Mod1 = [];
        for i = 1:prm.numSTS
            Pilot_Mod1 = [Pilot_Mod1 Pilot_Mod(1+(i-1)*prm.Msc:prm.Msc+(i-1)*prm.Msc)];
        end
        H_estim_Z{ind} = pinv(Pilot_Mod1)*tmp_dataMod_Pilot;              
        [F,W] = diagbfweights(H_estim_Z{ind});
        H_estim_Z1{ind} = pinv(Pilot_Mod1)*tmp_dataMod_Pilot1;                
        [F1,W1] = diagbfweights(H_estim_Z1{ind});      
        %% Цифровое Прекодирование
        % Apply precoding weights to the subcarriers, assuming perfect feedback
        dataMod1 = reshape(dataMod,prm.Nsymb*prm.Msc/prm.numSTS,prm.numSTS);       
        dataANDpilot_Mod = [Pilot_Mod1; dataMod1 ];
        if flag_prekod ==1            
            preData = dataANDpilot_Mod*F;
            preData1 = dataANDpilot_Mod;
        else
            preData = dataANDpilot_Mod;
            preData1 = dataANDpilot_Mod;
        end
        %% Модулятор OFDM       
        In_dataMod_STS = [];
        for i = 1:prm.numSTS
            In_dataMod_STS(:,i) = OFDMmod(preData(:,i),prm.Nsymb/prm.numSTS+prm.Nsymb_pilot,prm.Msc,prm.Ngi,prm.N_FFT);
        end
        In_dataMod_STS1 = [];
        for i = 1:prm.numSTS
            In_dataMod_STS1(:,i) = OFDMmod(preData1(:,i),prm.Nsymb/prm.numSTS+prm.Nsymb_pilot,prm.Msc,prm.Ngi,prm.N_FFT);
        end
        
        if flag_spectr ==1
            Sp = fft(In_dataMod_STS(prm.Ngi+1:prm.N_FFT,1),prm.N_FFT);
            plot(Fanaliz,20*log10(fftshift(abs(Sp))));
        end
        % Repeat over numTx
        In_dataMod = zeros(size(In_dataMod_STS,1),prm.numTx);
        for i = 1:prm.expFactorTx
            In_dataMod(:,(i-1)*prm.numSTS+(1:prm.numSTS)) = In_dataMod_STS;
        end
        %% Формируем луч на передачу
        if flag_Steering==1
            In_dataMod = In_dataMod.*conj(wT).';                   
        end
        In_dataMod1 = In_dataMod_STS1;
        %% Прохождение канала
        if flag_chanel == 1
            dataChan = chan(In_dataMod);
            dataChan1 = chan1(In_dataMod1);
        else
            dataChan = In_dataMod*H{ind};
            dataChan1 = In_dataMod1*H1{ind};
        end       
        %% Собственный  шум
        [dataNoise,sigma] = my_awgn(dataChan,SNR(indSNR));
%         dataNoise = awgn(dataChan,SNR(indSNR),'measured');
        [dataNoise1,sigma1] = my_awgn(dataChan1,SNR(indSNR));
%         dataNoise1 = awgn(dataChan1,SNR(indSNR),'measured');
        %% Прием
        if flag_Steering==1
            dataNoise = dataNoise.*conj(wR).';                  
        end       
        %% Демодулятор OFDM
        tmp_dataMod = [];
        for i = 1:prm.numRx
            tmp_dataMod(:,i) = OFDMdemod(dataNoise(:,i).',prm.N_symb_pack+prm.Nsymb_pilot,prm.Msc,prm.Ngi,...
                            prm.N_FFT,1);          
        end
        tmp_dataMod1 = [];
        for i = 1:prm.numSTS
            tmp_dataMod1(:,i) = OFDMdemod(dataNoise1(:,i).',prm.N_symb_pack+prm.Nsymb_pilot,prm.Msc,prm.Ngi,...
                            prm.N_FFT,1);          
        end  

        %% Оценка канальной матрицы по пилотам
        H_estim{ind} = pinv(Pilot_Mod1)*tmp_dataMod(1:size(Pilot_Mod1,1),:);                                       
        H_estim1{ind} = pinv(Pilot_Mod1)*tmp_dataMod1(1:size(Pilot_Mod1,1),:);
        if flag_cor ==0
            tmp_dataMod(1:size(Pilot_Mod1,1),:)= []; 
            tmp_dataMod1(1:size(Pilot_Mod1,1),:)= [];
        end
        %% Демодулятор MIMO(Эквалайзер MIMO) 
        %ZF
        tic        
        [Out_ZF] = My_MIMO_Equalize_ZF(tmp_dataMod,H_estim{ind},1);       
        time_ZF = toc;
%         %MMSE
%         tic
%         [Out_MMSE] = My_MIMO_Equalize_MMSE(tmp_dataMod,H_estim{ind},sigma,1);
%         time_MMSE = toc;  
%         %ML
%         tic
%         [Out_ML] = My_MIMO_Equalize_ML(tmp_dataMod,H_estim{ind},prm.M,prm.numSTS);
%         time_ML = toc;
%         %KB
%         num_k = 2;
%         tic
%         [Out_KB] = My_MIMO_Equalize_K_BEST(tmp_dataMod,num_k,H_estim{ind},prm.numSTS,prm.M);
%         time_KB = toc;
        %ZF
        tic        
        [Out_ZF1 ] = My_MIMO_Equalize_ZF(tmp_dataMod1,H_estim1{ind},1);
        time_ZF1 = toc;
%         %MMSE
%         tic
%         [Out_MMSE1] = My_MIMO_Equalize_MMSE(tmp_dataMod1,H_estim1{ind},sigma1,1);
%         time_MMSE1 = toc;
%         %ML
%         tic
%         [Out_ML1] = My_MIMO_Equalize_ML(tmp_dataMod1,H_estim1{ind},prm.M,prm.numSTS);
%         time_ML1 = toc;
%         %KB
%         tic
%         [Out_KB1] = My_MIMO_Equalize_K_BEST(tmp_dataMod1,num_k,H_estim1{ind},prm.numSTS,prm.M);
%         time_KB1 = toc;
        %%
        if flag_cor ==1
        % Оценка частотной характеристики для каждого потока
        for j = 1:prm.numSTS          
            Freq{ind}(:,j) = Out_ZF(1:size(Pilot_Mod1,1),j)./Pilot_Mod1(:,j);            
        end
        Out_ZF(1:size(Pilot_Mod1,1),:)= [];
        for j = 1:prm.numSTS 
            Freq1{ind}(:,j) = Out_ZF1(1:size(Pilot_Mod1,1),j)./Pilot_Mod1(:,j); 
        end
        Out_ZF1(1:size(Pilot_Mod1,1),:)= [];
        % Эквалайзер SISO
            for i = 1:prm.numSTS
                Out_ZF(:,i) = Out_ZF(:,i).*Freq{ind}(:,i);
                Out_ZF1(:,i) = Out_ZF1(:,i).*Freq1{ind}(:,i);
            end
        end
        %% Демодулятор
        Out_dataMod_ZF = Out_ZF(:);
%         scatterplot(Out_dataMod_ZF)
%         Out_dataMod_MMSE = Out_MMSE(:);
        Out_Symbols_ZF = pskdemod(Out_dataMod_ZF,prm.M );% Демодулятор 
%         Out_Symbols_MMSE = pskdemod(Out_dataMod_MMSE,prm.M);
        
        Out_dataMod_ZF1 = Out_ZF1(:);
%         Out_dataMod_MMSE1 = Out_MMSE1(:);
        Out_Symbols_ZF1 = pskdemod(Out_dataMod_ZF1,prm.M );% Демодулятор 
%         Out_Symbols_MMSE1 = pskdemod(Out_dataMod_MMSE1,prm.M);
%         Out_Symbols_ML = Out_ML(:);
%         Out_Symbols_ML1 = Out_ML1(:);
%         Out_Symbols_KB = Out_KB(:);
%         Out_Symbols_KB1 = Out_KB1(:);
        %% Выходные данные
        Out_Matrix_ZF = de2bi(Out_Symbols_ZF); %Разгруппировка из символов в поток бит
%         Out_Matrix_MMSE = de2bi(Out_Symbols_MMSE);
%         Out_Matrix_ML = de2bi(Out_Symbols_ML);
%         Out_Matrix_KB = de2bi(Out_Symbols_KB);
        Out_Matrix_ZF1 = de2bi(Out_Symbols_ZF1); %Разгруппировка из символов в поток бит
%         Out_Matrix_MMSE1 = de2bi(Out_Symbols_MMSE1);
%         Out_Matrix_ML1 = de2bi(Out_Symbols_ML1);
%         Out_Matrix_KB1 = de2bi(Out_Symbols_KB1);

        Out_data_ZF = Out_Matrix_ZF(:);   % Выходные данные в битах MIMO        
%         Out_data_MMSE = Out_Matrix_MMSE(:);
%         Out_data_ML = Out_Matrix_ML(:);
%         Out_data_KB = Out_Matrix_KB(:);
        Out_data_ZF1 = Out_Matrix_ZF1(:);   % Выходные данные в битах MIMO        
%         Out_data_MMSE1 = Out_Matrix_MMSE1(:);
%         Out_data_ML1 = Out_Matrix_ML1(:);
%         Out_data_KB1 = Out_Matrix_KB1(:);
        
        All_Out_data_ZF = [All_Out_data_ZF ; Out_data_ZF];
%         All_Out_data_MMSE = [All_Out_data_MMSE ; Out_data_MMSE];
%         All_Out_data_ML = [All_Out_data_ML ; Out_data_ML];
%         All_Out_data_KB = [All_Out_data_KB ; Out_data_KB];
        All_Out_data_ZF1 = [All_Out_data_ZF1 ; Out_data_ZF1];
%         All_Out_data_MMSE1 = [All_Out_data_MMSE1 ; Out_data_MMSE1];
%         All_Out_data_ML1 = [All_Out_data_ML1 ; Out_data_ML1];
%         All_Out_data_KB1 = [All_Out_data_KB1 ; Out_data_KB1];
    end
    [~,ber_ZF(indSNR)] = biterr(All_In_data,All_Out_data_ZF);
%     [~,ber_MMSE(indSNR)] = biterr(All_In_data,All_Out_data_MMSE);
%     [~,ber_ML(indSNR)] = biterr(All_In_data,All_Out_data_ML);
%     [~,ber_KB(indSNR)] = biterr(All_In_data,All_Out_data_KB);
    [~,ber_ZF1(indSNR)] = biterr(All_In_data,All_Out_data_ZF1);
%     [~,ber_MMSE1(indSNR)] = biterr(All_In_data,All_Out_data_MMSE1);
%     [~,ber_ML1(indSNR)] = biterr(All_In_data,All_Out_data_ML1);
%     [~,ber_KB1(indSNR)] = biterr(All_In_data,All_Out_data_KB1);
    indSNR
end
% Теоритическая BER 
EbN0_dB = 0:11;
EbN0 = 10.^(EbN0_dB/10);
tber = 0.5.*erfc(sqrt(EbN0));
%
load('PSK_16_awgn');
figure()
semilogy(SNR-(10*log10(prm.bps*K_norm)),ber_ZF,'r','LineWidth',1.5);%'-*k'
hold on
semilogy(SNR-(10*log10(prm.bps)),ber_ZF1,'b','LineWidth',1.5);%'-.k'
semilogy(ber_theor.data{1} ,ber_theor.data{2},'k','LineWidth',1.5)%-(10*log10(prm.numSTS))
% semilogy(SNR-(10*log10(prm.bps)),ber_MMSE1,'*b');%+10-ceil(10*log10(prm.bps*prm.numSTS))
% semilogy(SNR-(10*log10(prm.bps*K_norm)),ber_MMSE,'*r');
% semilogy(SNR,ber_ML,'g','LineWidth',1.5);
% semilogy(SNR,ber_ML1,'c','LineWidth',1.5);
% semilogy(SNR,ber_KB,'*g');
% semilogy(SNR,ber_KB1,'*c');
grid on
xlim([0 SNR(28)])
ylim([10^-5 10^0]);
% legend('MASSIV ZF','MASSIV MMSE','MIMO ZF','MIMO MMSE','theor')
str1 = ['MASSIV MIMO ' num2str(prm.numTx) 'x'  num2str(prm.numRx)];
str2 = ['MIMO ' num2str(prm.numSTS) 'x'  num2str(prm.numSTS)];
str11 = ['MM ' num2str(prm.numTx) 'x'  num2str(prm.numRx)];
str22 = ['M ' num2str(prm.numSTS) 'x'  num2str(prm.numSTS)];
% if flag_static==1
%    str3 = 'Статический канал';
% else
%    str3 = 'Динамический канал';
% end
% legend(['ZF ' str11],['ZF ' str22],'Теоретическая',['MMSE ' str11],['MMSE ' str22],...
%         ['ML ' str11],['ML ' str22],['KBEST ' str11],['KBEST ' str22]);
xlabel('E_b / N_0 , dB')
ylabel('BER')
% title(str3)
