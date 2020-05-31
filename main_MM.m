%% ---------Модель Massiv MIMO and MIMO-------- 
clear;clc;%close all;
%% Управление
flag_DN = 0; % Enable/disable построение ДН
flag_Steering = 1;       % Enable/disable steering
flag_chanel = 'RAYL_SPECIAL'; % 'AWGN' ,'RAYL','RIC','RAYL_SPECIAL','STATIC', 'BAD' 
flag_cor_MIMO = 1; % 1-коррекция АЧХ (эквалайзер для MIMO)
flag_cor_MM = 1; % 1-коррекция АЧХ (эквалайзер для M MIMO)
flag_wav_MIMO = 1; % вейвлет шумоподавление для MIMO
flag_wav_MM = 1; % вейвлет шумоподавление для M MIMO
%% Параметры системы 
prm.numTx = 4; % Кол-во излучающих антен 
prm.numRx = 4; % Кол-во приемных антен
prm.numSTS = 2; % Кол-во потоков  4/8/16/32/64
prm.M = 16;% Порядок модуляции
prm.bps = log2(prm.M); % Коль-во бит на символ в секунду
prm.LEVEL = 3;% Уровень декомпозиции вейвлет шумоподавления min(wmaxlev(N,'db4'),floor(log2(N)))
K_norm = prm.numSTS/prm.numRx; % Нормировка по энергии
%% Параметры OFDM 
prm.numSC = 450; % Кол-во поднессущих
prm.N_FFT = 512; % Длина FFT для OFDM
prm.Nsymb_ofdm = 10; % Кол-во символов OFDM от каждой антенны
prm.CyclicPrefixLength = 64;  % длина защитных интервалов = 2*Ngi
prm.tmp_NCI = prm.N_FFT - prm.numSC;
prm.NullCarrierIndices = [1:prm.tmp_NCI/2 prm.N_FFT-prm.tmp_NCI/2+1:prm.N_FFT]'; % Guards and DC
%% Расчет
prm.numSTSVec = prm.numSTS;
prm.numUsers = 1; % const
prm.n = prm.bps*prm.Nsymb_ofdm*prm.numSC*prm.numSTS;% Длина бинарного потока
[isTxURA,prm.expFactorTx,isRxURA,prm.expFactorRx] = helperArrayInfo(prm,true);
%Возвращает 1 если можно использовать URA для Tx и Rx антенн 
%expFactor во сколько раз антенн больше чем потоков.
%% Координаты базовой/мобильной станций
prm.fc = 28e9;               % 28 GHz system несущая
prm.cLight = physconst('LightSpeed');
prm.lambda = prm.cLight/prm.fc;

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
% Нормировка
wT = wT./sqrt(sum(abs(wT).^2));
wR = wR./sqrt(sum(abs(wR).^2));
%% Параметры канала
prm.KFactor = 1;% Для 'RIC'
prm.SEED = 122;% Для 'RAYL_SPECIAL' 586 122 12   
prm.SampleRate = 40e6;
dt = 1/prm.SampleRate;
switch flag_chanel
    case "RAYL"       
        prm.tau = [2*dt 5*dt 7*dt];
        prm.pdB = [-3 -9 -12];
        % prm.tau = [2*dt 7*dt 15*dt];
        % prm.pdB = [-3 -9 -12]
    otherwise
        prm.tau = 5*dt;
        prm.pdB = -10;
end

%% ---------Сам скрипт--------
if flag_cor_MIMO == 2
    ostbcEnc = comm.OSTBCEncoder('NumTransmitAntennas',prm.numTx);
    ostbcComb = comm.OSTBCCombiner('NumReceiveAntennas',prm.numRx);
    prm.n = prm.n/prm.numTx;
end
SNR_MAX = 40;
SNR = 0+floor(10*log10(prm.bps)):SNR_MAX+floor(10*log10(prm.bps*prm.numTx));
prm.MinNumErr = 100; % Порог ошибок для цикла 
prm.conf_level = 0.95; % Уровень достоверности
prm.MAX_indLoop = 5;% Максимальное число итераций в цикле while
Koeff = 1/15;%Кол-во процентов от BER  7%
Exp = 1;% Кол-во опытов
for indExp = 1:Exp
    %% Создание канала
    [H,~,H_STS] = create_chanel(flag_chanel,prm);
    for indSNR = 1:length(SNR)
        berconf_M = 0;
        berconf_MM = 0;
        ErrNum_M = 0; % кол-во ошибок MIMO 
        ErrNum_MM = 0; % кол-во ошибок MIMO 
        indLoop = 0;  % индикатор итераций цикла while
        LenIntLoop_MM = 100;
        LenIntLoop_M = 100;
        condition_M = ((LenIntLoop_M > berconf_M*Koeff)||(ErrNum_M < prm.MinNumErr));
        condition_MM = ((LenIntLoop_MM > berconf_MM*Koeff)||(ErrNum_MM < prm.MinNumErr));
        while (condition_MM || condition_M) && (indLoop < prm.MAX_indLoop)
            %% Формируем данные
            Inp_data = randi([0 1],prm.n,1); % Передаваемые данные    
            Inp_Mat = reshape(Inp_data,length(Inp_data)/prm.bps,prm.bps); %Группируем биты 
            Inp_Sym = bi2de(Inp_Mat);  % Входные данные в символах
            %% Модулятор
            % MIMO
            Mod_data_inp_tmp = qammod(Inp_Sym,prm.M);% Модулятор QAM-M для полезной инф
            Mod_data_inp = reshape(Mod_data_inp_tmp,prm.numSC,prm.Nsymb_ofdm,prm.numSTS);
            % Модулятор пилотов  MIMO
            [preambula,ltfSC] = My_helperGenPreamble(prm);
            %% Зондирование канала
            %% Модулятор OFDM
            OFDM_data_STS = ofdmmod(Mod_data_inp,prm.N_FFT,prm.CyclicPrefixLength,...
                         prm.NullCarrierIndices);                      
            OFDM_data_STS = [preambula ; OFDM_data_STS];
            % Repeat over numTx
            Inp_dataMod = zeros(size(OFDM_data_STS,1),prm.numTx);
            for i = 1:prm.expFactorTx
                Inp_dataMod(:,(i-1)*prm.numSTS+(1:prm.numSTS)) = OFDM_data_STS;
            end
            %% Формируем луч на передачу
            if flag_Steering==1
                Inp_dataMod = Inp_dataMod.*conj(wT).';                   
            end
            Inp_dataMod_M = OFDM_data_STS;
            %% Прохождение канала
            switch flag_chanel
                case {'RAYL','RIC','RAYL_SPECIAL'}
    %                 H.Visualization = 'Impulse and frequency responses';
    %                 H.AntennaPairsToDisplay = [2,2];
    %                 H_siso.Visualization = 'Impulse and frequency responses';
                    [Chanel_data, H_ist] = H(Inp_dataMod);
                    [Chanel_data_M, H_ist_M] = H_STS(Inp_dataMod_M);
                otherwise                  
                    Chanel_data  = Inp_dataMod*H;
                    Chanel_data_M  = Inp_dataMod_M*H_STS; 

            end
            %% Собственный  шум
%             Noise_data = awgn(Chanel_data,SNR(indSNR),'measured');
%             Noise_data_M = awgn(Chanel_data_M,SNR(indSNR),'measured');
            [Noise_data,sigma] = my_awgn(Chanel_data,SNR(indSNR));%SNR(indSNR)
            [Noise_data_M,sigma_M] = my_awgn(Chanel_data_M,SNR(indSNR));%SNR(indSNR)
            %% Прием
            if flag_Steering==1
                Noise_data = Noise_data.*conj(wR).';                  
            end
            %% Демодулятор OFDM
            Mod_data_out = ofdmdemod(Noise_data,prm.N_FFT,prm.CyclicPrefixLength,prm.CyclicPrefixLength, ...
                prm.NullCarrierIndices);
            Mod_data_out_M = ofdmdemod(Noise_data_M,prm.N_FFT,prm.CyclicPrefixLength,prm.CyclicPrefixLength, ...
                prm.NullCarrierIndices);
            %% Оценка канала  
            H_estim = My_helperMIMOChannelEstimate(Mod_data_out(:,1:prm.numSTS,:),ltfSC,prm);
            H_estim_STS = My_helperMIMOChannelEstimate(Mod_data_out_M(:,1:prm.numSTS,:),ltfSC,prm);
            %% Вейвлет шумоподавление
            if flag_wav_MM == 1
                H_estim = H_WAV_my_mimo(H_estim,prm.LEVEL);
            end
            if flag_wav_MIMO == 1
                H_estim_STS = H_WAV_my_mimo(H_estim_STS,prm.LEVEL);
            end
            %% Эквалайзер
            %ZF Massiv MIMO
            if flag_cor_MM == 1
                Mod_data_out_ZF_tmp= My_MIMO_Equalize_ZF_numSC(Mod_data_out(:,prm.numSTS+1:end,:),H_estim);
                Mod_data_out_ZF = reshape(Mod_data_out_ZF_tmp,prm.numSC*prm.Nsymb_ofdm,prm.numSTS);
            else
                Mod_data_out_ZF_tmp = Mod_data_out(:,prm.numSTS+1:end,:);
                H_tmp = repmat(H_STS,prm.expFactorRx,1);
                Mod_data_out_ZF1 = reshape(Mod_data_out_ZF_tmp,prm.numSC*prm.Nsymb_ofdm,prm.numRx);
                Mod_data_out_ZF = Mod_data_out_ZF1*H_tmp/prm.expFactorRx;
            end
            %ZF MIMO
            if flag_cor_MIMO == 1
                Mod_data_out_ZF_tmp_M = My_MIMO_Equalize_ZF_numSC(Mod_data_out_M(:,prm.numSTS+1:end,:),H_estim_STS);
                Mod_data_out_ZF_M = reshape(Mod_data_out_ZF_tmp_M,prm.numSC*prm.Nsymb_ofdm,prm.numSTS);
            else
                Mod_data_out_ZF_tmp_M = Mod_data_out_M(:,prm.numSTS+1:end,:);
                Mod_data_out_ZF_M = reshape(Mod_data_out_ZF_tmp_M,prm.numSC*prm.Nsymb_ofdm,prm.numSTS);
            end
            %% Демодулятор
            Mod_data_out_tmp = Mod_data_out_ZF(:);
            Out_Sym = qamdemod(Mod_data_out_tmp,prm.M);    
            Mod_data_out_M_tmp = Mod_data_out_ZF_M(:);
            Out_Sym_M = qamdemod(Mod_data_out_M_tmp,prm.M);
            %% Выходные данные
            Out_Mat = de2bi(Out_Sym);
            Out_data = Out_Mat(:);
            Out_Mat_M = de2bi(Out_Sym_M);
            Out_data_M  = Out_Mat_M(:);  
            ErrNum_MM = ErrNum_MM+sum(abs(Out_data-Inp_data));          
            ErrNum_M = ErrNum_M+sum(abs(Out_data_M-Inp_data));
            %%
            indLoop = indLoop+1;
            [berconf_MM,conf_int_MM] = berconfint(ErrNum_MM,indLoop*length(Inp_data),prm.conf_level);
            [berconf_M,conf_int_M] = berconfint(ErrNum_M,indLoop*length(Inp_data),prm.conf_level);
            LenIntLoop_MM = conf_int_MM(2)-conf_int_MM(1);
            LenIntLoop_M = conf_int_M(2)-conf_int_M(1);
            condition_MM = ((LenIntLoop_MM > berconf_MM/15)||(ErrNum_MM < prm.MinNumErr));
            condition_M = ((LenIntLoop_M > berconf_M/15)||(ErrNum_M < prm.MinNumErr));
        end
        ber(indExp,indSNR) = berconf_MM;
        ber_M(indExp,indSNR) = berconf_M;
        if ErrNum_M>ErrNum_MM
            ErrNum_disp = ErrNum_M;
            name = 'Er_MIMO';
        else
            ErrNum_disp = ErrNum_MM;
            name = 'Er_MMIMO';
        end
        fprintf(['Complete %d db ' name ' = %d, ind = %d\n'],SNR(indSNR),ErrNum_disp,indLoop);
    end
    fprintf('Exp %d  \n',indExp);
end
prm.bps = prm.bps*prm.numSTS;
ber_mean = mean(ber,1);
ber_M_mean = mean(ber_M,1);
Eb_N0_M = SNR(1:size(ber_mean,2))-(10*log10(prm.bps));
Eb_N0 = 0:60;
ther_ber = berawgn(Eb_N0,'qam',16);
figure() 
plot_ber(ther_ber,Eb_N0,1,'g',1.5,0)
plot_ber(ber_mean,SNR(1:size(ber_mean,2)),prm.bps,'k',1.5,0)
plot_ber(ber_M_mean,SNR(1:size(ber_M_mean,2)),prm.bps,'b',1.5,0)
str1 = ['MASSIV MIMO ' num2str(prm.numTx) 'x'  num2str(prm.numRx)];
str2 = ['MIMO ' num2str(prm.numSTS) 'x'  num2str(prm.numSTS)];
legend('Теоретическая',str1,str2);
% save(str,'ber_mean','ber_siso_mean','SNR','prm','ber','ber_siso')