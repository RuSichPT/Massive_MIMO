%% ---------Модель Massiv MIMO and MIMO-------- 
clear;clc;%close all;
rng(29)
%% Управление
flag_chanel_ber = 1; % Enable/disable chanel ber
flag_constel = 0;  % Enable/disable
flag_DN = 0; % Enable/disable построение ДН
flag_preCod = 1; % Enable/disable прекодирование не дает улучшений
flag_Steering = 0;       % Enable/disable steering не дает улучшений
flag_chanel = 'Scattering'; % 'AWGN' ,'RAYL','RIC','RAYL_SPECIAL','STATIC', 'BAD', 'Scattering' 'ScatteringFlat'
flag_cor_MIMO = 1; % 1-коррекция АЧХ (эквалайзер для MIMO)
flag_cor_MM = 1; % 1-коррекция АЧХ (эквалайзер для MAS MIMO)
flag_wav_MIMO = 0; % вейвлет шумоподавление для MIMO не работает при  prm.nRays>10
flag_wav_MM = 0; % вейвлет шумоподавление для MAS MIMO не работает при  prm.nRays>10
%% Параметры системы 
prm.numTx = 16; % Кол-во излучающих антен 
prm.numRx = 4; % Кол-во приемных антен
prm.numSTS = 4; % Кол-во потоков 2/4/8/16/32/64
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
% Account for channel filter delay
numPadSym = 3;          % number of symbols to zeropad
prm.numPadZeros = numPadSym*(prm.N_FFT+prm.CyclicPrefixLength); 
%% Расчет
prm.numSTSVec = prm.numSTS;
prm.numUsers = 1; % const
prm.n = prm.bps*prm.Nsymb_ofdm*prm.numSC;% Длина бинарного потока
%% Координаты базовой/мобильной станций
prm.fc = 28e9;               % 28 GHz system несущая
prm.fc_M = 28e9;
prm.cLight = physconst('LightSpeed');
prm.lambda = prm.cLight/prm.fc;
prm.lambda_M = prm.cLight/prm.fc_M;

prm.posTx = [0;0;0];         % BS/Transmit array position, [x;y;z], meters
% maxRange = 1000;            % all MSs within 1000 meters of BS
% prm.mobileRange = randi([1 maxRange],1,prm.numUsers);
prm.mobileRange = 300;
% Angles specified as [azimuth;elevation], az=[-180 180], el=[-90 90] elevation - угол места
% prm.mobileAngle = [rand(1,prm.numUsers)*360-180; ... в градусах
%                     rand(1,prm.numUsers)*180-90];
prm.mobileAngle = [33; 22];
prm.steeringAngle_Tx = [33; 22];%prm.mobileAngle;% Transmit steering angle (Близкие к mobileAngle)
[xRx,yRx,zRx] = sph2cart(deg2rad(prm.mobileAngle(1)),...
            deg2rad(prm.mobileAngle(2)),prm.mobileRange);
prm.posRx = [xRx;yRx;zRx];

[toRxRange,prm.steeringAngle_Rx] = rangeangle(prm.posTx,prm.posRx);
spLoss_db = fspl(toRxRange,prm.lambda); % db toRxRange = prm.mobileRange
spLoss_M_db = fspl(toRxRange,prm.lambda_M); % db toRxRange = prm.mobileRange
%% Весовые коэф для прд/прм ФАР            
[isTxURA,prm.expFactorTx,isRxURA,prm.expFactorRx] = helperArrayInfo(prm,true);
isTxURA_M = 0;isRxURA_M = 0; prm.expFactorTx_M = 1;  prm.expFactorRx_M = 1;
%Возвращает 1 если можно использовать URA для Tx и Rx антенн 
%expFactor во сколько раз антенн больше чем потоков.
[wT, prm.arrayTx] = Transmit_Beam_Steering(prm,prm.numTx,prm.fc,prm.lambda,prm.expFactorTx,isTxURA,flag_DN);
[wR, prm.arrayRx] = Receive_Beam_Steering(prm,prm.numRx,prm.fc,prm.lambda,prm.expFactorRx,isRxURA,flag_DN);
[wT_M, prm.arrayTx_M] = Transmit_Beam_Steering(prm,prm.numSTS,prm.fc_M,prm.lambda_M,prm.expFactorTx_M,isTxURA_M,flag_DN);
[wR_M, prm.arrayRx_M] = Receive_Beam_Steering(prm,prm.numSTS,prm.fc_M,prm.lambda_M,prm.expFactorRx_M,isRxURA_M,flag_DN);
prm.posTxElem = getElementPosition(prm.arrayTx)/prm.lambda;
prm.posRxElem = getElementPosition(prm.arrayRx)/prm.lambda;
prm.posTxElem_M = getElementPosition(prm.arrayTx_M)/prm.lambda_M;
prm.posRxElem_M = getElementPosition(prm.arrayRx_M)/prm.lambda_M;
%% Параметры канала
prm.nRays = 10; % Для 'Scattering'   
prm.KFactor = 1;% Для 'RIC'
prm.SEED = 86;% Для 'RAYL_SPECIAL'  'Scattering'    586 122 12   
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
SNR_MAX = 70;
SNR = 0+floor(10*log10(prm.bps)):SNR_MAX+floor(10*log10(prm.bps*prm.numTx));
if flag_constel == 1
    SNR = 100;
end
numSTS = prm.numSTS;
prm.MinNumErr = 100; % Порог ошибок для цикла 
prm.conf_level = 0.95; % Уровень достоверности
prm.MAX_indLoop = 1;% Максимальное число итераций в цикле while
prm.MaxNumZero = 2; %  max кол-во нулевых точек в цикле while
Koeff = 1/15;%Кол-во процентов от BER  7%
Exp = 1;% Кол-во опытов
for indExp = 1:Exp
    %% Создание канала
    [H,~,H_STS] = create_chanel(flag_chanel,prm);
    NumZero = 0; % кол-во нулевых точек
    for indSNR = 1:length(SNR)
        berconf_M = 0;
        berconf_MM = 0;
        ErrNum_M = 0; % кол-во ошибок MIMO
        if flag_chanel_ber ==1
            ErrNum_MM = zeros(1,prm.numSTS);% кол-во ошибок  M MIMO
        else
            ErrNum_MM = 0; % кол-во ошибок  M MIMO
        end
        indLoop = 0;  % индикатор итераций цикла while
        LenIntLoop_MM = 100;
        LenIntLoop_M = 100;
        condition_M = ((LenIntLoop_M > berconf_M*Koeff)||(ErrNum_M < prm.MinNumErr));
        condition_MM = ((LenIntLoop_MM > berconf_MM*Koeff)||(ErrNum_MM < prm.MinNumErr));
        if (NumZero >= prm.MaxNumZero)
            break;
        end
        while (condition_MM || condition_M) && (indLoop < prm.MAX_indLoop)
            %% Зондирование канала
            if flag_preCod == 1
                prm.numSTS = prm.numTx;
                [preambulaSTS,ltfSC_zond] = My_helperGenPreamble(prm);
                %Прохождение канала
                % добавляем интервал для задержки
                preambulaZond = [preambulaSTS ; zeros(prm.numPadZeros,prm.numTx)];
                switch flag_chanel
                    case {'RAYL','RIC','RAYL_SPECIAL'}
%                         H.Visualization = 'Impulse and frequency responses';
                        [Chanel_Zond, H_ist_Z] = H(preambulaZond);                       
                        chanDelay = 0;
                    case 'Scattering'                  
                        [Chanel_Zond, H_ist,tau] = H(preambulaZond);
                        chanDelay = floor(min(tau)*prm.SampleRate);
                    otherwise                  
                        Chanel_Zond  = preambulaZond*H;
                        chanDelay = 0;
                end
                %Собственный  шум
                Noise_Zond_tmp = awgn(Chanel_Zond,SNR(indSNR),'measured');%SNR(indSNR)
                %Демодулятор OFDM
                Noise_Zond = Noise_Zond_tmp(chanDelay+1:end-(prm.numPadZeros-chanDelay),:); % Убираем задержку
                Zond_out = ofdmdemod(Noise_Zond,prm.N_FFT,prm.CyclicPrefixLength,prm.CyclicPrefixLength, ...
                    prm.NullCarrierIndices);
                %Оценка канала  
                H_estim_zond = My_helperMIMOChannelEstimate(Zond_out,ltfSC_zond,prm);
                prm.numSTS =  prm.numSTS/prm.expFactorTx;
                [wP,wC] = diagbfweights(H_estim_zond);              
%                 squeeze(wP(1,:,:))*squeeze(H_estim_zond(1,:,:))*squeeze(wC(1,:,:))
            end
            %% Формируем данные
            Inp_data = randi([0 1],prm.n,numSTS); % Передаваемые данные
            Inp_data_tmp = reshape(Inp_data,prm.n*numSTS,1);
            %% Модулятор
            % MIMO
            Mod_data_inp_tmp = qammod(Inp_data,prm.M,'InputType','bit');% Модулятор QAM-M для полезной инф
            Mod_data_inp = reshape(Mod_data_inp_tmp,prm.numSC,prm.Nsymb_ofdm,numSTS);
            Mod_data_inp_M = Mod_data_inp;
            if flag_preCod ==1
                    for j = 1:prm.numSC
                        Q = squeeze(wP(j,1:numSTS,:));
                        Mod_data_inp_Tx(j,:,:) = squeeze(Mod_data_inp(j,:,:))*Q;       
                    end
                    Mod_data_inp = Mod_data_inp_Tx;
                % Модулятор пилотов  MIMO
                [preambula,ltfSC] = My_helperGenPreamble(prm,wP(:,1:numSTS,:));
            else
                % Модулятор пилотов  MIMO
                [preambula,ltfSC] = My_helperGenPreamble(prm);
            end
            %% Модулятор OFDM
            OFDM_data_STS = ofdmmod(Mod_data_inp,prm.N_FFT,prm.CyclicPrefixLength,...
                         prm.NullCarrierIndices);                            
            OFDM_data_STS = [preambula ; OFDM_data_STS];
            OFDM_data_STS = 1/sqrt(prm.expFactorTx)*OFDM_data_STS; % Нормируем мощность
            
            [preambula_M,ltfSC_M] = My_helperGenPreamble(prm);
            OFDM_data_STS_M = ofdmmod(Mod_data_inp_M,prm.N_FFT,prm.CyclicPrefixLength,...
                         prm.NullCarrierIndices);                      
            OFDM_data_STS_M = [preambula_M ; OFDM_data_STS_M];
            
            Inp_dataMod = OFDM_data_STS;
            Inp_dataMod_M = OFDM_data_STS_M;
            %% Формируем луч на передачу
            if flag_Steering==1
                Inp_dataMod = Inp_dataMod.*conj(wT).';
            end           
            %% Прохождение канала
            % добавляем интервал для задержки
            Inp_dataMod_tmp = [Inp_dataMod ; zeros(prm.numPadZeros,prm.numTx)];
            Inp_dataMod_tmp_M = [Inp_dataMod_M ; zeros(prm.numPadZeros,prm.numSTS)];
            switch flag_chanel
                case {'RAYL','RIC','RAYL_SPECIAL'}
    %                 H.Visualization = 'Impulse and frequency responses';
    %                 H.AntennaPairsToDisplay = [2,2];
    %                 H_siso.Visualization = 'Impulse and frequency responses';
                    [Chanel_data, H_ist] = H(Inp_dataMod_tmp);
                    [Chanel_data_M, H_ist_M] = H_STS(Inp_dataMod_tmp_M);
                    chanDelay = 0;
                    chanDelay_M = 0;
                case 'Scattering'                  
                    [Chanel_data, H_ist,tau] = H(Inp_dataMod_tmp);
                    [Chanel_data_M, H_ist_M,tau_M] = H_STS(Inp_dataMod_tmp_M);
                    chanDelay = floor(min(tau)*prm.SampleRate); 
                    chanDelay_M = floor(min(tau_M)*prm.SampleRate);
                otherwise                  
                    Chanel_data  = Inp_dataMod_tmp*H;
                    Chanel_data_M  = Inp_dataMod_tmp_M*H_STS;
                    chanDelay = 0;
                    chanDelay_M = 0;
            end
            %% Прием
            if flag_Steering==1
                Chanel_data = Chanel_data.*conj(wR).';
            end
            %% Собственный  шум
            Noise_data = awgn(Chanel_data,SNR(indSNR),'measured');
            Noise_data_M = awgn(Chanel_data_M,SNR(indSNR),'measured');
%             [Noise_data,sigma] = my_awgn(Chanel_data,SNR(indSNR));%SNR(indSNR)
%             [Noise_data_M,sigma_M] = my_awgn(Chanel_data_M,SNR(indSNR));%SNR(indSNR)

            % Убираем задержку и усиливаем
            Gain = db2pow(spLoss_db);
            Gain_M = db2pow(spLoss_M_db);
            Noise_data_tmp = Noise_data(chanDelay+1:end-(prm.numPadZeros-chanDelay),:)*Gain;
            Noise_data_tmp_M = Noise_data_M(chanDelay_M+1:end-(prm.numPadZeros-chanDelay_M),:)*Gain_M;          
            %% Демодулятор OFDM
            Mod_data_out = ofdmdemod(Noise_data_tmp,prm.N_FFT,prm.CyclicPrefixLength,prm.CyclicPrefixLength, ...
                prm.NullCarrierIndices);
            Mod_data_out_M = ofdmdemod(Noise_data_tmp_M,prm.N_FFT,prm.CyclicPrefixLength,prm.CyclicPrefixLength, ...
                prm.NullCarrierIndices);
%             if flag_preCod ==1
%                 for j = 1:prm.numSC
%                     Q = squeeze(wC(j,:,:));
%                     Mod_data_out(j,:,:) = squeeze(Mod_data_out(j,:,:))*Q;
%                 end
%             end
            %% Оценка канала  
            H_estim = My_helperMIMOChannelEstimate(Mod_data_out(:,1:prm.numSTS,:),ltfSC,prm);
            H_estim_STS = My_helperMIMOChannelEstimate(Mod_data_out_M(:,1:prm.numSTS,:),ltfSC_M,prm);      
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
                Mod_data_out_ZF_tmp = My_MIMO_Equalize_ZF_numSC(Mod_data_out(:,prm.numSTS+1:end,:),H_estim);
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
            if flag_constel == 1
                scatterplot(Mod_data_out_ZF(:));
            end
            Out_data = qamdemod(Mod_data_out_ZF,prm.M,'OutputType','bit');
            Out_data_tmp = reshape(Out_data,prm.n*prm.numSTS,1);
            if flag_constel == 1
                scatterplot(Mod_data_out_ZF_M(:));
            end
            Out_data_M = qamdemod(Mod_data_out_ZF_M,prm.M,'OutputType','bit');
            Out_data_M_tmp  = reshape(Out_data_M,prm.n*prm.numSTS,1);
            %% Выходные данные
            if flag_chanel_ber ==1
                for m = 1:prm.numSTS
                    ErrNum_MM(m) = ErrNum_MM(m)+sum(abs(Out_data(:,m)-Inp_data(:,m)));
                end
            else
                ErrNum_MM = ErrNum_MM+sum(abs(Out_data_tmp-Inp_data_tmp)); 
            end
            ErrNum_M = ErrNum_M+sum(abs(Out_data_M_tmp-Inp_data_tmp));
            %%
            indLoop = indLoop+1;
            if flag_chanel_ber ==1
                for m = 1:prm.numSTS
                    [berconf_MM(m),conf_int_MM(m,:)] = berconfint(ErrNum_MM(m),indLoop*length(Inp_data(:,m)),prm.conf_level);
                    LenIntLoop_MM(m) = conf_int_MM(m,2)-conf_int_MM(m,1);
                    condition_MM(m) = ((LenIntLoop_MM(m)  > berconf_MM(m)/15)||(ErrNum_MM(m) < prm.MinNumErr));
                end
                condition_MM = max(condition_MM);
                ErrNum_MM_max = max(ErrNum_MM);
            else
                [berconf_MM,conf_int_MM] = berconfint(ErrNum_MM,indLoop*length(Inp_data_tmp),prm.conf_level);
                LenIntLoop_MM = conf_int_MM(2)-conf_int_MM(1);
                condition_MM = ((LenIntLoop_MM > berconf_MM/15)||(ErrNum_MM < prm.MinNumErr));
                ErrNum_MM_max = ErrNum_MM;
            end
            [berconf_M,conf_int_M] = berconfint(ErrNum_M,indLoop*length(Inp_data_tmp),prm.conf_level);
            LenIntLoop_M = conf_int_M(2)-conf_int_M(1);
            condition_M = ((LenIntLoop_M > berconf_M/15)||(ErrNum_M < prm.MinNumErr));
        end
        if (ErrNum_MM_max == 0)&&(ErrNum_M==0)
            NumZero = NumZero+1;
        end
        if flag_chanel_ber ==1
            for m = 1:prm.numSTS
                ber{m}(indExp,indSNR) = berconf_MM(m);
            end
        else
            ber(indExp,indSNR) = berconf_MM;
        end
        ber_M(indExp,indSNR) = berconf_M;
        if ErrNum_M>ErrNum_MM_max
            ErrNum_disp = ErrNum_M;
            name = 'Er_MIMO';
        else
            ErrNum_disp = ErrNum_MM_max;
            name = 'Er_MMIMO';
        end
        fprintf(['Complete %d db ' name ' = %d, ind = %d NZ = %d\n'],...
            SNR(indSNR),ErrNum_disp,indLoop,NumZero);
    end
    fprintf('Exp %d  \n',indExp);
end
if flag_chanel_ber ==1
    for m = 1:prm.numSTS
        ber_mean(m,:) = mean(ber{m},1);
    end
else
    ber_mean = mean(ber,1);
end
ber_M_mean = mean(ber_M,1);
SNR = SNR(1:max(size(ber_M_mean,2),size(ber_mean,2)));
Eb_N0_M = SNR -(10*log10(prm.bps*prm.numSTS));
Eb_N0 = 0:60;
ther_ber = berfading(Eb_N0,'qam',256,1);
% ther_ber1 = berawgn(Eb_N0,'qam',4);
figure() 
% plot_ber(ther_ber,Eb_N0,1,'g',1.5,0)
% plot_ber(ther_ber1,Eb_N0,1,'r',1.5,0)
if flag_chanel_ber ==1
    color = {'r';'b';'g';'c';'m';'y';'--r';'--b';'--g';'--c';'--m';'--y';':r';':b';':g';':c';':m';':y'};
%     color = {'r';'b';'g';'c'};
    for m = 1:prm.numSTS
        plot_ber(ber_mean(m,:),SNR(1:size(ber_mean,2)),prm.bps,color{m},1.5,0)
    end
    mean_bear = mean(ber_mean,1);
    plot_ber(mean_bear,SNR(1:size(ber_mean,2)),prm.bps*prm.numSTS,'k',1.5,0)
%     legend('1','2','3','4','mean');
else
    plot_ber(ber_mean,SNR(1:size(ber_mean,2)),prm.bps*prm.numSTS,'k',1.5,0)
end
plot_ber(ber_M_mean,SNR(1:size(ber_M_mean,2)),prm.bps*prm.numSTS,'--k',3,0)
% str1 = ['MASSIV MIMO ' num2str(prm.numTx) 'x'  num2str(prm.numRx)];
% str2 = ['MIMO ' num2str(prm.numSTS) 'x'  num2str(prm.numSTS)];
% legend(str1,str2);
str = ['DataBase/' num2str(prm.numTx) 'x' num2str(prm.numRx) 'x' num2str(prm.numSTS)  '_pre ='...
    num2str(flag_preCod) '_ster =' num2str(flag_Steering) '_' flag_chanel '_Wmm=' ...
    num2str(flag_wav_MM) '_Wm=' num2str(flag_wav_MIMO) '_Exp=' num2str(Exp) '.mat'];
% save(str,'ber_mean','ber_M_mean','SNR','prm','ber','ber_M')