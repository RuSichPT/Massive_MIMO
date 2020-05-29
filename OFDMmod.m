function [dataOutput] = OFDMmod(dataInput,Nsymb,Msc,Ngi,N)
    dataOutput = [];
    for i = 0:Nsymb-1
        Spectr_OFDM = zeros(1,N); 
        Spectr_OFDM(1:Msc/2) = dataInput(i*Msc+1:i*Msc + Msc/2); % Помещаем символы на поднесущие левая часть
        Spectr_OFDM(N - Msc/2 + 1 :N) = dataInput(i*Msc + Msc/2+1:(i+1)*Msc);% Помещаем символы на поднесущие правая часть
        dataOFDM_signal = ifft(Spectr_OFDM,N);% Получаем временную реализацию
        OFDM_signal_guard = [dataOFDM_signal(1:Ngi) dataOFDM_signal dataOFDM_signal(N-9:N)];% Защитные интервалы
        dataOutput = [dataOutput OFDM_signal_guard];
    end
end

