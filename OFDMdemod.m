function [dataOutput,dataInput] = OFDMdemod(dataInput,Nsymb,Msc,Ngi,N,Freq_Char_Est)
% dataInput [1 x поток]
    dataOutput = [];
    for i = 0:Nsymb-1
        data_noisy = dataInput(1:N+2*Ngi);
        dataInput(1:N+2*Ngi) = [];
        data_noisy(1:Ngi)= []; % удаляем защитные интервалы
        data_noisy(end-Ngi+1:end)= [];% удаляем защитные интервалы
        Spectr_noisy = fft(data_noisy,N);
        dataMod_noisy = [Spectr_noisy(1:Msc/2) Spectr_noisy(N-Msc/2 + 1: N)]./ Freq_Char_Est;
        dataOutput = [dataOutput dataMod_noisy];
    end
end

