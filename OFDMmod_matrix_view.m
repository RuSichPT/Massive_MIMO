function [dataOutput] = OFDMmod_matrix_view(dataInput,Nsymb,Msc,Ngi,N)
% dataInput [1 x поток]    
dataOutput = [];
    for i = 0:Nsymb-1
        Spectr_OFDM = zeros(1,N);
        Spectr_OFDM(1:Msc/2) = dataInput(i*Msc+1:i*Msc + Msc/2);
        Spectr_OFDM(N - Msc/2 + 1 :N) = dataInput(i*Msc + Msc/2+1:(i+1)*Msc);
        dataOFDM_signal = ifft(Spectr_OFDM,N);
        OFDM_signal_guard = [dataOFDM_signal(1:Ngi) dataOFDM_signal dataOFDM_signal(N-9:N)];
        dataOutput = [dataOutput; OFDM_signal_guard];
    end
end

