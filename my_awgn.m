function [Y,sigma] = my_awgn(X,snr)
% X - [поток, numRx]
for i = 1:size(X,2)
    E_tmp(i) = sum(abs(X(:,i)).^2); % Энергия       
end
E = sum(E_tmp);
P = E/(size(X,1)*size(X,2));% Средняя Мощность
% P = E/size(X,1);% Средняя Мощность 
A = sqrt(P); % Средние амплитуды
sigma = A*10^((-snr)/20); % СКО
Y =  X+0.707 * sigma * (randn(size(X,1),size(X,2))+ 1i * randn(size(X,1),size(X,2)));
end

