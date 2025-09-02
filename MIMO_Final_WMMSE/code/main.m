%% (a)
clear;
clc;
close all;

K = 4; % BS
T = 1; % transmit antennas
R = 1; % receive antennas
epsilon = 1e-3;
max_iter = 100;
sigma2 = 1; % noise power
I = 4; % users per BS
alpha1 = ones(I,K); % weight
d = 4; % data stream
d = min([d, T, R]);
snr = 30;
P = db2pow(snr)*sigma2; % transimt power
H = cell(I,K,K); % channel
for i=1:I
    for k = 1:K
        for j=1:K
           H{i,k,j}=sqrt(1/2)*(randn(R,T)+1i*randn(R,T)); % H{i,k,j}=(R * T)
        end
    end
end
rate = []; 

% U = cell(I,K);
% U(:)={zeros(R,d)};

V = cell(I,K); 
for i=1:I
    for k=1:K
        v = randn(T,d)+1i*randn(T,d);
        V{i,k}=sqrt(P/I)*v/norm(v,"fro");
    end
end 

rate_old = sum_rate(H,V,sigma2,R,I,K,alpha1);
rate = [rate rate_old];

iter1 = 1;
while(1)
    U = find_U(H,V,sigma2,R,I,K,d); 
    W = find_W(U,H,V,I,K,d); 
    V = find_V(alpha1,H,U,W,T,I,K,P); 
    rate_new = sum_rate(H,V,sigma2,R,I,K,alpha1);
    rate = [rate rate_new];
    iter1 = iter1 + 1;
    if abs(rate_new-rate_old) / rate_old < epsilon || iter1 > max_iter
        break;
    end
    rate_old = rate_new;
end


plot(0:iter1-1,rate,'r-o')
grid on
xlabel('Iterations')
ylabel('Sum rate (bits per channel use)')
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1)
% title('MIMO-IFC, SNR=30, K=4, T=4, R=2, \epsilon=1e-3','Interpreter','tex')
title('SISO-IFC, SNR=30, K=4, T=1, R=1, \epsilon=1e-3','Interpreter','tex')

%% (b)
clear;
clc;
close all;

K = 3; % BS
T = 1; % transmit antennas
R = 1; % receive antennas
epsilon = 1e-3;
max_iter = 100;
sigma2 = 1; % noise power
I = 2; % users per BS
alpha1 = ones(I,K); % weight
d = 4; % data stream
d = min([d, T, R]);

snr_dB = 5:5:30;
rate_snr = zeros(size(snr_dB));
rate_snr_int10 = zeros(size(snr_dB));
num_realizations = 100;
for s = 1:length(snr_dB)
    snr = snr_dB(s);
    P = db2pow(snr)*sigma2; % transimt power
    total_rate = 0;
    total_rate_int10 = 0;

    for h = 1:num_realizations
        H = cell(I,K,K); % channel
        for i=1:I
            for k = 1:K
                for j=1:K
                   H{i,k,j}=sqrt(1/2)*(randn(R,T)+1i*randn(R,T)); % H{i,k,j}=(R * T)
                end
            end
        end
        
        % WMMSE_once ---------------------------------------------------------
        rate = []; 
        V = cell(I,K); 
        for i=1:I
            for k=1:K
                v = randn(T,d)+1i*randn(T,d);
                V{i,k}=sqrt(P/I)*v/norm(v,"fro");
            end
        end 
        
        rate_old = sum_rate(H,V,sigma2,R,I,K,alpha1);
        rate = [rate rate_old];
        
        iter1 = 1;
        while(1)
            U = find_U(H,V,sigma2,R,I,K,d); 
            W = find_W(U,H,V,I,K,d); 
            V = find_V(alpha1,H,U,W,T,I,K,P); 
            rate_new = sum_rate(H,V,sigma2,R,I,K,alpha1);
            rate = [rate rate_new];
            iter1 = iter1 + 1;
            if abs(rate_new-rate_old) / rate_old < epsilon || iter1 > max_iter
                break;
            end
            rate_old = rate_new;
        end

        total_rate = total_rate + rate(end);
        % WMMSE_once ---------------------------------------------------------

        % WMMSE_10 ---------------------------------------------------------
        best_rate = -inf;
        for t = 1:10
            rate = []; 
            V = cell(I,K); 
            for i=1:I
                for k=1:K
                    v = randn(T,d)+1i*randn(T,d);
                    V{i,k}=sqrt(P/I)*v/norm(v,"fro");
                end
            end 
            
            rate_old = sum_rate(H,V,sigma2,R,I,K,alpha1);
            rate = [rate rate_old];
            
            iter1 = 1;
            while(1)
                U = find_U(H,V,sigma2,R,I,K,d); 
                W = find_W(U,H,V,I,K,d); 
                V = find_V(alpha1,H,U,W,T,I,K,P); 
                rate_new = sum_rate(H,V,sigma2,R,I,K,alpha1);
                rate = [rate rate_new];
                iter1 = iter1 + 1;
                if abs(rate_new-rate_old) / rate_old < epsilon || iter1 > max_iter
                    break;
                end
                rate_old = rate_new;
            end

            if rate(end) > best_rate
                best_rate = rate(end);
            end
        end

        total_rate_int10 = total_rate_int10 + best_rate;
        % WMMSE_10 ---------------------------------------------------------

    end

    rate_snr(s) = total_rate / num_realizations;
    rate_snr_int10(s) = total_rate_int10 / num_realizations;
end



figure;
plot(snr_dB,rate_snr_int10,'r-o')
grid on
xlabel('SNR')
ylabel('Average sum rate (bits per channel use)')
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1)
% title('MIMO-IFC, K=3, T=2, R=2, \epsilon=1e-3','Interpreter','tex')
title('SISO-IFC, K=3, T=1, R=1, \epsilon=1e-3','Interpreter','tex')

hold on
plot(snr_dB,rate_snr,'b-*')
legend('WMMSE int10','WMMSE')