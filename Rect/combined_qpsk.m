clc
clear all
bits=10000;
data=randi(1,bits)>0.5;
%---debugging---
%data=[1 1 1]
%xxxxxxxxxx
SNR = -20:2:5;
BER=zeros(1,length(SNR));

colors={'b-*','k','r<'}; 
ind=1;
  
    %---Transmitter---------
    %Gray mapping of bits into symbols
    col=length(data)/2;
    I=zeros(1,col);
    Q=I;
   
    I=data(1:2:bits-1);
    Q=data(2:2:bits);
   
    I= -2.*I+1;
    Q= -2.*Q+1;
   
    symb=I+j.*Q;
   
           
            %----Filter
    psf=ones(1,1);
            %----
    M=length(psf);
for i=1:length(SNR)
            % inserting zeros between the bits
            % w.r.t number of coefficients of
            % PSF to pass the bit stream from the PSF
z=zeros(M-1,bits/2);

    upsamp=[symb;z];
    upsamp2=reshape(upsamp,1,(M)*bits/2);

    %Passing the symbols from PSF
    %tx_symb=conv(real(upsamp2),psf)+j*conv(imag(upsamp2),psf);
   
    tx_symb=conv(upsamp2,psf);
    %--------CHANNEL-----------
    %Random noise generation and addition to the signal
    npsd=10.^(SNR(i)/10);
    n_var=1/sqrt(2.*npsd);
    rx_symb=tx_symb+(n_var*randn(1,length(tx_symb))  +j*n_var*randn(1,length(tx_symb)) );
    %xxxxxxxxxxxxxxxxxxxxxxxxxx
   
    %-------RECEIVER-----------
    rx_match=conv(rx_symb,psf);   
    rx=rx_match(M:M:length(rx_match));
    rx=rx(1:1:bits/2);
    recv_bits=zeros(1,bits);
    %demapping
    k=1;
    for ii=1:bits/2
        recv_bits(k)=  -( sign(  real(  rx(ii)  )  )  -1)/2;
        recv_bits(k+1)=-( sign(  imag(  rx(ii)  )  )  -1)/2;
        k=k+2;
    end
       
       %sign(   real( rx )   )
       %sign(  imag(  rx )   )
        %data
        %tx_symb
        %rx_symb
       
        %recv_bits
   %xxxxxxxxxxxxxxxxxxxxxxxxxxx
   
   %---SIMULATED BIT ERROR RATE----
    errors=find(xor(recv_bits,data));   
    errors=size(errors,2);
    BER(i)=errors/bits;
    %xxxxxxxxxxxxxxxxxxxxxxxxxxx
end

fs=1;
n_pt=2^9;
tx_spec=fft(tx_symb,n_pt);
f= -fs/2:fs/n_pt:fs/2-fs/n_pt;
figure
plot(f,abs(fftshift(tx_spec)));
title('Signal Spectrum for Signal with Rectangular Pulse Shaping for QPSK');
xlabel('Frequency [Hz]');
ylabel('x(F)');
figure(1)


figure(2);
 plotHandle=plot(SNR,BER,char(colors(ind))); 
 set(plotHandle,'LineWidth',1.5); 
 hold on; 
 ind=ind+1;
 


 %%%%RC
T = 1/10e6;                 % Symbol time interval [s] 
r = 4;                      % Oversampling factor (r samples per pulse) 
N_symbols_per_pulse = 40;   % Duration of TX/RX filter in symbols 
alpha = 0.20;               % Roll-off factor (excess bandwidth) 

%Based on above we can define sampling frequency and sampling time interval as 
Fs = r/T;               % Sampling frequency 
Ts = 1/Fs;              % Sampling time interval 

%GENERATION OF BIT SEQUENCE AND QAM SYMBOLS 
N_symbols = 1000000;      % Number of symbols 

% Alphabet size 
M = 4;  % Number of symbols in the QAM alphabet (e.g. 16 means 16-QAM).

% Here qam_axis presents the symbol values in real/imaginary axis. 
qam_axis = -sqrt(M)+1:2:sqrt(M)-1; 
% qam_axis = [-7 -5 -3 -1 1 3 5 7]; % for 64-QAM 

% generation of a complex constellation: 
alphabet = bsxfun(@plus,qam_axis',1j*qam_axis); %help bsxfun 

% equivalent to alphabet = repmat(qam_axis', 1, sqrt(alphabet_size)) + ... 
%                          repmat(1j*qam_axis, sqrt(alphabet_size), 1);  
alphabet = alphabet(:).'; % alphabet symbols as a row vector   

% Scaling the constellation, so that the mean power of a transmitted symbol 
% is one (e.g., with QPSK this is 1/sqrt(2), and for 16-QAM 1/sqrt(10)) 
alphabet_scaling_factor = 1/sqrt(mean(abs(alphabet).^2));
alphabet = alphabet*alphabet_scaling_factor; 

% Number of bits, defined for a fixed alphabet size and number of symbols 
N_bits = log2(M)*N_symbols;  

% Random bit sequence 
bits = randi(2,N_bits,1)-1; 


% Block of bits in the columns 
B = reshape(bits,log2(M),[]); 

% Powers of two: ..., 2^3, 2^2, 2^1, 2^0 
q = 2.^(log2(M)-1:-1:0);    
 
% Symbol indices between 0...M-1, i.e. bit blocks to decimal numbers 
symbol_indices = q*B; 

% Gray coded symbol indices. This is basically an operation where the 
% symbol indices are mapped to Gray-coded symbol indices. Another option 
% would be reordering the original symbol alphabet, but the overall effects  
% would be exactly the same. 

[Gray_symbol_indices, ~] = bin2gray(symbol_indices, 'qpsk', M); 

symbols = alphabet(Gray_symbol_indices+1); 


%TRANSMITTER STRUCTURE 
%Implement the transit filter: Root-Raised-Cosine (RRC) and plot the pulse shape 
% Filter generation 
gt = rcosdesign(alpha,N_symbols_per_pulse,r,'normal'); 
gt1 = rcosdesign(alpha,N_symbols_per_pulse,r,'sqrt'); 

% Plot the pulse shape of the transmit/receive filter 
% figure 
% plot(-N_symbols_per_pulse*r/2*Ts:Ts:N_symbols_per_pulse*r/2*Ts,gt,'b') 
% hold on 
% stem(-N_symbols_per_pulse*r/2*Ts:T:N_symbols_per_pulse*r/2*Ts,gt(1:r:end),'ro') 
% xlabel('time [s]') 
% ylabel('Amplitude') 
% title('Transmit/receive RC filter (pulse shape)') 
% legend('Pulse shape','Ideal symbol-sampling locations') 


% Zero vector initilized for up-sampled symbol sequence 
symbols_upsampled = zeros(size(1:r*N_symbols));       
% symbol insertion 
symbols_upsampled(1:r: r*N_symbols) = symbols; 
% now the up-sampled sequence looks like {a1 0 0... a2 0 0... a3 0 0...}   
st = filter(gt,1,symbols_upsampled); % Transmitter filtering 
st = st(1+(length(gt)-1)/2:end);     % Filter delay correction 

st1 = filter(gt1,1,symbols_upsampled); % Transmitter filtering 
st1 = st1(1+(length(gt1)-1)/2:end);     % Filter delay correction 


%Plot the transmit signal s(t) in time and frequency domain 
% figure % zoom manually to see the signal better 
% plot(abs(st)) 
% xlabel('Time [s]') 
% ylabel('Amplitude (of a complex signal)') 
% title('Signal s(t) in time domain') 
NFFT = 2^14;                            %FFT size 
f = -Fs/2:1/(NFFT*Ts):Fs/2-1/(NFFT*Ts); %frequency vector 

% Plot the transmit signal in frequency domain 
% figure(4) 
% subplot(2,2,1); 
% plot(f/1e6, fftshift(abs(fft(st, NFFT)))); 
% xlabel('Frequency [MHz]') 
% ylabel('Amplitude ') 
% title('TX signal s(t)') 
% ylim([0 500]); 

%Generate the noise vector 
% Complex white Gaussian random noise 
n = (1/sqrt(2))*(randn(size(st)) + 1j*randn(size(st))); 
P_s = var(st);              % Signal power 
P_n = var(n);               % Noise power  


n1 = (1/sqrt(2))*(randn(size(st1)) + 1j*randn(size(st1))); 
P_s1 = var(st1);              % Signal power 
P_n1 = var(n1);               % Noise power  

% Defining noise scaling factor based on the desired SNR: 
noise_scaling_factor = sqrt(P_s/P_n./10.^(SNR./10)*(r/(1+alpha))); 
noise_scaling_factor1 = sqrt(P_s1/P_n1./10.^(SNR./10)*(r/(1+alpha))); 
%Add noise on top of the signal s(t). Remember that the variable “SNR” is now a vector.  
% Initialization for RX signal matrix, where each row represents the 
% received signal with a specific SNR value 
rt = zeros(length(SNR), length(st)); 
rt1 = zeros(length(SNR), length(st1)); 


% Received signal with different SNR values 
for ii = 1:1:length(SNR)     
rt(ii,:) = st + noise_scaling_factor(ii)*n; 
rt1(ii,:) = st1 + noise_scaling_factor1(ii)*n1; 
end 

% Plot the amplitude response of the noise with the SNR corresponding  
% to the last value in the SNR vector (just as an example) 
% figure(4) 
% subplot(2,2,2) 
% plot(f/1e6, fftshift(abs(fft(noise_scaling_factor(end)*n, NFFT)))); 
% xlabel('Frequency [MHz]') 
% ylabel('Amplitude') 
% title(['Noise (corresponding to SNR = ', num2str(SNR(end)), ' dB)']) 
% ylim([0 500]); 

% Received signal with the noise when the SNR is equal to the last value of 
% the SNR vector 
% figure(4) 
% subplot(2,2,3) 
% plot(f/1e6, fftshift(abs(fft(rt(end,:), NFFT)))); 
% xlabel('Frequency [MHz]') 
% ylabel('Amplitude') 
% title(['RX signal r(t) (SNR = ', num2str(SNR(end)), ' dB)']) 
% ylim([0 500]); 

%RECEIVER STRUCTURE 
% Creating the receive filter (it is the same as in the transmitter) 
ft = gt;   % Plotting the amplitude response of the receive filter figure(4) 
ft1 = gt1;
% subplot(2,2,4) 
% plot(f/1e6, fftshift(abs(fft(ft, NFFT)))); 
% xlabel('Frequency [MHz]') 
% ylabel('Amplitude') 
% title('RX filter f(t)') 

% Initialization for the received symbol matrix, where each row represents 
% the symbols with a specific SNR value 
qk = zeros(length(SNR), N_symbols - N_symbols_per_pulse); 
qk1 = zeros(length(SNR), N_symbols - N_symbols_per_pulse); 

% Filtering and sampling 
for ii = 1:1:length(SNR)     
    qt = filter(ft,1,rt(ii,:));  
    qt1 = filter(ft1,1,rt1(ii,:));  
    
    % Receiver filtering      
    qt = qt(1+(length(ft)-1)/2:end); 
    qt1 = qt1(1+(length(ft1)-1)/2:end);
    % Filter delay correction          
% Sampling the filtered signal. Remember that we used oversampling in     % the TX.     
qk(ii,:) = qt(1:r:end); 
qk1(ii,:) = qt1(1:r:end);
end

% Plot a few examples of the noisy samples and compare them with the 
% original symbol alphabet 
% figure(5) 
% subplot(3,1,1) 
% plot(qk(1,:),'b*') 
% hold on 
% plot(alphabet,'ro', 'MarkerFaceColor','r') 
% hold off 
% legend('Received samples', 'Original symbols') 
% 
% xlabel('Re') 
% ylabel('Im') 
% title(['Received samples with SNR = ', num2str(SNR(1)), ' dB']) 
% axis equal 

% figure(5) 
% subplot(3,1,2) 
mid_ind = ceil(length(SNR)/2);  % index for the entry in the middle 
% plot(qk(mid_ind,:),'b*') 
% hold on 
% plot(alphabet,'ro', 'MarkerFaceColor','r') 
% hold off 
% legend('Received samples ', 'Original symbols') 
% xlabel('Re') 
% ylabel('Im') 
% title(['Received samples with SNR = ', num2str(SNR(mid_ind)), ' dB']) 
% axis equal 

% Initialization 
BER1 = zeros(1,length(SNR));   
for ii = 1:1:length(SNR)    
     alphabet_error_matrix = abs(bsxfun(@minus,alphabet.',qk(ii,:))); 
     alphabet_error_matrix1 = abs(bsxfun(@minus,alphabet.',qk1(ii,:))); 
     % Searching for the indeces corresponding to the minimum distances    
     [~,estimated_Gray_symbol_ind] = min(alphabet_error_matrix); 
     [~,estimated_Gray_symbol_ind1] = min(alphabet_error_matrix1); 
     
     
      estimated_symbol_indices = ...         
          gray2bin(estimated_Gray_symbol_ind-1,'qam',M);  
      estimated_symbol_indices1 = ...         
          gray2bin(estimated_Gray_symbol_ind1-1,'qam',M);      
          % block of estimated bits in the columns     
          estimated_bit_blocks = ...         
              rem(floor((estimated_symbol_indices(:))*2.^(1-log2(M):0)),2)'; 
          estimated_bit_blocks1 = ...         
              rem(floor((estimated_symbol_indices1(:))*2.^(1-log2(M):0)),2)'; 
    % Bit blocks to bit vector     
    estimated_bits = estimated_bit_blocks(:);  
    estimated_bits1 = estimated_bit_blocks1(:); 
    % Finding out which bits were estimated incorrecly:     
    bit_errors = ...         
        estimated_bits ~= bits(1:length(estimated_bits));   
    bit_errors1 = ...         
        estimated_bits1 ~= bits(1:length(estimated_bits1));  
    % Bit error rate (0 means 0% of errors, 1 means 100% of errors)     
    BER1(1,ii) = mean(bit_errors);   
    BER2(1,ii) = mean(bit_errors1); 
end 

 



%----- Compute Theoretical Symbol Error Rates --------------------
 %EsN0lin = 10.^(EsN0dB/10); 
 EbN0linRC = 10.^(SNR/10); 
 symErrTheoryRC = 2*(1-1/sqrt(M))*erfc(sqrt(3/2*k*EbN0linRC/(M-1)));
 
plotHandle=plot(SNR,BER,char(colors(ind))); 
set(plotHandle,'LineWidth',1.5); 
ind=ind+1;


plotHandle=plot(SNR,BER1,char(colors(ind))); 
set(plotHandle,'LineWidth',1.5); 
 
 
 
 legend('RectQPSK','RCQPSK','RRCQPSK')
 axis([-20 5 0 1]); 
 set(gca,'XTick',-20:2:5); 
 xlabel('SNR'); 
 ylabel('Probability of BER Error');
 title('Probability of BER Error log10(Pb) Vs SNR'); 
 grid on; 
