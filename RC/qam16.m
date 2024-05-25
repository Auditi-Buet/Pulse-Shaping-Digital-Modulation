SNR = 0:2:24;               % Signal-to-noise ratio vector [dB] 
T = 1/10e6;                 % Symbol time interval [s] 
Fd=(1/T);
r = 6;                      % Oversampling factor (r samples per pulse) 
N_symbols_per_pulse = 60;   % Duration of TX/RX filter in symbols 
alpha = 0.50;               % Roll-off factor (excess bandwidth) 

%Based on above we can define sampling frequency and sampling time interval as 
Fs = r/T;               % Sampling frequency 
Ts = 1/Fs;              % Sampling time interval 

%GENERATION OF BIT SEQUENCE AND QAM SYMBOLS 
N_symbols = 1000000;      % Number of symbols 

% Alphabet size 
M = 16;  % Number of symbols in the QAM alphabet (e.g. 16 means 16-QAM).

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

[Gray_symbol_indices, ~] = bin2gray(symbol_indices, 'qam', M); 
symbols = alphabet(Gray_symbol_indices+1); 

%TRANSMITTER STRUCTURE 
%Implement the transit filter: Root-Raised-Cosine (RRC) and plot the pulse shape 
% Filter generation 
gt = rcosdesign(alpha,N_symbols_per_pulse,r,'normal'); 

% Plot the pulse shape of the transmit/receive filter 
figure 
plot(-N_symbols_per_pulse*r/2*Ts:Ts:N_symbols_per_pulse*r/2*Ts,gt,'b') 
hold on 
stem(-N_symbols_per_pulse*r/2*Ts:T:N_symbols_per_pulse*r/2*Ts,gt(1:r:end),'ro') 
xlabel('time [s]') 
ylabel('Amplitude') 
title('Transmit/receive RC filter (pulse shape)') 
legend('Pulse shape','Ideal symbol-sampling locations') 

%up-sampled symbol sequence 
symbols_upsampled = zeros(size(1:r*N_symbols));      
% Zero vector initilized for  
% symbol insertion 
symbols_upsampled(1:r: r*N_symbols) = symbols; 
% now the up-sampled sequence looks like {a1 0 0... a2 0 0... a3 0 0...}   
st = filter(gt,1,symbols_upsampled); % Transmitter filtering 
st = st(1+(length(gt)-1)/2:end);     % Filter delay correction 

%Plot the transmit signal s(t) in time and frequency domain 
figure % zoom manually to see the signal better 
plot(abs(st)) 
xlabel('Time [s]') 
ylabel('Amplitude (of a complex signal)') 
title('Signal s(t) in time domain') 
NFFT = 2^14;                            %FFT size 
f = -Fs/2:1/(NFFT*Ts):Fs/2-1/(NFFT*Ts); %frequency vector 

% Plot the transmit signal in frequency domain 
figure(4) 
subplot(2,2,1); 
plot(f/1e6, fftshift(abs(fft(st, NFFT)))); 
xlabel('Frequency [MHz]') 
ylabel('Amplitude ') 
title('TX signal s(t)') 
ylim([0 500]); 

%Generate the noise vector 
% Complex white Gaussian random noise 
n = (1/sqrt(2))*(randn(size(st)) + 1j*randn(size(st))); 
P_s = var(st);              % Signal power 
P_n = var(n);               % Noise power  

% Defining noise scaling factor based on the desired SNR: 
noise_scaling_factor = sqrt(P_s/P_n./10.^(SNR./10)*(r/(1+alpha))); 
%Add noise on top of the signal s(t). Remember that the variable “SNR” is now a vector.  
% Initialization for RX signal matrix, where each row represents the 
% received signal with a specific SNR value 
rt = zeros(length(SNR), length(st)); 

% Received signal with different SNR values 
for ii = 1:1:length(SNR)     
rt(ii,:) = st + noise_scaling_factor(ii)*n; 
end 

% Plot the amplitude response of the noise with the SNR corresponding  
% to the last value in the SNR vector (just as an example) 
figure(4) 
subplot(2,2,2) 
plot(f/1e6, fftshift(abs(fft(noise_scaling_factor(end)*n, NFFT)))); 
xlabel('Frequency [MHz]') 
ylabel('Amplitude') 
title(['Noise (corresponding to SNR = ', num2str(SNR(end)), ' dB)']) 
ylim([0 500]); 

% Received signal with the noise when the SNR is equal to the last value of 
% the SNR vector 
figure(4) 
subplot(2,2,3) 
plot(f/1e6, fftshift(abs(fft(rt(end,:), NFFT)))); 
xlabel('Frequency [MHz]') 
ylabel('Amplitude') 
title(['RX signal r(t) (SNR = ', num2str(SNR(end)), ' dB)']) 
ylim([0 500]); 

%RECEIVER STRUCTURE 
% Creating the receive filter (it is the same as in the transmitter) 
ft = gt;   % Plotting the amplitude response of the receive filter figure(4) 
subplot(2,2,4) 
plot(f/1e6, fftshift(abs(fft(ft, NFFT)))); 
xlabel('Frequency [MHz]') 
ylabel('Amplitude') 
title('RX filter f(t)') 

% Initialization for the received symbol matrix, where each row represents 
% the symbols with a specific SNR value 
qk = zeros(length(SNR), N_symbols - N_symbols_per_pulse); 

% Filtering and sampling 
for ii = 1:1:length(SNR)     
    qt = filter(ft,1,rt(ii,:));      
    % Receiver filtering      
    qt = qt(1+(length(ft)-1)/2:end); 
    % Filter delay correction          
% Sampling the filtered signal. Remember that we used oversampling in     % the TX.     
qk(ii,:) = qt(1:r:end); 
end

% Plot a few examples of the noisy samples and compare them with the 
% original symbol alphabet 
figure(5) 
subplot(3,1,1) 
plot(qk(1,:),'b*') 
hold on 
plot(alphabet,'ro', 'MarkerFaceColor','r') 
hold off 
legend('Received samples', 'Original symbols') 

xlabel('Re') 
ylabel('Im') 
title(['Received samples with SNR = ', num2str(SNR(1)), ' dB']) 
axis equal 

figure(5) 
subplot(3,1,2) 
mid_ind = ceil(length(SNR)/2);  % index for the entry in the middle 
plot(qk(mid_ind,:),'b*') 
hold on 
plot(alphabet,'ro', 'MarkerFaceColor','r') 
hold off 
legend('Received samples ', 'Original symbols') 
xlabel('Re') 
ylabel('Im') 
title(['Received samples with SNR = ', num2str(SNR(mid_ind)), ' dB']) 
axis equal 

% Initialization 
BER = zeros(1,length(SNR));   
for ii = 1:1:length(SNR)    
     alphabet_error_matrix = abs(bsxfun(@minus,alphabet.',qk(ii,:))); 
     % Searching for the indeces corresponding to the minimum distances    
     [~,estimated_Gray_symbol_ind] = min(alphabet_error_matrix); 
     
      estimated_symbol_indices = ...         
          gray2bin(estimated_Gray_symbol_ind-1,'qam',M);          
          % block of estimated bits in the columns     
          estimated_bit_blocks = ...         
              rem(floor((estimated_symbol_indices(:))*2.^(1-log2(M):0)),2)'; 
    % Bit blocks to bit vector     
    estimated_bits = estimated_bit_blocks(:);          
    % Finding out which bits were estimated incorrecly:     
    bit_errors = ...         
        estimated_bits ~= bits(1:length(estimated_bits));              
    % Bit error rate (0 means 0% of errors, 1 means 100% of errors)     
    BER(1,ii) = mean(bit_errors);   
end 

% Find out the minimum distance between two constellation points 
d = min(abs(alphabet(1)-alphabet(2:end)));
sigma = sqrt(0.5 * P_n * noise_scaling_factor.^2); 
% Theoretical symbol error probability applicable to all M-QAM alphabets. 
P_sym_error = (4*qfunc(d./(2*sigma)) - 4*qfunc(d./(2*sigma)).^2) * ...     
    ((M-4-4*(sqrt(M)-2))/M) + ...     
    (2*qfunc(d./(2*sigma))- qfunc(d./(2*sigma)).^2) * (4/M) + ...     
    (3*qfunc(d./(2*sigma)) - 2*qfunc(d./(2*sigma)).^2) * (4*(sqrt(M)-2)/M); 
% Compare the simulated and theoretical results. 
figure 
semilogy(SNR, BER,'r-'); 
hold on; 
semilogy(SNR, (1/log2(M))*P_sym_error,'b*'); 
title('Bit error rate') 
xlabel('SNR [dB]') 
ylabel('BER') 
legend('Simulated BER with Gray coding',...     
    'Theoretical bit error probability (approx.)','Location', 'SouthWest');