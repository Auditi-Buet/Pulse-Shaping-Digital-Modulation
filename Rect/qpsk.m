clc
clear all
bits=10000;
data=randi(1,bits)>0.5;
%---debugging---
%data=[1 1 1]
%xxxxxxxxxx
ebno=0:10;
BER=zeros(1,length(ebno));
  
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
for i=1:length(ebno)
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
    npsd=10.^(ebno(i)/10);
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
figure


semilogy(ebno,BER,'b.-');
hold on

thr=0.5*erfc(sqrt(10.^(ebno/10)));
semilogy(ebno,thr,'rx-');
xlabel('Eb/No (dB)')
ylabel('Bit Error rate')
title('Simulated Vs Theoritical Bit Error Rate for QPSK')
legend('Simulation','Theory')
grid on



