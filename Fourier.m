% System Identification
% Fast Fourier transform (FFT)
% 
% Author: Jungwei Huang, NSYSU

clear all; close all; clc;

SIGNAL_START = 0; %(sec)
SIGNAL_TIME = 30; %(sec)
SAMPLE_TIME = 59; %(sec) 
SAMPLE_RATE = 2048; % (hz)

CUTOFF_FREQUENCY = 64; %(hz)
MODEL_FREQUENCY = 256; %(hz)

DATA_DIR = "C:/Users/Public/Downloads/Fitplus-7355/Impulse/";
MODEL_DIR = "C:/Users/Public/Downloads/Fitplus-7355/Model/";

for n = 1
  
    %%%% Load data %%%%
    raw = load([DATA_DIR, num2str(n), ".txt"]);
    output = raw(1+SAMPLE_RATE*SIGNAL_START:SAMPLE_RATE*(SIGNAL_START+SAMPLE_TIME),1)+ raw(1+SAMPLE_RATE*SIGNAL_START:SAMPLE_RATE*(SIGNAL_START+SAMPLE_TIME),2) + raw(1+SAMPLE_RATE*SIGNAL_START:SAMPLE_RATE*(SIGNAL_START+SAMPLE_TIME),3)+ raw(1+SAMPLE_RATE*SIGNAL_START:SAMPLE_RATE*(SIGNAL_START+SAMPLE_TIME),4);
    input = raw(1+SAMPLE_RATE*SIGNAL_START:SAMPLE_RATE*(SIGNAL_START+SAMPLE_TIME),5) + raw(1+SAMPLE_RATE*SIGNAL_START:SAMPLE_RATE*(SIGNAL_START+SAMPLE_TIME),6) + raw(1+SAMPLE_RATE*SIGNAL_START:SAMPLE_RATE*(SIGNAL_START+SAMPLE_TIME),7) + raw(1+SAMPLE_RATE*SIGNAL_START:SAMPLE_RATE*(SIGNAL_START+SAMPLE_TIME),8);   

    %%%% Modeling with FFT %%%%
    input = input - mean(input);
    input_fft = fft(input(1+15*SAMPLE_RATE:45*SAMPLE_RATE),SIGNAL_TIME*SAMPLE_RATE);
    [input_theta input_rho]=cart2pol(real(input_fft),imag(input_fft));
    input_bandwidth = bandwidth(input, SAMPLE_RATE, 95);
    
    output = output - mean(output);
    output_fft = fft(output(1+15*SAMPLE_RATE:45*SAMPLE_RATE),SIGNAL_TIME*SAMPLE_RATE);    
    [output_theta output_rho]=cart2pol(real(output_fft),imag(output_fft));
    output_bandwidth = bandwidth(output, SAMPLE_RATE, 95);
        
    H_fft = output_fft./input_fft; % transfer function
    [H_theta H_rho] = cart2pol(real(H_fft),imag(H_fft));  
    H_Dnmag = 20*log10(H_rho);  
    
    %%%% Spectral analysis %%%%
    fig11 =figure(11) %time domain
    %set( fig101 , 'Visible' , 'off' )
    subplot(2,1,1); plot((1:length(input(1:SIGNAL_TIME*(SAMPLE_RATE))))/SAMPLE_RATE, input(1:SIGNAL_TIME*(SAMPLE_RATE)));
    title('Input Signal'); xlabel('time(sec)'); ylabel('amplitude');      
    subplot(2,1,2); plot((1:length(output(1:SIGNAL_TIME*(SAMPLE_RATE))))/SAMPLE_RATE, output(1:SIGNAL_TIME*(SAMPLE_RATE)));
    title('Output Signal'); xlabel('time(sec)'); ylabel('amplitude');    
    %saveas( fig11 , [MODEL_DIR, num2str(n), '_signal.png' ])
     
    fig12 = figure(12); %frequency domain
    %set( fig12 , 'Visible' , 'off' )
    subplot(2,2,1)
    plot((1:length(input_rho)/2)/length(input_rho)*SAMPLE_RATE, input_rho(2:length(input_rho)/2+1))
    plot((1:length(input_rho)/2)/length(input_rho)*SAMPLE_RATE, input_rho(2:length(input_rho)/2+1))
    axis([0, SAMPLE_RATE/2,-inf,inf]); xlabel('Hz'); ylabel('amplitude'); title('Input Spectrum (Amplitude)','FontSize',14);
    subplot(2,2,2)
    plot((1:length(input_theta)/2)/length(input_theta)*SAMPLE_RATE,(input_theta(2:length(input_theta)/2+1))/pi)
    axis([0, SAMPLE_RATE/2,-inf,inf]); xlabel('Hz'); ylabel('phase(rad)'); title('Input Spectrum (Phase)','FontSize',14);     
    subplot(2,2,3)
    plot((1:length(output_rho)/2)/length(output_rho)*SAMPLE_RATE, output_rho(2:length(output_rho)/2+1))
    axis([0,SAMPLE_RATE/2,-inf,inf]); xlabel('Hz'); ylabel('amplitude'); title('Input Spectrum (Amplitude)','FontSize',14);
    subplot(2,2,4)
    plot((1:length(output_theta)/2)/length(output_theta)*SAMPLE_RATE, (output_theta(2:length(output_theta)/2+1))/pi)
    axis([0, SAMPLE_RATE/2,-inf,inf]); xlabel('Hz'); ylabel('phase(rad)'); title('Output Spectrum (Amplitude)','FontSize',14)
    %saveas( fig12 , [MODEL_DIR, num2str(n), '_spectrum.png' ])

    fig13 = figure(13) %transfer function
    %set( fig13 , 'Visible' , 'off' )
    subplot(2,1,1)
    plot((1:length(H_Dnmag)/2)/length(H_Dnmag)*SAMPLE_RATE, H_Dnmag(2:length(H_Dnmag)/2+1)); 
    axis([0, SAMPLE_RATE/2, -inf, inf]); xlabel('Hz'); ylabel('dB'); title('Transfer Function (Amplitude)');
    subplot(2,1,2)
    plot((1:length(H_theta)/2)/length(H_theta)*SAMPLE_RATE, (H_theta(2:length(H_theta)/2+1))/pi)
    axis([0, SAMPLE_RATE/2, -inf, inf]); xlabel('Hz'); ylabel('rad');
    title('Transfer Function (Phase)')
    %saveas( fig13 , [MODEL_DIR, num2str(n), '_transfer-function.png' ])
 
 
 


     
    input_pred = input_pred(1+12*SAMPLE_RATE:48*SAMPLE_RATE);
    input_actual = (input(1+12*SAMPLE_RATE:48*SAMPLE_RATE))';

    [input_pred, cutoff_freq, si] = butterworthIIR(input_pred, 25, 30, SAMPLE_RATE);
    [input_actual, cutoff_freq, ci] = butterworthIIR(input_actual, 25, 30, SAMPLE_RATE);    
 
    input_pred = input_pred(1+3*SAMPLE_RATE:33*SAMPLE_RATE); % truncate the signal between 15sec~45sec
    input_actual = input_actual(1+3*SAMPLE_RATE:33*SAMPLE_RATE);
     
    input_pred_bandwidth = bandwidth(input_pred, SAMPLE_RATE, 95);
    input_actual_bandwidth = bandwidth(input_actual, SAMPLE_RATE, 95);
     
    input_pred = (input_pred - mean(input_pred))/std(input_pred);    
    input_actual = (input_actual - mean(input_actual))/std(input_actual);
     
    harmonic_distor = sum((input_actual - input_pred).^2)./sum(input_actual.^2);
    signal_corr = corrcoef(input_actual, input_pred);
     
  
    figure(21)
    subplot(2,1,1)
    plot((1:length(input_pred_rho)/2)/length(input_pred_rho)*SAMPLE_RATE, input_pred_rho(2:length(input_pred_rho)/2+1))
    axis([0,SAMPLE_RATE/2,-inf,inf]); xlabel('Hz'); ylabel('amplitude'); title('Predict Input Spectrum (Amplitude)');
    subplot(2,1,2)
    plot((1:length(input_pred_theta)/2)/length(input_pred_theta)*SAMPLE_RATE,(input_pred_theta(2:length(input_pred_theta)/2+1))/pi)
    axis([0,SAMPLE_RATE/2,-inf,inf]); xlabel('Hz'); ylabel('phase(rad)'); title('Predict Input Spectrum (Phase)');
 
    figure(22)
    subplot(2,1,1)
    plot((1:length(input_actual_rho)/2)/length(input_actual_rho)*SAMPLE_RATE, input_actual_rho(2:length(input_actual_rho)/2+1))
    axis([0,SAMPLE_RATE/2,-inf,inf]); xlabel('Hz'); ylabel('amplitude'); title('Actual Input Spectrum (Amplitude)');
    subplot(2,1,2)
    plot((1:length(input_actual_theta)/2)/length(input_actual_theta)*SAMPLE_RATE,(input_actual_theta(2:length(input_actual_theta)/2+1))/pi)
    axis([0,SAMPLE_RATE/2,-inf,inf]); xlabel('Hz'); ylabel('phase(rad)'); title('Actual Input Spectrum (Phase)')
     
end