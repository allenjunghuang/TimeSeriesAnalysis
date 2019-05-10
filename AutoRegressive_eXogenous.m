% System Identification
% AutoRegressive eXogenous (ARX)
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
    output = raw(1+SAMPLE_RATE*SIGNAL_START:SAMPLE_RATE*(SIGNAL_START+SAMPLE_TIME),1) + raw(1+SAMPLE_RATE*SIGNAL_START:SAMPLE_RATE*(SIGNAL_START+SAMPLE_TIME),2) + raw(1+SAMPLE_RATE*SIGNAL_START:SAMPLE_RATE*(SIGNAL_START+SAMPLE_TIME),3) + raw(1+SAMPLE_RATE*SIGNAL_START:SAMPLE_RATE*(SIGNAL_START+SAMPLE_TIME),4);
    input = raw(1+SAMPLE_RATE*SIGNAL_START:SAMPLE_RATE*(SIGNAL_START+SAMPLE_TIME),5) + raw(1+SAMPLE_RATE*SIGNAL_START:SAMPLE_RATE*(SIGNAL_START+SAMPLE_TIME),6) + raw(1+SAMPLE_RATE*SIGNAL_START:SAMPLE_RATE*(SIGNAL_START+SAMPLE_TIME),7) + raw(1+SAMPLE_RATE*SIGNAL_START:SAMPLE_RATE*(SIGNAL_START+SAMPLE_TIME),8);   

    %%%%% Preprocessing %%%%%
    [input, cutoff_freq, ai] = buterworthIIR(input, CUTOFF_FREQUENCY, CUTOFF_FREQUENCY + 5, SAMPLE_RATE); %butterworth filter     
    [output, cutoff_freq, ao] = buterworthIIR(output, CUTOFF_FREQUENCY, CUTOFF_FREQUENCY + 5, SAMPLE_RATE);   
    
    model_input = resample(output, MODEL_FREQUENCY, SAMPLE_RATE);
    model_output = resample(input, MODEL_FREQUENCY, SAMPLE_RATE); 
    
    %%%%% ARX inverse %%%%%
    input_initial = mean(model_input(1+3*MODEL_FREQUENCY:8*MODEL_FREQUENCY)); #remove dc component
    output_initial = mean(model_output(1+3*MODEL_FREQUENCY:8*MODEL_FREQUENCY));
    model_input = (model_input - input_initial);
    model_output = (model_output - output_initial);  
    
    Data = iddata(model_output', model_input', 1/MODEL_FREQUENCY);
    arx_est = struc(1:32,1:32,1:16);
    
    %%%% Modeling with estimation and validation data %%%%
     ze = detrend(Data(1+9*MODEL_FREQUENCY:29*MODEL_FREQUENCY)); %estimation data
     zr = detrend(Data(1+29*MODEL_FREQUENCY:49*MODEL_FREQUENCY)); %validation data
     arx_loss = arxstruc(ze,zr,arx_est);
     arx_order = selstruc(arx_loss,0);
     arx_model = arx(ze,arx_order,'Focus','Stability');
    
    %%%% Modeling ONLY with estimation data %%%%
    ze2 = Data(1+15*MODEL_FREQUENCY:45*MODEL_FREQUENCY); %estimation data
    arx_loss2 = arxstruc(ze2,ze2,arx_est);
    arx_order2 = selstruc(arx_loss2,'AIC');
    arx_model2 = arx(ze2,arx_order2,'Focus','Stability');
    
    %%%% Modeling ONLY with estimation data, select struct with FPE %%%%
    fpe_min = 10000000;
    arx_est = struc(16:32,16:32,1:8); 
    for kk = 1:size(arx_est, 1)
        models = arx(ze2,arx_est(kk,:),'Focus','Stability');
        arx_fpe(kk)=fpe(models);
        if arx_fpe(kk) < fpe_min
            arx_model3 = models;
            fpe_min = arx_fpe(kk);
        end   
    end
    
    
    %%%% Estimate model %%%%
    input_initial = mean(input(1+5*SAMPLE_RATE:9*SAMPLE_RATE));
    output_initial = mean(output(1+5*SAMPLE_RATE:9*SAMPLE_RATE));
    input_normalized = (input - input_initial)/std(input);
    output_normalized = (output - output_initial)/std(output);
    Data=iddata(output_normalized, input_normalized, 1/SAMPLE_RATE);
%    ze = detrend(Data(1+9*SAMPLE_RATE:29*SAMPLE_RATE)); % estimation data
%    zr = detrend(Data(1+29*SAMPLE_RATE:49*SAMPLE_RATE)); % validation data
    ze = Data(1+19*SAMPLE_RATE:39*SAMPLE_RATE); % estimation data
    zr = Data(1+29*SAMPLE_RATE:49*SAMPLE_RATE); % validation data
    
    arx_est = struc(1:12,1:12,0:12);  
    arx_loss = arxstruc(ze,ze,arx_est);
    arx_order = selstruc(arx_loss,0);
    arx_model = arx(ze, arx_order); 
    arx_err = resid(ze, arx_model);    
    an_len = length(arx_model.a);
    bn_len = length(arx_model.b);

    input_pred = zeros(1,8*SAMPLE_RATE);
    for p=(1+8*SAMPLE_RATE):48*SAMPLE_RATE
        an = 0;
        bn = 0;
        for k = 2:an_len
            an = an + input_pred(p-k+1)*arx_model.a(k);
        end        
        for s = 1:bn_len
            bn = bn + output_normalized(p-s+1)*arx_model.b(s);
        end
        input_pred(p) = bn - an;
        q_A(p) = an;
        q_B(p) = bn;
    end
     
%     input_pred = zeros(1,8*SAMPLE_RATE);
%     for p=(1+8*SAMPLE_RATE):48*SAMPLE_RATE
%         an = 0;
%         bn = 0;
%         for k = 1:an_len
%             an = an + output_normalized(p-k+1)*arx_model.a(k);
%         end        
%         for s = 2:bn_len
%             bn = bn + input_pred2(p-s+1)*arx_model.b(s);
%         end
%         input_pred(p)=(an - bn)/arx_model.b(1);
%         q_A(p) = an;
%         q_B(p) = bn;
%     end

    input_pred = input_pred(1+15*SAMPLE_RATE:35*SAMPLE_RATE);
    input_actual = input_normalized(1+15*SAMPLE_RATE:35*SAMPLE_RATE);
     
    input_pred = (input_pred - mean(input_pred))/std(input_pred);    
    input_actual = (input_actual - mean(input_actual))/std(input_actual);
    
    harmonic_distor = sum((input_actual - input_pred).^2)./sum(input_actual.^2);
    signal_corr = corrcoef(input_actual, input_pred);
     
         
    figure(31)
    subplot(2,1,1); plot(q_A); title('A(q) Recursive');
    subplot(2,1,2); plot(q_B); title('B(q) Recursive');
     
    figure(32)
    plot((15*SAMPLE_RATE+1:45*SAMPLE_RATE)/SAMPLE_RATE, input_actual, 'r'); 
    hold on;
    plot((15*SAMPLE_RATE+1:45*SAMPLE_RATE)/SAMPLE_RATE, input_pred, 'b'); 
    axis([30,33,-2,4]);xlabel('time(sec)'); ylabel('z-score'); title('Simulation');
    legend('Actual input signal', 'Predict input signal');
     
    figure(41)
    compare(zr, arx_model)
     
    figure(42)
    ffplot(arx_model);    
     
    figure(43)
    bode(arx_model);
     
    figure(44)
    plot(arx_err(:,1,[])), title('Residuals') 


    %%%% Save model %%%%
    save([MODEL_DIR, num2str(n), "_ARX.mat"], 'arx_model2','arx_model3','arx_est','arx_fpe');
    
end



