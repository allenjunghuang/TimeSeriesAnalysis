% System Identification
% Estimate models (Cross-validation)
% 
% Author: Jungwei Huang, NSYSU

clear all; close all; clc;

for kk=1
    ss=num2str(kk);
    arx_model=load(['H:\Lab\lab_data\calibration_data\step_modeling\64hz_arx\' ss '.mat']);   
    
    an_leng=length(arx_model.arx_model.a);
    bn_leng=length(arx_model.arx_model.b);
    %%%step response
    output_S=ones(1,256*5);
    simu_S=zeros(1,1*256);
    for p=(1+1*256):5*256
        temp_a=0;
        temp_b=0;
        for k=2:an_leng
            temp_a=temp_a+simu_S(p-k+1)*arx_model.arx_model.a(k);
        end
        for s=1:bn_leng
            temp_b=temp_b+output_S(p-s+1)*arx_model.arx_model.b(s);
        end
        simu_S(p)=temp_b-temp_a;
        qx_A(p)=temp_a;
        qx_B(p)=temp_b;
    end    
    temp_S=simu_S-simu_S(end)*0.98; %check settling time
    ts_indx=find(temp_S(1:end-1).*temp_S(2:end)<0);
    ts_samp=max(ts_indx)-(1+1*256);
    ts_time=ts_samp/256;
    
    figure(1)
    plot(simu_S); 
    axis([256,512,-inf,inf]);xlabel('index'); ylabel('amplitude'); title('step response');
    
%     for pp=1
%         clf;
%         rr=num2str(pp);
%         %verify data
%         signal_start=30;
%         sample_rate=512;
%         sample_time=59;
% %         data_A=load(['H:\Lab\lab_data\calibration_data\step_modeling\' rr '.txt']);
% %         output_A=data_A(1+sample_rate*signal_start:sample_rate*(signal_start+sample_time),1)+data_A(1+sample_rate*signal_start:sample_rate*(signal_start+sample_time),2)+data_A(1+sample_rate*signal_start:sample_rate*(signal_start+sample_time),3)+data_A(1+sample_rate*signal_start:sample_rate*(signal_start+sample_time),4);
% %         input_A=data_A(1+sample_rate*signal_start:sample_rate*(signal_start+sample_time),5)+data_A(1+sample_rate*signal_start:sample_rate*(signal_start+sample_time),6)+data_A(1+sample_rate*signal_start:sample_rate*(signal_start+sample_time),7)+data_A(1+sample_rate*signal_start:sample_rate*(signal_start+sample_time),8);   
%         data_A=load(['H:\Lab\lab_data\calibration_data\stepping_1310\gary_simurun_1001.txt']);
%         input_A=data_A(1+sample_rate*signal_start:sample_rate*(signal_start+sample_time),1)+data_A(1+sample_rate*signal_start:sample_rate*(signal_start+sample_time),2)+data_A(1+sample_rate*signal_start:sample_rate*(signal_start+sample_time),3)+data_A(1+sample_rate*signal_start:sample_rate*(signal_start+sample_time),4);
%         output_A=data_A(1+sample_rate*signal_start:sample_rate*(signal_start+sample_time),5)+data_A(1+sample_rate*signal_start:sample_rate*(signal_start+sample_time),6)+data_A(1+sample_rate*signal_start:sample_rate*(signal_start+sample_time),7)+data_A(1+sample_rate*signal_start:sample_rate*(signal_start+sample_time),8);
%         input_A=resample(input_A, 256,sample_rate);
%         output_A=resample(output_A, 256,sample_rate);
%         
%         sample_rate=256;
%         signal_time=30;
%         intial_input_A=mean(input_A(1+3*sample_rate:8*sample_rate));
%         intial_output_A=mean(output_A(1+3*sample_rate:8*sample_rate));
%         input_A=(input_A-intial_input_A);
%         output_A=(output_A-intial_output_A);
%         
%         an_leng=length(arx_model.arx_model.a);
%         bn_leng=length(arx_model.arx_model.b);
%         simu_A=zeros(1,8*sample_rate);
%         for p=(1+8*sample_rate):48*sample_rate
%             temp_a=0;
%             temp_b=0;
%             for k=2:an_leng
%                 temp_a=temp_a+simu_A(p-k+1)*arx_model.arx_model.a(k);
%             end
%             for s=1:bn_leng
%                 temp_b=temp_b+output_A(p-s+1)*arx_model.arx_model.b(s);
%             end
%             simu_A(p)=temp_b-temp_a;
%         end
%         
%         simu_A=simu_A(1+12*sample_rate:48*sample_rate);
%         real_A=(input_A(1+12*sample_rate:48*sample_rate))';
%         
%         simu_A_fft=fft(simu_A(1+3*sample_rate:33*sample_rate),signal_time*sample_rate);
%         [simu_A_theta simu_A_rho]=cart2pol(real(simu_A_fft),imag(simu_A_fft)); %transform Cartesian to polar
%         simu_A_rho_energy=sum(simu_A_rho(2:end/2+1).^2);
%         energy95_simu_A=0;
%         for k=1:length(simu_A_rho)
%             if energy95_simu_A<=0.95*simu_A_rho_energy
%                 energy95_simu_A=energy95_simu_A+simu_A_rho(1+k)^2;
%                 BW_simu_A=k/signal_time;
%             end
%         end
%         %     simu_A_Dnmag=20*log10(simu_A_rho);
%         
%         real_A_fft=fft(real_A(1+3*sample_rate:33*sample_rate),signal_time*sample_rate);
%         [real_A_theta real_A_rho]=cart2pol(real(real_A_fft),imag(real_A_fft)); %transform Cartesian to polar
%         real_A_rho_energy=sum(real_A_rho(2:end/2+1).^2);
%         energy95_real_A=0;
%         for k=1:length(real_A_rho)
%             if energy95_real_A<=0.95*real_A_rho_energy
%                 energy95_real_A=energy95_real_A+real_A_rho(1+k)^2;
%                 BW_real_A=k/signal_time;
%             end
%         end
%         %     real_A_Dnmag=20*log10(real_A_rho);
%         
%         corner_hz=25;
%         [simu_B,cutfreq_si,temp_si]=buter_filter3(simu_A,corner_hz,sample_rate);
%         [real_B,cutfreq_ci,temp_ci]=buter_filter3(real_A,corner_hz,sample_rate);
%         simu_C=simu_B(1+3*sample_rate:33*sample_rate); %原始60sec的訊號取15sec~45sec間的訊號做比較
%         real_C=real_B(1+3*sample_rate:33*sample_rate);
%         
%         simu_C_fft=fft(simu_C,signal_time*sample_rate);
%         [simu_C_theta simu_C_rho]=cart2pol(real(simu_C_fft),imag(simu_C_fft)); %transform Cartesian to polar
%         simu_C_rho_energy=sum(simu_C_rho(2:end/2+1).^2);
%         energy95_simu_C=0;
%         for k=1:length(simu_C_rho)
%             if energy95_simu_C<=0.95*simu_C_rho_energy
%                 energy95_simu_C=energy95_simu_C+simu_C_rho(1+k)^2;
%                 BW_simu_C=k/signal_time;
%             end
%         end
%         %     simu_C_Dnmag=20*log10(simu_C_rho);
%         
%         real_C_fft=fft(real_C,signal_time*sample_rate);
%         [real_C_theta real_C_rho]=cart2pol(real(real_C_fft),imag(real_C_fft)); %transform Cartesian to polar
%         real_C_rho_energy=sum(real_C_rho(2:end/2+1).^2);
%         energy95_real_C=0;
%         for k=1:length(real_C_rho)
%             if energy95_real_C<=0.95*real_C_rho_energy
%                 energy95_real_C=energy95_real_C+real_C_rho(1+k)^2;
%                 BW_real_C=k/signal_time;
%             end
%         end
%         %     real_C_Dnmag=20*log10(real_C_rho);
%         
%         simu_C=(simu_C-mean(simu_C))/std(simu_C);
%         real_C=(real_C-mean(real_C))/std(real_C);
%         real_C_energy=sum(real_C.^2);
%         harmon_dist_x=sum((real_C-simu_C).^2)./real_C_energy
%         signal_corr_x=corrcoef(real_C,simu_C)
%         
%         corr_data(pp,kk)=signal_corr_x(1,2);
%         distor_data(pp,kk)=harmon_dist_x;
%         
% %         figure(21)
% %         %     subplot(2,1,1)
% %         plot((1:length(simu_A_rho)/2)/length(simu_A_rho)*sample_rate,simu_A_rho(2:length(simu_A_rho)/2+1))
% %         axis([0,sample_rate/2,-inf,inf]); xlabel('Hz'); ylabel('amplitude'); title('濾波前模擬輸入頻譜(振幅)');
% %         %     subplot(2,1,2)
% %         %     plot((1:length(simu_A_theta)/2)/length(simu_A_theta)*sample_rate,(simu_A_theta(2:length(simu_A_theta)/2+1))/pi)
% %         %     axis([0,sample_rate/2,-inf,inf]); xlabel('Hz'); ylabel('phase(rad)'); title('濾波前模擬輸入頻譜(相位)');
% %         
% %         figure(22)
% %         %     subplot(2,1,1)
% %         plot((1:length(real_A_rho)/2)/length(real_A_rho)*sample_rate,real_A_rho(2:length(real_A_rho)/2+1))
% %         axis([0,sample_rate/2,-inf,inf]); xlabel('Hz'); ylabel('amplitude'); title('濾波前真實輸入頻譜(振幅)');
% %         %     subplot(2,1,2)
% %         %     plot((1:length(real_A_theta)/2)/length(real_A_theta)*sample_rate,(real_A_theta(2:length(real_A_theta)/2+1))/pi)
% %         %     axis([0,sample_rate/2,-inf,inf]); xlabel('Hz'); ylabel('phase(rad)'); title('濾波前真實輸入頻譜(相位)')
% %         
% %         figure(23)
% %         %     subplot(2,1,1)
% %         plot((1:length(simu_C_rho)/2)/length(simu_C_rho)*sample_rate,simu_C_rho(2:length(simu_C_rho)/2+1))
% %         axis([0,sample_rate/2,-inf,inf]); xlabel('Hz'); ylabel('amplitude'); title('濾波後模擬輸入頻譜(振幅)');
% %         %     subplot(2,1,2)
% %         %     plot((1:length(simu_C_theta)/2)/length(simu_C_theta)*sample_rate,(simu_C_theta(2:length(simu_C_theta)/2+1))/pi)
% %         %     axis([0,sample_rate/2,-inf,inf]); xlabel('Hz'); ylabel('phase(rad)'); title('濾波後模擬輸入頻譜(相位)');
% %         
% %         figure(24)
% %         %     subplot(2,1,1)
% %         plot((1:length(real_C_rho)/2)/length(real_C_rho)*sample_rate,real_C_rho(2:length(real_C_rho)/2+1))
% %         axis([0,sample_rate/2,-inf,inf]); xlabel('Hz'); ylabel('amplitude'); title('濾波後真實輸入頻譜(振幅)');
% %         %     subplot(2,1,2)
% %         %     plot((1:length(real_C_theta)/2)/length(real_C_theta)*sample_rate,(real_C_theta(2:length(real_C_theta)/2+1))/pi)
% %         %     axis([0,sample_rate/2,-inf,inf]); xlabel('Hz'); ylabel('phase(rad)'); title('濾波後真實輸入頻譜(相位)')
% %         
% %         figure(31)
% %         subplot(2,1,1); plot(qx_A); title('A(q) Recursive');
% %         subplot(2,1,2); plot(qx_B); title('B(q) Recursive');
%         
%         fig32=figure(32);
%         set( fig32 , 'Visible' , 'off' )
%         plot((15*sample_rate+1:45*sample_rate)/sample_rate,real_C,'r'); hold on;
%         plot((15*sample_rate+1:45*sample_rate)/sample_rate,simu_C,'b');
%         axis([36,37,-2.5,2.5]);xlabel('t(sec)'); ylabel('Z-score'); title('模擬結果');
%         legend('真實輸入訊號','模擬輸入訊號');
%         saveas( fig32 , ['H:\Lab\lab_data\calibration_data\fig' ss '' rr '.png' ])
%         
% %         figure(41)
% %         ffplot(arx_model);
% %         figure(42)
% %         bode(arx_model);
% %         figure(43)
% %         plot(arx_err(:,1,[])), title('The residuals')
%         
% %         clear data_A output_A input_A intial_input_A intial_output_A
% %         clear simu_A simu_B simu_C real_A real_B real_C
% %         clear simu_A_fft simu_C_fft real_A_fft real_C_fft
%     end   
% %     clear arx_model

end
% save(['H:\Lab\lab_data\calibration_data\part1.mat' ],'corr_data','distor_data');