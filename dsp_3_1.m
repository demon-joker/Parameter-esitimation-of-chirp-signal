clear all;

rng(44);

% 参数设置
fs = 1e7;          % 采样率为10MHz
time = 1e-3;       % 整段时间为1ms
t_all=0:1/fs:time-1/fs; 

% 生成观测信号

    t_pulse=3e-4;   % 脉冲持续时间为0.3ms
    start_freq = 2e6;  % 起始频率为2MHz
    end_freq = 4e6;    % 结束频率为4MHz
    bandwidth = end_freq-start_freq;  % 带宽为2MHz
    t = 0:1/fs:t_pulse-1/fs;    % 生成时间向量
    
    chirp_signal = exp(1i * 2 * pi * (start_freq * t + 0.5*1/t_pulse * bandwidth * t.^2));
    
    start_estimate=[];              %  起始时间估计RMSE/SNR
    continue_estimate=[];           %  持续时间估计RMSE/SNR
    bandwidth_estimate=[];          %  带宽估计RMSE/SNR
    freq_estimate=[];               %  起始频率估计RMSE/SNR

    SNR=0:1:35;
    for snr = SNR
        % 对于每个信噪比进行50次蒙洛卡特实验
        start_estimate1=[];              %  起始时间估计误差/kkk
        continue_estimate1=[];           %  持续时间估计误差/kkk
        bandwidth_estimate1=[];          %  带宽估计误差/kkk
        freq_estimate1=[];               %  起始频率估计误差/kkk
        kkkkk=50;
        for kkk=1:kkkkk
            % 在信号前后补零以满足时间为1ms
            lower_limit = 0;  % 下限
            upper_limit = round((time - t_pulse)*fs);  % 上限
            begin_index=randi([lower_limit, upper_limit])+1; % chirp信号的起始Index
            front_zero = zeros(1 , begin_index-1);
            back_zero = zeros(1 , upper_limit-begin_index+1 );
            extended_signal=[front_zero,chirp_signal,back_zero];

            observe_signal=awgn(extended_signal,snr,"measured");  % 观测时间信号
            
            %figure;
            %plot(t_all, real(observe_signal));
            %title('观测信号时域波形');
            %xlabel('时间（秒）');
            %ylabel('幅度');
        
        % 找到chirp信号的起始时间
        
            window_size =128;  % 窗口大小
            overlap = 0.9;      % 重叠比例为90%，则时频图时间两点之差为(window_size-round(overlap*window_size))
            nfft = 2^(nextpow2(window_size)+1);    % FFT点数
            [S, F, T] = spectrogram(observe_signal, window_size, round(overlap*window_size), nfft, fs);
            
%             figure;
%             imagesc(T, F, 10*log10(abs(S)));
%             axis xy;
%             xlabel('时间 (s)');
%             ylabel('频率 (Hz)');
%             title('Chirp信号的短时傅里叶变换');
%             colormap('jet');
        
            energy = sum(abs(S).^2) / nfft;
            diff_energy = diff(energy);    % 能量差分值
            [peaks,peaks_index]=findpeaks(abs(diff_energy));
            threshold = max(5*mean(abs(diff_energy)),0.7*max(peaks));  % 设置差分阈值，用于检测脉冲
            pulse_indices = find(abs(diff_energy) > threshold);        % 检测脉冲出现的时间

            signal_index=round(T(pulse_indices)*fs);   
            if ~isempty(signal_index) && ~(max(signal_index)-min(signal_index)<100)
                start_index=min(signal_index);    % chirp信号起始Index
                end_index=max(signal_index);      % chirp信号结束Index
            else
                % 较低信噪比
                [maxx,maxx_index]=max(energy);
                threshold= 0.85*maxx;
                pulse_indices = find(abs(energy) > threshold); 
                signal_index=round(T(pulse_indices)*fs);
                start_index=min(signal_index);    % chirp信号起始Index
                end_index=max(signal_index);      % chirp信号结束Index
            end
            
            start_index=start_index+17;
            end_index=end_index-3;              %对结果进行微调

        % 估计起始频率及带宽，先将信号取出，变为xn=Aexp(j*2*pi*(1/2αn^2+fn))+en，其中B=N*α*fs，start_freq=f*fs
        % 参考论文Parameter Estimation of Chirp Signals,使用相位展开对α与f进行估计
        
            discrete_signal=observe_signal(start_index:end_index);
            N=length(discrete_signal);
            xn=discrete_signal;
            yn=xn.*conj([0,xn(1:end-1)]);
            zn=yn.*conj([0,yn(1:end-1)]);
            phi_2=angle(zn);
            phi=zeros(size(phi_2));
            phi(1)=angle(xn(1));
            phi(2)=angle(yn(2))+phi(1);
            for i=3:N
                phi(i)=phi_2(i)+2*phi(i-1)-phi(i-2);
            end
            G = [ones(N,1),(1:N)',((1:N).^2)'];
            theta = (inv(G'*G)*G')*phi';
            B = mod(N*theta(3)*fs/(pi),fs);
            start_freq_estimate = mod(theta(2)*fs/(2*pi),fs);  % 初始频率需要限制在fs内
                       
            start_estimate1=[start_estimate1,(start_index-begin_index)/fs];
            continue_estimate1=[continue_estimate1,(end_index-start_index)/fs-t_pulse];
            bandwidth_estimate1=[bandwidth_estimate1,B-bandwidth];          
            freq_estimate1=[freq_estimate1,start_freq_estimate-start_freq];   
%             start_estimate1=[start_estimate1,(start_index-begin_index)/begin_index];
%             continue_estimate1=[continue_estimate1,((end_index-start_index)/fs-t_pulse)/t_pulse];
%             bandwidth_estimate1=[bandwidth_estimate1,(B-bandwidth)/bandwidth];          
%             freq_estimate1=[freq_estimate1,(start_freq_estimate-start_freq)/start_freq];   
%             start_estimate1=[start_estimate1,(start_index-begin_index)];
%             continue_estimate1=[continue_estimate1,(end_index-start_index)];
%             bandwidth_estimate1=[bandwidth_estimate1,B];          
%             freq_estimate1=[freq_estimate1,start_freq_estimate];  
        
        
        end
        start_estimate=[start_estimate,sqrt(sum(start_estimate1.^2)/kkkkk)];              %  起始时间估计RMSE/SNR
        continue_estimate=[continue_estimate,sqrt(sum(continue_estimate1.^2)/kkkkk)];           %  持续时间估计RMSE/SNR
        bandwidth_estimate=[bandwidth_estimate,sqrt(sum(bandwidth_estimate1.^2)/kkkkk)];          %  带宽估计RMSE/SNR
        freq_estimate=[freq_estimate,sqrt(sum(freq_estimate1.^2)/kkkkk)];               %  起始频率估计RMSE/SNR

%         start_estimate=[start_estimate,(sum(start_estimate1)/kkkkk)];             
%         continue_estimate=[continue_estimate,(sum(continue_estimate1)/kkkkk)];          
%         bandwidth_estimate=[bandwidth_estimate,(sum(bandwidth_estimate1)/kkkkk)];        
%         freq_estimate=[freq_estimate,(sum(freq_estimate1)/kkkkk)];              
    end

    time_points = SNR;

    % 画出 start_estimate 的图
    figure;
    plot(time_points, start_estimate, '-o');
    title('起始时间估计');
    xlabel('SNR');
    ylabel('RMSE');
    
    % 画出 continue_estimate 的图
    figure;
    plot(time_points, continue_estimate, '-o');
    title('持续时间估计');
    xlabel('SNR');
    ylabel('RMSE');
    
    % 画出 bandwidth_estimate 的图
    figure;
    plot(time_points, bandwidth_estimate, '-o');
    title('带宽估计');
    xlabel('SNR');
    ylabel('RMSE');
    
    % 画出 freq_estimate 的图
    figure;
    plot(time_points, freq_estimate, '-o');
    title('起始频率估计');
    xlabel('SNR');
    ylabel('RMSE');
        


