clear all;

rng(111);

% 参数设置
fs = 1e7;          % 采样率为10MHz
time = 1e-3;       % 整段时间为1ms
t_all=0:1/fs:time-1/fs; 

% 生成观测信号，在上一题中，设置信号起始时间为随机值
% 在本题中固定信号数与起始时间，且对于估计算法来说信号数已知,信号数为2
% 使用分数阶傅里叶变换(FRFT)在时频平面搜索
% FRFT 函数参考 https://github.com/nadavleh/FrFT/blob/master/frft.m

    t_pulse_1 = 3e-4;   % 脉冲持续时间为0.3ms
    start_freq_1 = 1e6;  % 起始频率为5MHz
    end_freq_1 = 4e6;    % 结束频率为2MHz
    bandwidth_1 = end_freq_1-start_freq_1;  % 带宽为-3MHz
    t_1 = 0:1/fs:t_pulse_1-1/fs;    % 生成时间向量
    chirp_signal_1 = 0.12*exp(1i * 2 * pi * (start_freq_1 * t_1 + 0.5*1/t_pulse_1 * bandwidth_1 * t_1.^2));

    t_pulse_2 = 4.5e-4;   % 脉冲持续时间为0.45ms
    start_freq_2 = 5e6;  % 起始频率为2MHz
    end_freq_2 = 2e6;    % 结束频率为4MHz
    bandwidth_2 = end_freq_2-start_freq_2;  % 带宽为2MHz
    t_2 = 0:1/fs:t_pulse_2-1/fs;    % 生成时间向量
    chirp_signal_2 = 0.12*exp(1i * 2 * pi * (start_freq_2 * t_2 + 0.5*1/t_pulse_2 * bandwidth_2 * t_2.^2));

    begin_index_1 = 4000;         % chirp_1信号的起始Index
    front_zero = zeros(1 , begin_index_1-1);
    back_zero = zeros(1 , round((time - t_pulse_1)*fs)-begin_index_1+1 );
    extended_signal_1=[front_zero,chirp_signal_1,back_zero];

    begin_index_2 = 3000;         % chirp_2信号的起始Index
    front_zero = zeros(1 , begin_index_2-1);
    back_zero = zeros(1 , round((time - t_pulse_2)*fs)-begin_index_2+1 );
    extended_signal_2=[front_zero,chirp_signal_2,back_zero];

    % 由于识别并估计的参数并不需要和原信号的序号对应，
    % 采用两个信号的参数估计值的和与差的绝对值计算RMSE，这样也可以解出相同的两个参数
    % 后缀1代表和，2代表差的绝对值

    start_estimate_RMSE_1=[];              %  起始时间估计RMSE/SNR
    continue_estimate_RMSE_1=[];           %  持续时间估计RMSE/SNR
    bandwidth_estimate_RMSE_1=[];          %  带宽估计RMSE/SNR
    freq_estimate_RMSE_1=[];               %  起始频率估计RMSE/SNR
    start_estimate_RMSE_2=[];              %  起始时间估计RMSE/SNR
    continue_estimate_RMSE_2=[];           %  持续时间估计RMSE/SNR
    bandwidth_estimate_RMSE_2=[];          %  带宽估计RMSE/SNR
    freq_estimate_RMSE_2=[];               %  起始频率估计RMSE/SNR

    SNR=0:1.5:15;
for snr = SNR

    start_estimate_1=[];              %  起始时间估计误差
    continue_estimate_1=[];           %  持续时间估计误差
    bandwidth_estimate_1=[];          %  带宽估计误差
    freq_estimate_1=[];               %  起始频率估计误差
    start_estimate_2=[];              %  起始时间估计误差
    continue_estimate_2=[];           %  持续时间估计误差
    bandwidth_estimate_2=[];          %  带宽估计误差
    freq_estimate_2=[];               %  起始频率估计误差

    
    kkkkk = 3;
    for kkk = 1:kkkkk % 由于计算时间较长，对于每个信噪比进行3次蒙特卡洛实验
    observe_signal=awgn(extended_signal_1+extended_signal_2,snr,"measured");  % 观测时间信号
    
%     figure;
%     plot(t_all, real(observe_signal));
%     title('观测信号时域波形');
%     xlabel('时间（秒）');
%     ylabel('幅度');
%
%                         observe_signal作为算法输入，以下为估计算法
%

    window_size =128;  % 窗口大小
    overlap = 0.9;      % 重叠比例为90%
    nfft = 2^(nextpow2(window_size)+1);    % FFT点数
    [S, F, T] = spectrogram(observe_signal, window_size, round(overlap*window_size), nfft, fs);

%     figure;
%     surf(T, F, abs(S), 'EdgeColor', 'none');
%     axis xy;
%     xlabel('时间 (s)');
%     ylabel('频率 (Hz)');
%     title('Chirp信号的短时傅里叶变换');
%     colormap('jet');
    
    k=mean(abs(S));
    S(abs(S)<3*k)=0.01*mean(k);  % 直接滤除大部分噪声

    % 若两信号在时域上存在重叠，估计最早的起始时间和总持续时间
    % 若两信号在时域上不重叠，则直接估计两个信号的起始时间与持续时间
    % 对于-5及以上的高信噪比，两信号时域上不重合要达到200-300点才可从能量上分辨出来
    % 对于更低的信噪比，很可能会出现干扰峰，无法完美分辨

    energy = sum(abs(S).^2) / nfft;   
    energy1=energy;
    energy=10*log10(energy);  % 让能量与噪声之间差别变大
    middle_window = 15;
    energy = medfilt1(energy, middle_window); % 中值滤波滤除噪声产生的脉冲分量
    
    diff_energy = diff(energy);    % 能量差分值
    [peaks,peaks_index]=findpeaks(abs(diff_energy));
    threshold = 0.5*max(peaks);                     % 设置差分阈值，用于检测脉冲,
    % 在高信噪比(如大于20dB)时，该阈值可能会过高而出现时间分割错误
    pulse_indices = peaks_index(find( peaks > threshold));        % 检测脉冲出现的时间
    pulse_value = diff_energy(pulse_indices);
    
%     figure;
%     subplot(3,1,1);
%     plot(T, energy1, 'b');
%     xlabel('时间 (s)');
%     ylabel('能量');
%     title('信号能量随时间的变化');
%     subplot(3,1,2);
%     plot(T, energy, 'b');
%     xlabel('时间 (s)');
%     ylabel('能量');
%     title('信号能量随时间的变化(滤波后)');
%     subplot(3,1,3);
%     plot(T(2:end), diff_energy, 'r');
%     hold on;
%     plot(T(pulse_indices), diff_energy(pulse_indices), 'go', 'MarkerSize', 8, 'DisplayName', '脉冲');
%     xlabel('时间 (s)');
%     ylabel('能量差分值');
%     title('能量差分值随时间的变化');
    
    signal_overlap = 0; % 判断两信号是否在时域重叠
    start_index_of_two = 0;
    end_index_of_two = 0;
    end_of_signal_1 = 0;
    start_of_siganl_2 = 0;

    signal_index=round(T(pulse_indices)*fs);   
    lkl = length(signal_index);
    if  lkl > 2 
        % 确定起始点
        this_index = 0;
        for jj=1:lkl-1
            if abs(signal_index(jj)-signal_index(jj+1))>300 && (pulse_value(jj)>0) && (pulse_value(jj+1)<0)      % 否则认为该脉冲为噪声
                start_index_of_two = signal_index(jj)+2*middle_window;
                this_index = jj;
                break;
            end
        end
        signal_index = signal_index(this_index:end);
        pulse_value = pulse_value(this_index:end);
        lkl = length(signal_index);
        % 确定终点
        for jj=1:lkl-1
            if abs(signal_index(lkl-jj+1)-signal_index(lkl-jj))>300 && (pulse_value(lkl-jj+1)<0) && (pulse_value(lkl-jj)>0)      % 否则认为该脉冲为噪声
                end_index_of_two = signal_index(lkl-jj+1)-2*middle_window;
                this_index = lkl-jj+1;
                break;
            end
        end
        signal_index = signal_index(1:this_index);
        pulse_value = pulse_value(1:this_index);
        lkl = length(signal_index); 
        if lkl >= 4 % 认为两信号分离
            end_of_signal_1 = signal_index(2)-2*middle_window;
            start_of_siganl_2 = signal_index(end-1)+2*middle_window;
            signal_overlap = 0;
        else   % 认为两信号存在时域重叠
            signal_overlap = 1;
        end
    else  % 只有两个脉冲，两信号存在时域重叠
        start_index_of_two = signal_index(1)+2*middle_window;
        end_index_of_two = signal_index(end)-2*middle_window;
        signal_overlap = 1;
    end
    
    if signal_overlap % 存在时域重叠，进一步分离两个信号的时间
        discrete_signal = observe_signal(start_index_of_two:end_index_of_two);

        window_size =128;  % 窗口大小
        overlap = 0.9;      % 重叠比例为90%
        nfft = 2^(nextpow2(window_size)+1);    % FFT点数
        [S, F, T] = spectrogram(discrete_signal, window_size, round(overlap*window_size), nfft, fs);
    
%         figure;
%         surf(T, F, abs(S), 'EdgeColor', 'none');
%         axis xy;
%         xlabel('时间 (s)');
%         ylabel('频率 (Hz)');
%         title('Chirp信号的短时傅里叶变换');
%         colormap('jet');

        k=mean(abs(S));
        S(abs(S)<3*k)=0;  % 直接滤除大部分噪声
 
        % 找出每个T对应的信号数
        num_of_signal=zeros(size(T));
        for jkjk=1:length(T)
            [peaks]=findpeaks(abs(S(:,jkjk)));
            num_of_signal(jkjk) = length(peaks);  
        end
        middle_window = 15;
        num_of_signal = medfilt1(num_of_signal, middle_window); % 中值滤波
        diff_num = diff(num_of_signal);    % 差分值
        [peaks,peaks_index]=findpeaks(abs(diff_num));
        peak_value = diff_num(peaks_index);
        sds = abs(length(peaks)-2); 

        % 分离重叠部分与非重叠部分，分开估计参数
        if sds  % 信号数只有一次突变
            if peak_value(1)>0
                overlap_begin_index=round(T(peaks_index(1))*fs)+start_index_of_two+2*middle_window;
                overlap_end_index=end_index_of_two;
                unoverlap_signal = observe_signal(start_index_of_two:overlap_begin_index-1);
            else
                overlap_begin_index=start_index_of_two;
                overlap_end_index=round(T(peaks_index(1))*fs)+start_index_of_two-2*middle_window;
                unoverlap_signal = observe_signal(overlap_end_index+1:end_index_of_two);
            end
            % 对没有重叠的部分估计参数
            [start_freq_estimate_unoverlap,...
             end_freq_estimate_unoverlap,...
             B_estimate_unoverlap]=estimate_overlap_chirp(unoverlap_signal,fs,0);
        
        else  % 信号数有两次突变
            overlap_begin_index=round(T(peaks_index(1))*fs)+start_index_of_two+2*middle_window;
            overlap_end_index=round(T(peaks_index(2))*fs)+start_index_of_two-2*middle_window;

            unoverlap_signal_1 = observe_signal(start_index_of_two:overlap_begin_index-1);
            unoverlap_signal_2 = observe_signal(overlap_end_index+1:end_index_of_two);
            % 估计未重叠信号的参数
            [start_freq_estimate_unoverlap_1,...
                end_freq_estimate_unoverlap_1,...
             B_estimate_unoverlap_1]=estimate_overlap_chirp(unoverlap_signal_1,fs,0);

            [start_freq_estimate_unoverlap_2,...
             end_freq_estimate_unoverlap_2,...
             B_estimate_unoverlap_2]=estimate_overlap_chirp(unoverlap_signal_2,fs,0);
        end

        % 使用FRFT估计重叠部分的参数
        % B_estimate = -fs*cot(alpha) =- fs*cot(p*pi/2),
        % middle_freq_estimate = mod(u_estimate*csc(alpha*pi/2),fs)
        overlap_signal = observe_signal(overlap_begin_index:overlap_end_index);
         
        [overlap_begin_freq_estimate,...
         overlap_end_freq_estimate,...
         B_estimate_overlap]=estimate_overlap_chirp(overlap_signal,fs,1);
        
        alpha_1 = B_estimate_overlap(1)/(length(overlap_signal)*fs);
        alpha_2 = B_estimate_overlap(2)/(length(overlap_signal)*fs);

        if sds   % 信号数只有一次突变，,这里信号1为重叠信号(1)，信号2为重叠信号(2)
            alpha_3 = B_estimate_unoverlap(1)/(length(unoverlap_signal)*fs);
            if peak_value(1)>0 % 突变为正，未重叠信号在重叠信号左边
                error_1 = abs(end_freq_estimate_unoverlap(1) - overlap_begin_freq_estimate(1));
                error_2 = abs(end_freq_estimate_unoverlap(1) - overlap_begin_freq_estimate(2));

                if ( abs((abs(alpha_3-alpha_2) - abs(alpha_3-alpha_1))) < 0.05*abs(alpha_3-alpha_2)...
                     && error_2 < error_1 ) ||  abs(alpha_3-alpha_1) > 1.05*abs(alpha_3-alpha_2)% 结束频率与开始频率差不多且调频斜率差不多
                    % 连信号2
                    B_estimate_2 = B_estimate_unoverlap(1) + B_estimate_overlap(2);
                    B_estimate_1 = B_estimate_overlap(1);
                    start_freq_estimate_2 = start_freq_estimate_unoverlap(1);
                    start_freq_estimate_1 = overlap_begin_freq_estimate(1);
                    start_index_estimate_2 = start_index_of_two;
                    start_index_estimate_1 = overlap_begin_index;
                    continue_index_estimate_2 = overlap_end_index - start_index_of_two;
                    continue_index_estimate_1 = overlap_end_index - overlap_begin_index;
                else %连信号1
                    B_estimate_1 = B_estimate_unoverlap(1) + B_estimate_overlap(1);
                    B_estimate_2 = B_estimate_overlap(2);
                    start_freq_estimate_1 = start_freq_estimate_unoverlap(1);
                    start_freq_estimate_2 = overlap_begin_freq_estimate(2);
                    start_index_estimate_1 = start_index_of_two;
                    start_index_estimate_2 = overlap_begin_index;
                    continue_index_estimate_1 = overlap_end_index - start_index_of_two;
                    continue_index_estimate_2 = overlap_end_index - overlap_begin_index;
                end
           
            else    %  突变为负，未重叠信号在重叠信号右边
                error_1 = abs(start_freq_estimate_unoverlap(1) - overlap_end_freq_estimate(1));
                error_2 = abs(start_freq_estimate_unoverlap(1) - overlap_end_freq_estimate(2));
                    
                if ( abs((abs(alpha_3-alpha_2) - abs(alpha_3-alpha_1))) < 0.05*abs(alpha_3-alpha_2)...
                     && error_2 < error_1 ) ||  abs(alpha_3-alpha_1) > 1.05*abs(alpha_3-alpha_2)% 结束频率与开始频率差不多且调频斜率差不多
                    %连信号2
                    B_estimate_2 = B_estimate_unoverlap(1) + B_estimate_overlap(2);
                    B_estimate_1 = B_estimate_overlap(1);
                    start_freq_estimate_1 = overlap_begin_freq_estimate(1);
                    start_freq_estimate_2 = overlap_begin_freq_estimate(2);
                    start_index_estimate_1 = overlap_begin_index;
                    start_index_estimate_2 = overlap_begin_index;
                    continue_index_estimate_2 = end_index_of_two - overlap_begin_index;
                    continue_index_estimate_1 = overlap_end_index - overlap_begin_index;
                else  %连信号1
                    B_estimate_1 = B_estimate_unoverlap(1) + B_estimate_overlap(1);
                    B_estimate_2 = B_estimate_overlap(2);
                    start_freq_estimate_1 = overlap_begin_freq_estimate(1);
                    start_freq_estimate_2 = overlap_begin_freq_estimate(2);
                    start_index_estimate_1 = overlap_begin_index;
                    start_index_estimate_2 = overlap_begin_index;
                    continue_index_estimate_1 = end_index_of_two - overlap_begin_index;
                    continue_index_estimate_2 = overlap_end_index - overlap_begin_index;
                end
            end

        else % 信号数有两次突变，可视为上述两种情况的级联,这里信号1为重叠信号(1)，信号2为重叠信号(2)

            % 先是突变为正,
                alpha_3 = B_estimate_unoverlap_1(1)/(length(unoverlap_signal_1)*fs);
                error_1 = abs(end_freq_estimate_unoverlap_1(1) - overlap_begin_freq_estimate(1));
                error_2 = abs(end_freq_estimate_unoverlap_1(1) - overlap_begin_freq_estimate(2));

                if ( abs((abs(alpha_3-alpha_2) - abs(alpha_3-alpha_1))) < 0.05*abs(alpha_3-alpha_2)...
                     && error_2 < error_1 ) ||  abs(alpha_3-alpha_1) > 1.05*abs(alpha_3-alpha_2)% 结束频率与开始频率差不多且调频斜率差不多
                    
                    B_estimate_2 = B_estimate_unoverlap_1(1) + B_estimate_overlap(2);
                    B_estimate_1 = B_estimate_overlap(1);
                    start_freq_estimate_2 = start_freq_estimate_unoverlap_1(1);
                    start_freq_estimate_1 = overlap_begin_freq_estimate(1);
                    start_index_estimate_1 = overlap_begin_index;
                    start_index_estimate_2 = start_index_of_two;
                    continue_index_estimate_1 = overlap_end_index - overlap_begin_index;
                    continue_index_estimate_2 = overlap_end_index - start_index_of_two;
                else
                    B_estimate_1 = B_estimate_unoverlap_1(1) + B_estimate_overlap(1);
                    B_estimate_2 = B_estimate_overlap(2);
                    start_freq_estimate_1 = start_freq_estimate_unoverlap_1(1);
                    start_freq_estimate_2 = overlap_begin_freq_estimate(2);
                    start_index_estimate_1 = start_index_of_two;
                    start_index_estimate_2 = overlap_begin_index;
                    continue_index_estimate_1 = overlap_end_index - start_index_of_two;
                    continue_index_estimate_2 = overlap_end_index - overlap_begin_index;
                end
                    
           
            % 然后突变为负，则可把前面的估计参数视为重叠信号的估计参数
                alpha_3 = B_estimate_unoverlap_2(1)/(length(unoverlap_signal_2)*fs);
                error_1 = abs(start_freq_estimate_unoverlap_2(1) - overlap_end_freq_estimate(1));
                error_2 = abs(start_freq_estimate_unoverlap_2(1) - overlap_end_freq_estimate(2));

                if ( abs((abs(alpha_3-alpha_2) - abs(alpha_3-alpha_1))) < 0.05*abs(alpha_3-alpha_2)...
                     && error_2 < error_1 ) ||  abs(alpha_3-alpha_1) > 1.05*abs(alpha_3-alpha_2)  % 结束频率与开始频率差不多且调频斜率差不多
                    % 连信号2
                    B_estimate_2 = B_estimate_unoverlap_2(1) + B_estimate_2;
                    B_estimate_1 = B_estimate_1;                            % 不变
                    start_freq_estimate_1 = start_freq_estimate_1;          % 不变
                    start_freq_estimate_2 = start_freq_estimate_2;          % 不变
                    start_index_estimate_1 = start_index_estimate_1;        % 不变
                    start_index_estimate_2 = start_index_estimate_2;        % 不变
                    continue_index_estimate_1 = continue_index_estimate_1;  % 不变
                    continue_index_estimate_2 = continue_index_estimate_2 + end_index_of_two - overlap_end_index;
                else  %连信号1
                    B_estimate_1 = B_estimate_unoverlap(1) + B_estimate_1;
                    B_estimate_2 = B_estimate_overlap(2);
                    start_freq_estimate_1 = start_freq_estimate_1;          % 不变
                    start_freq_estimate_2 = start_freq_estimate_2;          % 不变
                    start_index_estimate_1 = start_index_estimate_1;        % 不变
                    start_index_estimate_2 = start_index_estimate_2;        % 不变
                    continue_index_estimate_1 = continue_index_estimate_1 + end_index_of_two - overlap_end_index;
                    continue_index_estimate_2 = continue_index_estimate_2;  % 不变
                end
        end
    else % 两信号时域分离，分别估计参数

        discrete_signal_1=observe_signal(start_index_of_two:end_of_signal_1);
        discrete_signal_2=observe_signal(start_of_siganl_2:end_index_of_two);
        
        [start_freq_estimate_1,...
         end_freq_estimate_1,...
         B_estimate_1]=estimate_overlap_chirp(discrete_signal_1,fs,0);
    
        [start_freq_estimate_2,...
         end_freq_estimate_2,...
         B_estimate_2]=estimate_overlap_chirp(discrete_signal_2,fs,0);

        start_index_estimate_1 = start_index_of_two;      
        start_index_estimate_2 = start_of_siganl_2;        
        continue_index_estimate_1 = end_of_signal_1 -start_index_of_two;
        continue_index_estimate_2 = end_index_of_two-start_of_siganl_2;  

    end
    % 已得到对于两个信号的参数估计值，计算RMSE
    
    
    start_estimate_1=[start_estimate_1,(start_index_estimate_1+start_index_estimate_2-(begin_index_2+begin_index_1))/fs];              %  起始时间估计误差
    continue_estimate_1=[continue_estimate_1,(continue_index_estimate_1+continue_index_estimate_2)/fs-(t_pulse_1+t_pulse_2)];           %  持续时间估计误差
    bandwidth_estimate_1=[bandwidth_estimate_1,(B_estimate_1+B_estimate_2)-(bandwidth_1+bandwidth_2)];          %  带宽估计误差
    freq_estimate_1=[freq_estimate_1,(start_freq_estimate_1+start_freq_estimate_2)-(start_freq_1+start_freq_2)];               %  起始频率估计误差
    start_estimate_2=[start_estimate_2,(abs(start_index_estimate_1-start_index_estimate_2)-abs(begin_index_2-begin_index_1))/fs];              %  起始时间估计误差
    continue_estimate_2=[continue_estimate_2,abs(continue_index_estimate_1-continue_index_estimate_2)/fs-abs(t_pulse_1-t_pulse_2)];           %  持续时间估计误差
    bandwidth_estimate_2=[bandwidth_estimate_2,abs(B_estimate_1-B_estimate_2)-abs(bandwidth_1-bandwidth_2)];          %  带宽估计误差
    freq_estimate_2=[freq_estimate_2,abs(start_freq_estimate_1-start_freq_estimate_2)-abs(start_freq_1-start_freq_2)];               %  起始频率估计误差

    end

    start_estimate_RMSE_1=[start_estimate_RMSE_1,sqrt(sum(start_estimate_1.^2)/kkkkk)];              %  起始时间估计RMSE/SNR
    continue_estimate_RMSE_1=[continue_estimate_RMSE_1,sqrt(sum(continue_estimate_1.^2)/kkkkk)];           %  持续时间估计RMSE/SNR
    bandwidth_estimate_RMSE_1=[bandwidth_estimate_RMSE_1,sqrt(sum(bandwidth_estimate_1.^2)/kkkkk)];          %  带宽估计RMSE/SNR
    freq_estimate_RMSE_1=[freq_estimate_RMSE_1,sqrt(sum(freq_estimate_1.^2)/kkkkk)];               %  起始频率估计RMSE/SNR
    start_estimate_RMSE_2=[start_estimate_RMSE_2,sqrt(sum(start_estimate_2.^2)/kkkkk)];              %  起始时间估计RMSE/SNR
    continue_estimate_RMSE_2=[continue_estimate_RMSE_2,sqrt(sum(continue_estimate_2.^2)/kkkkk)];           %  持续时间估计RMSE/SNR
    bandwidth_estimate_RMSE_2=[bandwidth_estimate_RMSE_2,sqrt(sum(bandwidth_estimate_2.^2)/kkkkk)];          %  带宽估计RMSE/SNR
    freq_estimate_RMSE_2=[freq_estimate_RMSE_2,sqrt(sum(freq_estimate_2.^2)/kkkkk)];               %  起始频率估计RMSE/SNR

end

    time_points = SNR;

    % 画图
    figure;
    subplot(2, 2, 1);
    plot(time_points, start_estimate_RMSE_1, 'o-', 'LineWidth', 1);
    title('起始时间估计RMSE/SNR - 和');
    subplot(2, 2, 2);
    plot(time_points, continue_estimate_RMSE_1, 'o-', 'LineWidth', 1);
    title('持续时间估计RMSE/SNR - 和');
    subplot(2, 2, 3);
    plot(time_points, bandwidth_estimate_RMSE_1, 'o-', 'LineWidth', 1);
    title('带宽估计RMSE/SNR - 和');
    subplot(2, 2, 4);
    plot(time_points, freq_estimate_RMSE_1, 'o-', 'LineWidth', 1);
    title('起始频率估计RMSE/SNR - 和');
    
    figure;
    subplot(2, 2, 1);
    plot(time_points, start_estimate_RMSE_2, 'o-', 'LineWidth', 1);
    title('起始时间估计RMSE/SNR - 差的绝对值');
    subplot(2, 2, 2);
    plot(time_points, continue_estimate_RMSE_2, 'o-', 'LineWidth', 1);
    title('持续时间估计RMSE/SNR - 差的绝对值');
    subplot(2, 2, 3);
    plot(time_points, bandwidth_estimate_RMSE_2, 'o-', 'LineWidth', 1);
    title('带宽估计RMSE/SNR - 差的绝对值');
    subplot(2, 2, 4);
    plot(time_points, freq_estimate_RMSE_2, 'o-', 'LineWidth', 1);
    title('起始频率估计RMSE/SNR - 差的绝对值');





function [overlap_begin_freq_estimate,...
          overlap_end_freq_estimate,...
          B_estimate_overlap]=estimate_overlap_chirp(overlap_signal,fs,overlap)
    % 只能识别符合奈奎斯特采样率的信号，更高频的信号会发生混叠
    fff = linspace(-fs/2,fs/2-1,length(overlap_signal)); 
    the_max=[];
    P=[-1:0.001:-0.5,0.5:0.001:1];
    for p=P
        frft_signal = frft(overlap_signal, p);
        [maxx,~]=max(abs(frft_signal));
        the_max=[the_max,maxx];
    end
%     figure;
%     plot(P,the_max);
    [peaks,peak_index] = findpeaks(the_max);
    if overlap
        ashgah = peak_index(peaks>0.6*max(peaks));
    else
        ashgah = peak_index(peaks == max(peaks));
    end
    p_chosen=P(ashgah);
    B_estimate_overlap=[];
    for pp = p_chosen
        B_estimate_overlap =[B_estimate_overlap, - fs*cot(pp*pi/2)];
    end
    f_middle = [];
    for ssss=p_chosen
        frft_signal = frft(overlap_signal, ssss);
%         figure;
%         plot(fff,frft_signal);
%         title("frft");
        [peaks,peaks_index]=findpeaks(abs(frft_signal));
        maxx_index = peaks_index(peaks>0.7*max(peaks));    % 调频斜率相同时，可能有两个多个peak
        u_estimate=fff(maxx_index);
        f_middle=[f_middle, mod(u_estimate*csc(ssss*pi/2),fs)];
    end
    disp(f_middle);
    disp(B_estimate_overlap);
    if length(B_estimate_overlap) < length(f_middle) % 只可能为调频斜率一样，但中心频率不同
        B_estimate_overlap=[B_estimate_overlap(1),B_estimate_overlap(1)];
    end
        overlap_end_freq_estimate = f_middle + 1/2*B_estimate_overlap;
        overlap_begin_freq_estimate = f_middle - 1/2*B_estimate_overlap;
end


% function [B,start]=estimate_single_chirp(signal,fs)
%     N=length(signal);
%     xn=signal;
%     yn=xn.*conj([0,xn(1:end-1)]);
%     zn=yn.*conj([0,yn(1:end-1)]);
%     phi_2=angle(zn);
%     phi=zeros(size(phi_2));
%     phi(1)=angle(xn(1));
%     phi(2)=angle(yn(2))+phi(1);
%     for i=3:N
%         phi(i)=phi_2(i)+2*phi(i-1)-phi(i-2);
%     end
%     G = [ones(N,1),(1:N)',((1:N).^2)'];
%     theta = (inv(G'*G)*G')*phi';
%     B =  mod(N*theta(3)*fs/(pi),fs);
%     start = mod(theta(2)*fs/(2*pi),fs);
% end




         