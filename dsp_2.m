clear all;

data= load('wave_data.mat');
wave_data=data.wave_data;


% 参数设置
fs = 50e6;          % 采样率为50MHz
t_pulse = 10e-6;    % 脉冲持续时间为10μs
t = 0:1/fs:t_pulse-1/fs;  % 生成时间向量

start_freq = 15e6;  % 起始频率为15MHz
end_freq = 35e6;    % 结束频率为35MHz
bandwidth = 20e6;  % 带宽为20MHz

start_freq1 = 20e6;  % 起始频率为20MHz
end_freq1 = 32e6;    % 结束频率为32MHz
bandwidth1 = 12e6;  % 带宽为12MHz

% 生成复数“chirp”信号
chirp_signal = exp(1i * 2 * pi * (start_freq * t + 0.5*1/t_pulse * bandwidth * t.^2));
chirp_signal1 = exp(1i * 2 * pi * (start_freq1 * t + 0.5*1/t_pulse * bandwidth1 * t.^2));

new_data= chirp_signal1+chirp_signal+wave_data;

% 进行STFT
window_size =34;  % 窗口大小
overlap = 0.9;      % 重叠比例为50%
nfft = 256;         % FFT点数

[S, F, T] = spectrogram(new_data, window_size, round(overlap*window_size), nfft, fs);

% 绘制STFT结果

%  在固定重叠比例的情况下，winodw_size在30-40既能分辨wave_data的变化趋势可以比较清晰地分辨，
%  也在较长的时间里可以分辨两个复指数信号，窗更小时由于采样点不够，会发生较明显的栅栏效应，使得
%  难以分辨两个复指数信号的频谱，窗更大时两个复指数信号的宽度更窄，更易区分，但由于窗长增加，
%  不同时间的频率被乘以一个相位再相加，在不同时间的相同频率会发生混叠，使得信号的时频特性无法区分

%  wave_data的

figure;
imagesc(T, F, 10*log10(abs(S)));
axis xy;
xlabel('时间 (s)');
ylabel('频率 (Hz)');
title('Chirp信号的短时傅里叶变换');
colormap('jet');
