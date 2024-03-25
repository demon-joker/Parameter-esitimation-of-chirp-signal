% 参数设置
fs = 50e6;          % 采样率为50MHz
t_pulse = 10e-6;    % 脉冲持续时间为10μs
t = 0:1/fs:t_pulse-1/fs;  % 生成时间向量

start_freq = 15e6;  % 起始频率为15MHz
end_freq = 35e6;    % 结束频率为35MHz
bandwidth = 20e6;  % 带宽为20MHz

% 生成复数“chirp”信号
chirp_signal = exp(1i * 2 * pi * (start_freq * t + 0.5*1/t_pulse * bandwidth * t.^2));

% 进行STFT
window_size = 60;  % 窗口大小
overlap = 55;      % 重叠比例为50%
nfft = 256;         % FFT点数

[S, F, T] = spectrogram(chirp_signal, window_size, overlap, nfft, fs);

% 绘制STFT结果
figure;


plot(t, real(chirp_signal));
title('Chirp信号时域波形');
xlabel('时间（秒）');
ylabel('幅度');


figure;
imagesc(T, F, 10*log10(abs(S)));
axis xy;
xlabel('时间 (s)');
ylabel('频率 (Hz)');
title('Chirp信号的短时傅里叶变换');
colormap('jet');
