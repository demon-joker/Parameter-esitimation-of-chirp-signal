DSP大作业报告

### 实现chirp信号并用TFT验证

chirp信号的形式为
$$
x(t)=A\exp({j(2\pi(ft+\frac{1}{2}wt^2)+\phi)})\;\;\;t\in[0,T]
$$
其中 $f$ 为初始频率，$w$ 为调频速率，有带宽 $B=wT$ 成立

本次作业中不考虑初相位 $\phi$ ，故 $\phi=0$

则可用以下代码生成复数chirp信号，

~~~matlab
start_freq = 15e6;  % 起始频率为15MHz
bandwidth = 20e6;  % 带宽为20MHz
chirp_signal = exp(1i * 2 * pi * (start_freq * t + 0.5*1/t_pulse * bandwidth * t.^2));
~~~

 使用spectrogram函数进行STFT，得

![image-20240126051820839](C:\Users\ZJ\AppData\Roaming\Typora\typora-user-images\image-20240126051820839.png)

参数与要求相同。

### 调整STFT窗长使三个信号在时频图上较好分辨

将三个信号值相加并做STFT，经调整窗长至34，重叠长度为窗长的0.9倍，得到的较好结果如下（该问与上一问的STFT幅度均为取对数的结果）

![image-20240126052133666](C:\Users\ZJ\AppData\Roaming\Typora\typora-user-images\image-20240126052133666.png)

+ 取较小窗长，如15

  ![image-20240126052419582](C:\Users\ZJ\AppData\Roaming\Typora\typora-user-images\image-20240126052419582.png)

+ 取较大窗长，如70

  ![image-20240126052543859](C:\Users\ZJ\AppData\Roaming\Typora\typora-user-images\image-20240126052543859.png)

+ 分析：

  + 在固定重叠比例的情况下，winodw_size在30-40既可以比较清晰地分辨wave_data的变化趋势，也在较长的时间里可以分辨两个复指数信号。
  + 窗更小时由于采样点不够，会发生较明显的栅栏效应，虽然提高了时域分辨率，但频域分辨率降低了，使得难以分辨两个复指数信号的频谱。
  + 窗更大时两个复指数信号的宽度更窄，更易区分，频域分辨率增大，但由于窗长增加，不同时间的频率被乘以一个相位再相加，在不同时间的相同频率会发生混叠，使得信号的时频特性无法区分，时域分辨率降低。

### 加噪声的chirp信号参数估计

chirp信号只存在在观测时间的一段时间内，且整个观测时间都存在高斯白噪声，需要估计chirp信号的加入时间、持续时间、起始频率、带宽。

+ 解决该问题分为两步：

  + 估计信号起始时间与持续时间（结束时间）
  + 在持续时间内估计chirp信号的起始频率与带宽

+ **估计信号起始时间与持续时间**

  对观测到的进行STFT后，注意到存在chirp信号对应的时频点幅度较高，以故STFT结果的不同时间的能量值为指标，进行信号的时间划分。

  ~~~matlab
  [S, F, T] = spectrogram(observe_signal, window_size, round(overlap*window_size), nfft, fs);
  energy = sum(abs(S).^2) / nfft;
  diff_energy = diff(energy);    % 能量差分值
  [peaks,peaks_index]=findpeaks(abs(diff_energy));
  threshold = max(5*mean(abs(diff_energy)),0.7*max(peaks));  % 设置差分阈值，用于检测脉冲
  pulse_indices = find(abs(diff_energy) > threshold);        % 检测脉冲出现的时间
  ~~~

  经实验，该结果在较高的信噪比（>15dB）下准确率较高，在信噪比较低时，能量曲线的差分值将会无法分辨信号的位置，此时采用能量值的大小进行估计。

  ~~~matlab
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
  ~~~

+ **估计信号起始频率与带宽**

  确定信号起始时间与持续时间后，将时域信号从观测信号中取出。经查阅文献，采用相位展开（phase unwarapping）的方法进行估计。（参考论文《Parameter estimation of chirp signals》）

  + 原理：

    chirp信号经抽样后形式可变为
    $$
    x_n=A\exp(j2\pi(\frac{1}{2}\alpha n^2+fn))+e_n \;,\;n_0\le n\le N+n_0-1
    $$
    其中 $B=N\alpha f_s,\text{start_freq}=f*f_s$

    在信噪比较高时，有
    $$
    x_n\approx A\exp(j2\pi(\frac{1}{2}\alpha n^2+fn)+w_n) \;,\;n_0\le n\le N+n_0-1
    $$
    其中 $w_n$ 为实高斯白噪声信号

    则有
    $$
    \phi_n=\pi \alpha n^2+2\pi fn+w_n\\
    \Delta^2\phi_n=\Delta\phi_n-\Delta\phi_{n-1}=\phi_n-2\phi_{n-1}+\phi_{n-2}=2\pi\alpha +\Delta^2w_n
    $$
    若能限制 $-0.5<\alpha<0.5$ ，使得$\Delta^2\phi_n$可以以大概率落在 $(-\pi,\pi)$内，则可以通过 $\Delta^2\phi_n$ 得到原来的连续相位，求解得到 $\alpha,f$ 

  + 实现的流程图：
    $$
    y_n=x_nx^*_{n-1}\rightarrow z_n=y_ny^*_{n-1}\rightarrow\Delta^2\phi_n=\angle z_n\rightarrow\phi_n=\Delta^2\phi_n+2\phi_{n-1}-\phi_{n-2}\\
    \phi_{n_0}=\angle x_{n_0},\phi_{n_0+1}=\angle x_{n_0+1}x^*_{n_0}+\phi_{n_0}
    $$
    估计 $\theta=[h \;\;2\pi f\;\;\pi\alpha]=(G^TG)^{-1}G^T\phi$ 

    其中 $h$ 为初相位，本题不关心，$G=[(\text{ones}(1,N))^T,(1:N)^T,((1:N).\text{^}2)^T]$

    再根据 $B=N\alpha f_s,\text{start_freq}=f*f_s$ 即可得到起始频率与带宽的估计

    由于相位展开得到的是连续相位，在代码实现中，需要把带宽与起始频率的估计结果限制在 $f_s$ 内

    ~~~matlab
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
    ~~~

  + 实验

    对每个SNR做50次蒙特卡洛实验，具体实验参数见代码dsp3_1.m，结果如下

    ![image-20240126062202927](C:\Users\ZJ\AppData\Roaming\Typora\typora-user-images\image-20240126062202927.png)

    ![image-20240126062208371](C:\Users\ZJ\AppData\Roaming\Typora\typora-user-images\image-20240126062208371.png)

    ![image-20240126062213969](C:\Users\ZJ\AppData\Roaming\Typora\typora-user-images\image-20240126062213969.png)

    ![image-20240126062219119](C:\Users\ZJ\AppData\Roaming\Typora\typora-user-images\image-20240126062219119.png)

    其中起始时间和持续时间的量级为 $10^{-4}s$ ，均方误差都在信噪比为5dB时都达到了 $10^{-6}s$ 量级，较为准确。起始频率与带宽的量级均为 $10^6\text{Hz}$ ，在信噪比为5dB往后都维持在 $10^3\text{Hz}$ 量级的误差，估计较为准确。

### 加噪声的多chirp信号参数估计

考虑几个不同参数的chirp信号，都在观测时间的某一段时间里存在，且加有高斯白噪声，对这些chirp信号进行参数估计。

+ 同上一问，首先需要确定chirp信号在观测信号中的位置，再对chirp参数进行估计。但多chirp信号的叠加带来了更多的麻烦，上一问的相位展开方法无法估计同一段时间内两个chirp信号叠加的参数，且多个chirp信号的位置排列更多，需要考虑更全面的情况。

+ 为方便与启发性起见，本题固定信号数为2，且估计算法知道信号数。

+ 经过文献调研，对于多chirp叠加信号的参数估计方法主要有极大似然估计、分数阶傅里叶变换（FRFT）等算法及各种改进快速算法，由于时间精力原因，在本题中使用FRFT对叠加信号进行分析，并在时频平面上搜索最优变换阶数以估计chirp信号参数，其中FRFT的函数代码参考 https://github.com/nadavleh/FrFT/blob/master/frft.m

+ **估计信号起始时间与持续时间**

  同上一问，我也打算使用STFT结果计算能量并进行时间分割，但多chirp信号的起伏与排列更多，不易直接实现。经过分析，采用以下方法：

  + 对STFT结果的大部分噪声时频点直接滤除

    ~~~matlab
    k=mean(abs(S));
    S(abs(S)<3*k)=0.01*mean(k);  % 直接滤除大部分噪声
    ~~~

  + 对滤过噪声后的信号计算每一时间窗的能量，取对数以增大噪声与信号值幅度差，减小信号叠加带来的能量突起，再用中值滤波器滤除剩余的噪声脉冲，最后对能量进行差分。信号出现与结束的地方将出现显著的尖端脉冲，并以以下逻辑估计叠加信号的最早的起始时间和总持续时间：

    + 由于信号起始点的差分脉冲必大于0，且经过实验分析，如果不存在时域重叠，两个信号需要相隔200-300个点才能在窗能量谱中分开，故要找到第一个大于0的脉冲满足其下一个脉冲小于0且两者之间相距点数大于300，作为信号的起始点
    + 使用相同的过程反过来找到信号的截止点
    + 若在起始点与截止点之间仍存在超过两个脉冲，则认为该叠加信号的两个chirp信号时域不重叠，找到离起始点和截止点最近的两个符合条件的脉冲，作为两个信号的时间划分。
    + 否则认为存在时域重叠，需要进一步划分两个信号的时间

    ~~~matlab
     diff_energy = diff(energy);    % 能量差分值
     [peaks,peaks_index]=findpeaks(abs(diff_energy));
     threshold = 0.5*max(peaks); % 设置差分阈值，用于检测脉冲
     % 在高信噪比(如大于22dB)时，该阈值可能会过高而出现时间分割错误
     pulse_indices = peaks_index(find( peaks > threshold));
     pulse_value = diff_energy(pulse_indices);
    ~~~

  + 对于时域重叠的情况，先根据信号的起始点和截止点取出有效信号段，再重新做STFT，滤除大部分噪声（事实上由于时间缩小，信号功率增大，基本上所有的噪声都能被滤除）。再对时频图的每个时间找出脉冲的个数，即为信号数，经过中值滤波滤除噪声脉冲后，得到一个服从1，2两点分布的信号，很好地描述了两个chirp信号重叠的位置。

  + 由此得到了存在时域重叠的信号的三或四个分隔点，将叠加信号拆分成重叠段与非重叠段分别估计参数，重点对重叠段进行参数估计。

+ **估计信号起始频率与带宽**

  对于时域重叠段，在 $[t_1,t_2]$ 内，存在两个chirp信号并充满了整个时间。经过文献调研，发现分数阶傅里叶变换可以用于处理这类信号，且有比较好的几何直观。（参考：https://github.com/realman01/FRFT/blob/main/wanghuaijin.pdf）

  + FRFT的原理： 

    ![image-20240126072009052](C:\Users\ZJ\AppData\Roaming\Typora\typora-user-images\image-20240126072009052.png)

    其中 $\alpha = p\pi/2$ 相当于FRFT对时频平面的旋转角。

    ![image-20240126072152805](C:\Users\ZJ\AppData\Roaming\Typora\typora-user-images\image-20240126072152805.png)

    可见在将时频平面的时轴旋转到与chirp信号垂直时，chirp信号表现为冲激信号，而在别的角度会投影在一段区间内。

    故使用FRFT估计叠加chirp参数的过程就是：

    + 搜索旋转角使得在某个旋转角，chirp信号的最大值达到最大

    + 以该旋转角对应的分数阶做FRFT，得到chirp在该阶FRFT的冲击信号

    + 测量其在对应坐标轴上的偏移，代入参数估计公式

      ![image-20240126073043957](C:\Users\ZJ\AppData\Roaming\Typora\typora-user-images\image-20240126073043957.png)

      其中 $s=\sqrt{TF_S}$ 为归一化时频系数。

    在具体实现中，有
    $$
    B_{estimate} = -f_s\cot(\alpha) =- f_s\cot(p\pi/2)\\
    freq_{中心} = \mod(u_{estimate}\csc(p\pi/2),f_s)\\
    freq_{起始} = freq_{中心}-B_{estimate}
    $$


          ~~~matlab
  function [overlap_begin_freq_estimate,... % 也可以用在非时域重叠信号中
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
          [peaks,peaks_index]=findpeaks(abs(frft_signal));
          maxx_index = peaks_index(peaks>0.7*max(peaks));    % 调频斜率相同时，可能有两个或多个peak
          u_estimate=fff(maxx_index);
          f_middle=[f_middle, mod(u_estimate*csc(ssss*pi/2),fs)];
      end
      if length(B_estimate_overlap) < length(f_middle) % 只可能为调频斜率一样，但中心频率不同
          B_estimate_overlap=[B_estimate_overlap(1),B_estimate_overlap(1)];
      end
          overlap_end_freq_estimate = f_middle + 1/2*B_estimate_overlap;
          overlap_begin_freq_estimate = f_middle - 1/2*B_estimate_overlap;
  end
          ~~~

  下图为测试算法时生成，为带噪分辨两个叠加chirp信号，可见FRFT分辨不同chirp信号的效果十分显著

  ![image-20240126074501204](C:\Users\ZJ\AppData\Roaming\Typora\typora-user-images\image-20240126074501204.png)

  ![image-20240126074506144](C:\Users\ZJ\AppData\Roaming\Typora\typora-user-images\image-20240126074506144.png)

  ![image-20240126074510931](C:\Users\ZJ\AppData\Roaming\Typora\typora-user-images\image-20240126074510931.png)

  由于未采用预处理产生精确范围的阶数的方法，每估计一个重叠信号需要做1002次FRFT，带来了极高的时间复杂度。

+ 对于未重叠部分，我也采用了该方法，具体在调用estimate_overlap_chirp(overlap_signal,fs,overlap)时将overlap置零。

+ 得到叠加信号的重叠部分与非重叠部分的信号参数估计后，只需比较不同参数的相似程度，即可将两到三个分割的时间拼接起来，得到完整的信号参数估计（起始时间、持续时间、起始频率、带宽），具体见dsp_3_2.m。

+ 结果：

  由于极高的时间复杂度，对于每个SNR只进行了三次蒙洛卡特实验，且由于无需也无法对估计信号与原信号的序号匹配，计算RMSE时采用两信号参数的和与差的绝对值而非参数本身，实验结果如下，实验参数见dsp_3_2.m

  ![image-20240126075304149](C:\Users\ZJ\AppData\Roaming\Typora\typora-user-images\image-20240126075304149.png)

  ![image-20240126075313015](C:\Users\ZJ\AppData\Roaming\Typora\typora-user-images\image-20240126075313015.png)

  参数的量级与上一问相同，但实验结果与实验前测试结果存在误差。实验中对于持续时间和带宽的估计产生了相较实验前测试较大的误差，原因可能是时间分隔算法在不同的噪声环境下不够鲁棒，进而导致带宽的估计也产生了一定的偏差。

  总体而言估计参数的效果良好，但误差都不随信噪比提高而减小，而是维持在一个比较稳定的值，也说明FRFT的参数估计方法可以在低信噪比下达到一个比较好的估计效果。

