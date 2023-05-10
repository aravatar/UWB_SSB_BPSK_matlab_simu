%% 模拟TR-UWB信号发生器和接收器
%SSB调制

% 定义变量
close all;
clear;
Ep = 1;                     %振幅
fs = 3e10;      %Hz   采样率
Fn = fs/2;      %Hz   最小奈奎斯特率
TD = 50e-9;     %seconds     参考脉冲和信息脉冲之间的距离
TS = 300e-9;    %seconds     参考脉冲之间的距离
tp = 2e-9;      %seconds     脉冲宽度
t=-1e-9:1/fs:1.85e-6;       %信号时间
t_ref = t;
t_pulse = t;
a=tp/2.5;       %经验系数

%% 生成单个脉冲
Y=(1-(4*pi.*(t.^2))/a^2) .* exp(-2*pi.*(t.^2)/a^2) / sqrt(Ep);
plot(t,Y,'-r');
axis([-1e-9 1e-9 -0.6 1.2])
ylabel('幅度', 'fontSize',14)
xlabel('时间/s','fontSize',14)
title('单脉冲','fontSize',14)

%% 创建延迟TS的参考脉冲
binary_seq = [1 0 1 1 0 1 0];
for i = 1:length(binary_seq)       
        t_ref = t_ref - TS;
        Y = Y + (1-(4*pi.*(t_ref.^2))/a^2) .* exp(-2*pi.*(t_ref.^2)/a^2) / sqrt(Ep);
end

%% 产生延迟TD和TS的信息脉冲
for i = 1:length(binary_seq)
    if i == 1
        t_pulse = t_pulse - TD;
        Y = Y + (1-(4*pi.*(t_pulse.^2))/a^2) .* exp(-2*pi.*(t_pulse.^2)/a^2) / sqrt(Ep); 
    end
    if i > 1
        if binary_seq(i) == 1
            t_pulse = t_pulse - TS;
            Y = Y + (1-(4*pi.*(t_pulse.^2))/a^2) .* exp(-2*pi.*(t_pulse.^2)/a^2) / sqrt(Ep);
        end
        if binary_seq(i) == 0
            t_pulse = t_pulse - TS;
            Y = Y + (1-(4*pi.*(t_pulse.^2))/a^2) .* exp(-2*pi.*(t_pulse.^2)/a^2) / sqrt(100);
        end
    end
end

%% 生成的信号及其与时间的关系
binary_signal = Y;
figure
subplot(421);
plot(t,binary_signal,'-r')
axis([-40e-9 1.9e-6 -0.6 1.2])
ylabel('幅度', 'fontSize',8)
xlabel('时间/s','fontSize',8)
title('二进制编码脉冲信号图','fontSize',8)

% figure
% plot(t,binary_signal,'-r')
% axis([-3e-9 6e-8 -0.6 1.2])
% ylabel('幅度', 'fontSize',14)
% xlabel('时间/s','fontSize',14)
% title('区域放大：1','fontSize',14)
% 
% figure
% plot(t,binary_signal,'-r')
% axis([2.9e-7 3.6e-7 -0.6 1.2])
% ylabel('幅度', 'fontSize',14)
% xlabel('时间/s','fontSize',14)
% title('区域放大：0','fontSize',14)


%% 生成并绘制频谱

NFFY=2.^(ceil(log(length(binary_signal))/log(2)));
FFTY=fft(binary_signal,NFFY);
NumUniquePts=ceil((NFFY+1)/2); 
FFTY=FFTY(1:NumUniquePts);
MY=abs(FFTY);
MY=MY*2;
MY(1)=MY(1)/2;
MY(length(MY))=MY(length(MY))/2;
MY=MY/length(binary_signal);
f=(0:NumUniquePts-1)*2*Fn/NFFY;
% figure
subplot(422); 
plot(f,20*log10(MY));
xlabel('频率/Hz','fontSize',8); 
title('功率谱图','fontSize',8);
ylabel('功率密度(dB/Hz)','fontSize',8);

%% SSB调制基带信号并绘制
fp =4.492e9;
[modulated_signal,tm] = modulate(binary_signal,fp,fs,'amssb');
% figure
subplot(423); 
plot(tm,modulated_signal,'-r')
axis([-4e-8 1.9e-6 -1.2 1.2])
ylabel('幅度', 'fontSize',8)
xlabel('时间/s','fontSize',8)
title('调制信号图','fontSize',8)

% figure
% plot(tm,modulated_signal,'-r')
% axis([2.95e-7 3.55e-7 -1.2 1.2])
% ylabel('幅度', 'fontSize',14)
% xlabel('时间/s','fontSize',14)
% title('区域放大','fontSize',14)

%% 生成并绘制调制信号的频谱

NFFY=2.^(ceil(log(length(modulated_signal))/log(2)));
FFTY=fft(modulated_signal,NFFY);
NumUniquePts=ceil((NFFY+1)/2); 
FFTY=FFTY(1:NumUniquePts);
MY=abs(FFTY);
MY=MY*2;
MY(1)=MY(1)/2;
MY(length(MY))=MY(length(MY))/2;
MY=MY/length(modulated_signal);
f=(0:NumUniquePts-1)*2*Fn/NFFY;
% figure
subplot(424); 
plot(f,20*log10(MY));
xlabel('频率/Hz','fontSize',8); 
title('调制信号功率谱图','fontSize',8);
ylabel('功率密度(dB/Hz)','fontSize',8);


%% 对接收到的信号进行平方并绘制
modulated_signal=awgn(modulated_signal,30);
received_signal_sqr = modulated_signal.^2;
% figure
subplot(425); 
plot(tm,received_signal_sqr,'-r')
axis([-4e-8 1.9e-6 -0.2 1.2])
ylabel('幅度', 'fontSize',8)
xlabel('时间/s','fontSize',8)
title('接收的含噪信号平方','fontSize',8)
% 
% figure
% plot(tm,received_signal_sqr,'-r')
% axis([2.95e-7 3.55e-7 -0.2 1.2])
% ylabel('幅度', 'fontSize',14)
% xlabel('时间/s','fontSize',14)
% title('区域放大','fontSize',14)

%% 生成并绘制平方接收信号的频谱
NFFY=2.^(ceil(log(length(received_signal_sqr))/log(2)));
FFTY=fft(received_signal_sqr,NFFY);
NumUniquePts=ceil((NFFY+1)/2); 
FFTY=FFTY(1:NumUniquePts);
MY=abs(FFTY);
MY=MY*2;
MY(1)=MY(1)/2;
MY(length(MY))=MY(length(MY))/2;
MY=MY/length(received_signal_sqr);
f=(0:NumUniquePts-1)*2*Fn/NFFY;
% figure
subplot(426); 
plot(f,20*log10(MY));
xlabel('频率/Hz','fontSize',8); 
title('接收含噪信号平方功率谱图','fontSize',8);
ylabel('功率密度(dB/Hz)','fontSize',8);


%% 按100ns的批次对信号进行时间积分并绘制
integrated_signal = zeros(1,length(received_signal_sqr));
start_time = min(tm);
next_time = 0;
while next_time < max(tm)
    next_time = next_time + 100e-9; %每100ns
    index = find((start_time <= tm) & (tm < next_time));
    integrated_signal(index) = cumtrapz(received_signal_sqr(index)); %每100ns积分一次
    start_time = next_time;
end

% figure
subplot(427); 
plot(tm,integrated_signal,'-r')
ylabel('幅度', 'fontSize',8)
xlabel('时间/s','fontSize',8)
title('积分后信号','fontSize',8)

%% 设计阈值决策01
recovered_data = zeros(1,7);
threshold = 15;
start_time = min(tm);
next_time = 0;
i = 1;
result = zeros(1,length(binary_signal));
while next_time < max(tm)
    next_time = next_time + 100e-9; %每100ns
    index = find((start_time <= tm) & (tm < next_time));
    if(max(integrated_signal(index)) > threshold && max(integrated_signal(index) > 0))
        recovered_data(i) = 1;
        result(index) = integrated_signal(index);
    else
        recovered_data(i) = 0;
    end
    next_time = next_time + 200e-9; % 延迟200ns
    start_time = next_time;
    i = i + 1;
end
% figure
subplot(428); 
stem(recovered_data,'-r')
ylabel('幅度', 'fontSize',8)
title('结果','fontSize',8)
axis([0 8 -0.1 1.1])
disp('恢复信号: ')
disp(recovered_data)



