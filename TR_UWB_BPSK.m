%% 模拟TR-UWB信号发生器和接收器
%BPSK调制

% 定义变量
close all;
clear;
Ep = 1;                     %振幅
fs = 3e10;      %Hz   采样率
Fn = fs/2;      %Hz   最小奈奎斯特率
TD = 50e-9;     %seconds     参考脉冲和信息脉冲之间的距离
TS = 300e-9;    %seconds     参考脉冲之间的距离
tp = 2e-9;      %seconds     脉冲宽度
t=0:1/fs:2e-6-1/fs;       %信号时间
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
            Y = Y + (1-(4*pi.*(t_pulse.^2))/a^2) .* exp(-2*pi.*(t_pulse.^2)/a^2) / (-1*sqrt(Ep));
        end
    end
end

%% BPSK调制
binary_signal = Y;
figure
subplot(221);
plot(t,binary_signal,'-r')
ylabel('幅度', 'fontSize',8)
xlabel('时间/s','fontSize',8)
title('BPSK调制脉冲信号图','fontSize',8)

% figure
% plot(t,binary_signal,'-r')
% ylabel('幅度', 'fontSize',14)
% xlabel('时间/s','fontSize',14)
% title('区域放大：1','fontSize',14)
% 
% figure
% plot(t,binary_signal,'-r')
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
subplot(222); 
plot(f,20*log10(MY));
xlabel('频率/Hz','fontSize',8); 
title('功率谱图','fontSize',8);
ylabel('功率密度(dB/Hz)','fontSize',8);

%% 解调
received_signal = awgn(binary_signal,25);
shift_signal=[received_signal(1,1501:60000),zeros(1,1500)];
demot_signal=received_signal.*shift_signal;
% figure
subplot(223); 
plot(demot_signal,'-r')
ylabel('幅度', 'fontSize',8)
xlabel('时间/s','fontSize',8)
title('解调信号','fontSize',8)
% 
% figure
% plot(tm,received_signal_sqr,'-r')
% ylabel('幅度', 'fontSize',14)
% xlabel('时间/s','fontSize',14)
% title('区域放大','fontSize',14)

%% 按100ns的批次对信号进行时间积分并绘制
integrated_signal = zeros(1,length(demot_signal));
start_time = min(t);
next_time = 0;
while next_time < max(t)
    next_time = next_time + 100e-9; %每100ns
    index = find((start_time <= t) & (t < next_time));
    integrated_signal(index) = cumtrapz(demot_signal(index)); %每100ns积分一次
    start_time = next_time;
end

% figure
subplot(224); 
plot(t,integrated_signal,'-r')
ylabel('幅度', 'fontSize',8)
xlabel('时间/s','fontSize',8)
title('积分后信号','fontSize',8)
%% 设计阈值决策01
% recovered_data = zeros(1,7);
% n=1;
% threshold = 0.3;
% i=1;
% while(i<=length(demot_signal))
%     if (demot_signal(i)>0 && abs(demot_signal(i))>threshold)
%         recovered_data(n)=1;
%     elseif (demot_signal(i)<0 && abs(demot_signal(i))<threshold)
%         recovered_data(n)=0;
%     end
%     n=n+1;
%     i=i+9000;
% end
% % figure
% subplot(224); 
% stem(recovered_data,'-r')
% ylabel('幅度', 'fontSize',8)
% title('结果','fontSize',8)
% 
% disp('恢复信号: ')
% disp(recovered_data)



