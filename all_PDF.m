cc
load('时变冲激响应624_0.1_3.mat')
n_time = length(time_sec);
n_delay = length(delay_ms_all);

fs = 16000;
NFFT = 2^nextpow2(size(h_matrix, 2));
H = fft(h_matrix, NFFT, 2) / NFFT; % 对每一行做FFT，得到每个信道的频率响应
H = fftshift(H,1);
H_sig = H(:, 1:NFFT/2+1);
f_fft = fs/2*linspace(0,1,NFFT/2+1);
imagesc(f_fft, time_sec, abs(H_sig));
H_trun = H_sig(:,5374:7668);
% peaks_freq = abs(H_sig(abs(H_sig) > 10e-4));
% plot(abs(H_sig(1,:)))
% 
% plotfft(h_matrix(1,:),fs)
peaks_freq = abs(H_trun);

% figure
% % 归一化
% data_correct = peaks_freq / sqrt(mean(peaks_freq.^2));
% % 绘制 PDF 对比
% x_axis = linspace(0, 1, 200);
% histogram(data_correct, 1000, 'Normalization', 'pdf'); hold on;
% % 拟合 Rician
% pd_corr = fitdist(data_correct, 'Rician');
% plot(x_axis, pdf(pd_corr, x_axis), 'r-', 'LineWidth', 2);
% title('方法A：局部峰值提取 (正确)');
% subtitle(['拟合 K 因子: s=' num2str(pd_corr.s, '%.2f')]);
% legend('样本直方图', 'Rician 拟合');
% grid on;

figure;
imagesc(delay_ms_all, time_sec, abs(h_matrix));
ylabel('时间 (s)');
xlabel('时延 (ms)');
% title(sprintf('滑动步长为%g,窗长度为%g的时变冲激响应', step_duration,window_duration));colorbar;
colorbar;        % 显示颜色条
axis xy;         % 将y轴的原点设置在左下角
colormap('jet'); % 设置颜色图
xlim([-50, 500]) % 根据需要设置时延显示范围

delay_power = mean(h_matrix, 1);
% figure
% plot(delay_ms_all, delay_power)
% xlim([0,500]);

% figure
% plot(abs(h_matrix(50,:)))
% xlim([0,500]);

% h_abs = abs(h_matrix(:,47860:48580));
h_abs = abs(h_matrix);

% h_abs(h_abs < 0.3) = 0;  % 将 matrix 中所有小于 0.3 的元素置零


threshold_val = max(h_abs(:)) * 10^(-20/20); % -20dB 阈值
min_peak_dist = 200; % 最小峰值距离 (防止把同一个sinc波形的旁瓣当成多径)
valid_peaks = []; % 用于存储提取出来的样本


for t = 1:n_time
    % 取出当前时刻的一列数据
    current_cir = h_abs(t, :);
    
    % 使用 findpeaks 寻找局部极大值
    % 'MinPeakHeight': 过滤噪声
    % 'MinPeakDistance': 过滤同一路径的宽度效应
    % [pks, ~] = findpeaks(current_cir, ...
    %                      'MinPeakHeight', threshold_val, ...
    %                      'MinPeakDistance', min_peak_dist);
    
    pks = maxk(current_cir, 2);

    % 将找到的峰值加入统计池
    valid_peaks = [valid_peaks; pks'];
end



figure
% 归一化
data_correct = valid_peaks / sqrt(mean(valid_peaks.^2));
% 绘制 PDF 对比
x_axis = linspace(0, 4, 100);
histogram(data_correct, 30, 'Normalization', 'pdf'); hold on;
% 拟合 Rician
pd_corr = fitdist(data_correct, 'Rician');
plot(x_axis, pdf(pd_corr, x_axis), 'r-', 'LineWidth', 2);
title('方法A：局部峰值提取 (正确)');
subtitle(['拟合 K 因子: s=' num2str(pd_corr.s, '%.2f')]);
legend('样本直方图', 'Rician 拟合');
grid on;

seg_length = (size(h_matrix,2) + 1)/2;
h = h_matrix(:, seg_length:end);     % 去除负时延
figure;imagesc(abs(h));xlim([1,500])
valid_peaks = abs(h(14:66,1)); %(15:56,351)
plotrician(valid_peaks, 10, 't');  % exp Rician gamma weibull beta normal











