function plotrician(valid_peaks, num, type)
figure
% 归一化
data_correct = valid_peaks / sqrt(mean(valid_peaks.^2));
% data_correct = valid_peaks / max(valid_peaks);

% 绘制 PDF 对比
x_axis = linspace(0, 1.5, 100);
histogram(data_correct, num, 'Normalization', 'pdf'); hold on;
% 拟合 Rician
pd_corr = fitdist(data_correct, type);
plot(x_axis, pdf(pd_corr, x_axis), 'r-', 'LineWidth', 2);
title('方法A：局部峰值提取 (正确)');
subtitle(['拟合 K 因子: s=' num2str(pd_corr.s, '%.2f')]);
legend('样本直方图', 'Rician 拟合');
grid on;
end