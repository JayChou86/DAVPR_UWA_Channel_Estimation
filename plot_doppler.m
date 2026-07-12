% 定义文件名和对应的深度
filenames = {
    '475doppler_power_spectrum.mat',
    '624doppler_power_spectrum.mat',
    '824doppler_power_spectrum.mat',
    '924doppler_power_spectrum.mat'
    };

depths = [475, 624, 824, 924]; % 对应的深度值

% 预定义不同的线型和颜色
line_styles = {'-', '--', ':', '-.'};
colors = ['#0072BD', '#D95319', '#EDB120', '#7E2F8E'];
colors = [
        0.00, 0.45, 0.74;  % #0072BD 的RGB值
        0.85, 0.33, 0.10;  % #D95319 的RGB值
        0.93, 0.69, 0.13;  % #EDB120 的RGB值
        0.49, 0.18, 0.56   % #7E2F8E 的RGB值
    ];
% 创建图形
figure('Position', [100, 100,  500, 350]);
hold on;

% 用于图例的字符串
legend_labels = cell(length(filenames), 1);

% 循环加载每个文件并绘制
for i = 1:length(filenames)
    % 检查文件是否存在
    if ~exist(filenames{i}, 'file')
        fprintf('警告: 文件 %s 不存在，跳过\n', filenames{i});
        continue;
    end

    % 加载数据
    data = load(filenames{i});

    % 检查必要的变量是否存在
    if ~isfield(data, 'Doppler_axis') || ~isfield(data, 'S_fd')
        fprintf('警告: 文件 %s 中缺少必要的变量，跳过\n', filenames{i});
        continue;
    end

    % 归一化功率谱（使其最大值为0 dB）
    normalized_power = data.S_fd / max(data.S_fd);

    % 绘制多普勒功率谱
    plot(data.Doppler_axis, 10*log10(normalized_power), ...
        'LineStyle', line_styles{i}, 'Color', colors(i,:), 'LineWidth', 2);

    % 创建图例标签
    legend_labels{i} = sprintf('Node %s (%dm)', char('A'+i-1), depths(i));
end

% 图形美化
xlabel('多普勒频率 (Hz)', 'FontSize', 12);
ylabel('归一化功率 (dB)', 'FontSize', 12);
title('Doppler Power Spectra: Per-Node Comparison (UASN Deployment)', 'FontSize', 14);
grid on;

% 添加图例
 legend(legend_labels, 'Location', 'best', 'FontSize', 10);

% 设置坐标轴范围
xlim([-max(abs(xlim)), max(abs(xlim))]);
xlim([-0.1/2,0.1/2]);

% 添加零频率参考线
y_limits = ylim;
plot([0, 0], y_limits, 'k--', 'LineWidth', 0.5, 'HandleVisibility', 'off');

hold off;

% 保存合并图形
% saveas(gcf, 'combined_doppler_spectra.png');
% saveas(gcf, 'combined_doppler_spectra.fig');
% 
% fprintf('合并图形已保存为 combined_doppler_spectra.png 和 .fig\n');
