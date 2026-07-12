% cc

% 定义文件名和对应的深度
filenames = {
    '475time_cor.mat',
    '624time_cor.mat',
    '824time_cor.mat',
    '924time_cor.mat'
};
depths = [475, 624, 824, 924]; % 对应的深度值

% 预定义不同的线型和颜色
line_styles = {'-', '--', ':', '-.'};
colors = [
    0.00, 0.45, 0.74;  % #0072BD (蓝色)
    0.85, 0.33, 0.10;  % #D95319 (橙色)
    0.93, 0.69, 0.13;  % #EDB120 (黄色)
    0.49, 0.18, 0.56   % #7E2F8E (紫色)
];

% --- 定义平滑样条拟合参数 ---
smoothing_param = 0.999; 

% 创建图形
figure('Position', [100, 100, 500, 350]);
hold on;

% 用于图例的字符串
legend_labels = cell(length(filenames), 1);

% 循环加载每个文件并进行拟合与绘制
for i = 1:length(filenames)
    % 检查文件是否存在
    if ~exist(filenames{i}, 'file')
        fprintf('警告: 文件 %s 不存在，跳过\n', filenames{i});
        continue;
    end

    % 加载数据
    data = load(filenames{i});

    % 检查必要的变量是否存在
    if ~isfield(data, 'tau_axis_slow') || ~isfield(data, 'rho_tau')
        fprintf('警告: 文件 %s 中缺少必要的变量，跳过\n', filenames{i});
        continue;
    end

    % 归一化功率谱（使其最大值为1）
    normalized_power = data.rho_tau / max(data.rho_tau);

    % 扩展时间轴到负方向
    extended_tau_axis = [-fliplr(data.tau_axis_slow(2:end)), data.tau_axis_slow];
    extended_power = [flip(normalized_power(2:end))', normalized_power'];

    % 使用平滑样条进行曲线拟合
    [fit_curve, gof] = fit(extended_tau_axis(:), extended_power(:), 'smoothingspline', ...
                           'SmoothingParam', smoothing_param);

    % --- 核心修改：分步绘制拟合曲线 ---
    % 1. 使用拟合对象 fit_curve 计算出在原始x轴位置上的平滑y值
    smoothed_y_values = fit_curve(extended_tau_axis);

    % 2. 使用标准的 plot(x, y, ...) 函数进行绘制，这样就可以自由指定属性了
    plot(extended_tau_axis, smoothed_y_values, ...
        'LineStyle', line_styles{i}, 'Color', colors(i,:), 'LineWidth', 1.5);
    
    % 创建图例标签
    legend_labels{i} = sprintf('Node %s (%dm)', char('A'+i-1), depths(i));
end

% 图形美化
xlabel('时间差 (s)', 'FontSize', 12);
ylabel('归一化幅度', 'FontSize', 12);
title('时间间隔相关函数 (平滑样条拟合)', 'FontSize', 14);
grid on;

% 添加图例（这种方式下，图例处理更简单）
legend(legend_labels, 'Location', 'best', 'FontSize', 10);

% 设置坐标轴范围
if ~isempty(extended_tau_axis)
    xlim([-max(abs(extended_tau_axis)), max(abs(extended_tau_axis))]);
end

% 添加零频率参考线
y_limits = ylim;
plot([0, 0], y_limits, 'k--', 'LineWidth', 0.5, 'HandleVisibility', 'off');

% 添加水平参考线
yline(exp(-1), '--k', 'e^{-1}', 'HandleVisibility', 'off', 'LabelVerticalAlignment', 'bottom');
yline(0.5, '--k', '0.5', 'HandleVisibility', 'off', 'LabelVerticalAlignment', 'bottom');
yline(0.9, '--k', '0.9', 'HandleVisibility', 'off', 'LabelVerticalAlignment', 'bottom');

hold off;
