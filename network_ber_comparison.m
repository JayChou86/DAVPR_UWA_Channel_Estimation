%% =================================================================
%  network_ber_comparison.m — 4节点网络级BER汇总对比
%  目的: 在一张图上展示4个传感器节点的BER性能
%  回答审稿人: 展示824.8m为网络瓶颈 + 深度感知部署的有效性
% =================================================================
cc

%% 配置
depths = [475, 624, 824, 924];
depth_labels = {'Node A (475.8m)', 'Node B (624.8m)', ...
                'Node C (824.8m)', 'Node D (924.8m)'};
colors = [
    0.00, 0.45, 0.74;  % 蓝
    0.85, 0.33, 0.10;  % 橙
    0.93, 0.69, 0.13;  % 黄
    0.49, 0.18, 0.56;  % 紫
    ];
markers = {'o', 's', '^', 'd'};
SNR_vec = 0:5:30;

%% 尝试加载各深度的BER数据
% 优先加载新格式 (_fb*_lambda*), 回退到旧格式 (_%g_%g)
data_dir = '..\五月份实测数据\误码率数据0到30';
all_BER_RLS = [];
all_BER_LMS = [];
loaded_depths = [];
loaded_params = cell(4, 1);

fprintf('加载 BER 数据...\n');
for d = 1:4
    depth = depths(d);
    % 搜索匹配的 .mat 文件
    pattern = sprintf('%d_*.mat', depth);
    files = dir(fullfile(data_dir, pattern));

    found = false;
    for f_idx = 1:length(files)
        fname = files(f_idx).name;
        % 跳过非BER结果文件 (如旧格式的 475_%g_%g.mat 不含 params)
        data = load(fullfile(data_dir, fname));
        if isfield(data, 'BER_vec_RLS') && isfield(data, 'params')
            all_BER_RLS = [all_BER_RLS; data.BER_vec_RLS(:)'];
            if isfield(data, 'BER_vec_LMS')
                all_BER_LMS = [all_BER_LMS; data.BER_vec_LMS(:)'];
            end
            loaded_depths = [loaded_depths, depth];
            loaded_params{d} = data.params;
            fprintf('  加载: %s (RLS BER range: %.1e - %.1e)\n', ...
                fname, min(data.BER_vec_RLS), max(data.BER_vec_RLS));
            found = true;
            break;
        elseif isfield(data, 'BER_vec') && ~found
            % 旧格式: 仅含 BER_vec
            all_BER_RLS = [all_BER_RLS; data.BER_vec(:)'];
            loaded_depths = [loaded_depths, depth];
            fprintf('  加载(旧格式): %s\n', fname);
            found = true;
            break;
        end
    end
    if ~found
        fprintf('  深度 %d: 未找到BER数据\n', depth);
    end
end

%% 绘图1: 网络级BER汇总 (RLS-DFE)
if ~isempty(all_BER_RLS)
    figure('Position', [100, 100, 750, 550]);

    % 标记824.8m为瓶颈
    bottleneck_idx = find(loaded_depths == 824);

    for d = 1:size(all_BER_RLS, 1)
        lw = 1.5;
        ms = 8;
        style = '-';
        if loaded_depths(d) == 824
            lw = 2.5;  % 加粗瓶颈曲线
            ms = 10;
            style = '-';
        end
        h(d) = semilogy(SNR_vec, all_BER_RLS(d, :), [style markers{d}], ...
            'Color', colors(d, :), 'LineWidth', lw, 'MarkerSize', ms);
        hold on;
    end

    % 标注瓶颈
    if ~isempty(bottleneck_idx)
        % 在824曲线旁加文字
        x_pos = 15;
        y_pos = all_BER_RLS(bottleneck_idx, 4);  % SNR=15dB
        text(x_pos, y_pos * 1.5, '← Network Bottleneck', ...
            'Color', colors(bottleneck_idx, :), 'FontSize', 11, 'FontWeight', 'bold');
    end

    grid on;
    xlabel('SNR (dB)', 'FontSize', 13);
    ylabel('Bit Error Rate (BER)', 'FontSize', 13);
    title('Network-Level BER Performance: RLS-DFE Per-Node Evaluation', 'FontSize', 14);

    % 图例
    leg_str = cell(size(all_BER_RLS, 1), 1);
    for d = 1:size(all_BER_RLS, 1)
        if ~isempty(loaded_params{d})
            leg_str{d} = sprintf('%s (fb=%d)', depth_labels{d}, loaded_params{d}.eq_fb_taps);
        else
            leg_str{d} = depth_labels{d};
        end
    end
    legend(h, leg_str, 'Location', 'SouthWest', 'FontSize', 10);
    axis([0 30 1e-5 1]);

    % 保存
    saveas(gcf, 'network_ber_summary.png');
    saveas(gcf, 'network_ber_summary.fig');
    fprintf('网络级BER汇总图已保存。\n');
end

%% 绘图2: RLS vs LMS 对比 (如有数据)
if ~isempty(all_BER_LMS) && size(all_BER_LMS, 1) == size(all_BER_RLS, 1)
    figure('Position', [100, 100, 1200, 500]);

    for d = 1:size(all_BER_RLS, 1)
        subplot(1, size(all_BER_RLS, 1), d);
        semilogy(SNR_vec, all_BER_RLS(d, :), '-o', 'Color', [0.2 0.4 0.8], ...
            'LineWidth', 1.5, 'MarkerSize', 6); hold on;
        semilogy(SNR_vec, all_BER_LMS(d, :), '-s', 'Color', [0.8 0.4 0.2], ...
            'LineWidth', 1.5, 'MarkerSize', 6);
        grid on;
        xlabel('SNR (dB)'); ylabel('BER');
        title(depth_labels{d}, 'FontSize', 10);
        legend('RLS-DFE', 'LMS-DFE', 'Location', 'SouthWest', 'FontSize', 8);
        axis([0 30 1e-5 1]);
    end
    sgtitle('RLS-DFE vs LMS-DFE: Per-Node Performance Comparison', 'FontSize', 14);
    saveas(gcf, 'network_rls_vs_lms.png');
    saveas(gcf, 'network_rls_vs_lms.fig');
end

%% 输出网络级统计
if ~isempty(all_BER_RLS)
    fprintf('\n========== 网络级统计 ==========\n');
    % 找出每个SNR的最差BER (网络瓶颈)
    worst_ber_per_snr = max(all_BER_RLS, [], 1);
    fprintf('SNR(dB) | 最差BER (瓶颈节点)\n');
    for i = 1:length(SNR_vec)
        [~, worst_node] = max(all_BER_RLS(:, i));
        fprintf('  %2d     | %.2e  (%s)\n', ...
            SNR_vec(i), worst_ber_per_snr(i), depth_labels{worst_node});
    end
end

%% 保存汇总数据
if ~isempty(all_BER_RLS)
    save('network_ber_summary.mat', 'all_BER_RLS', 'all_BER_LMS', ...
        'SNR_vec', 'loaded_depths', 'depth_labels', 'loaded_params');
    disp('网络BER汇总数据已保存到 network_ber_summary.mat');
end
