%% =================================================================
%  parameter_config_table.m — 部署配置表导出
%  目的: 自动生成深度感知传感器网络节点的配置表
%  输出: deployment_config.json + deployment_config.csv
%  用于论文中的 Table: Deployment Configuration Table
% =================================================================
cc
addpath('utils');

%% 配置
depths = [475, 624, 824, 924];
depth_labels = {'Node A (475.8m, Upper Thermocline)', ...
                'Node B (624.8m, Within Thermocline)', ...
                'Node C (824.8m, Lower Thermocline)', ...
                'Node D (924.8m, Below Thermocline)'};
data_files = {
    '时变冲激响应475_0.1_3.mat',
    '时变冲激响应624_0.1_3.mat',
    '时变冲激响应824_0.1_3.mat',
    '时变冲激响应924_0.1_3.mat'
    };

fs = 16000;
symbol_rate = 4000;
T_step = 0.1;

%% 为每个深度计算信道统计量 + 自动推导最优参数
fprintf('========== 生成部署配置表 ==========\n\n');

config_table = cell(5, 13);
config_table(1, :) = {'Depth', 'Node', 'tau_rms (ms)', 'tau_mean (ms)', ...
    'T_c (ms)', 'B_c (Hz)', 'f_d_rms (Hz)', ...
    'Optimal ff_taps', 'Optimal fb_taps', 'Optimal forget', ...
    'RLS Mul/Sym', 'LMS Mul/Sym', 'Complexity Ratio'};

for d = 1:4
    fprintf('处理深度 %d...\n', depths(d));
    load(data_files{d});

    % 预处理 (标准流程)
    avg_pdp = mean(abs(h_matrix).^2, 1);
    [~, peak_idx] = max(avg_pdp);
    start_idx = max(1, peak_idx - 100);
    end_idx = min(size(h_matrix,2), peak_idx + 2500);
    h_cut1 = h_matrix(:, start_idx:end_idx);

    [m, n] = size(h_cut1);
    result_matrix = zeros(m, n);
    for row = 1:m
        row_data = h_cut1(row, 1:min(150, n));
        [~, max_idx] = max(row_data);
        shift_amount = 100 - max_idx;
        if shift_amount > 0
            result_matrix(row, shift_amount+1:end) = h_cut1(row, 1:end-shift_amount);
        elseif shift_amount < 0
            shift_amount = abs(shift_amount);
            result_matrix(row, 1:end-shift_amount) = h_cut1(row, shift_amount+1:end);
        else
            result_matrix(row, :) = h_cut1(row, :);
        end
    end

    % 自动多径提取
    [h_cut, keep_indices] = extract_multipath_by_threshold(abs(result_matrix), -20, 5);

    % 自动参数选择
    params = auto_parameter_selector(h_cut, fs, symbol_rate, T_step, depth_labels{d});
    stats = params.stats;

    % 填充表格
    config_table{d+1, 1}  = depths(d);
    config_table{d+1, 2}  = depth_labels{d};
    config_table{d+1, 3}  = round(stats.tau_rms, 2);
    config_table{d+1, 4}  = round(stats.tau_mean, 2);
    config_table{d+1, 5}  = round(stats.T_c * 1e3, 1);
    config_table{d+1, 6}  = round(stats.B_c, 1);
    config_table{d+1, 7}  = round(stats.f_d_rms, 4);
    config_table{d+1, 8}  = params.eq_ff_taps;
    config_table{d+1, 9}  = params.eq_fb_taps;
    config_table{d+1, 10} = params.rls_forget;
    config_table{d+1, 11} = round(params.rls_flops_per_sym);
    config_table{d+1, 12} = round(params.lms_flops_per_sym);
    config_table{d+1, 13} = round(100 * params.lms_flops_per_sym / params.rls_flops_per_sym, 1);

    fprintf('  完成: fb=%d, forget=%.4f, %d径\n', ...
        params.eq_fb_taps, params.rls_forget, length(keep_indices));
end

%% 输出到控制台
fprintf('\n========== 部署配置表 ==========\n');
fprintf('%-6s | %-8s | %-8s | %-6s | %-6s | %-6s | %-5s | %-5s | %-6s | %-8s | %-8s\n', ...
    'Depth', 'tau_rms', 'tau_mean', 'T_c', 'B_c', 'f_d', ...
    'ff', 'fb', 'lambda', 'RLS ops', 'LMS ops');
fprintf('%s\n', repmat('-', 1, 95));
for d = 2:5
    fprintf('%-6d | %-8.2f | %-8.2f | %-6.1f | %-6.1f | %-6.4f | %-5d | %-5d | %-6.4f | %-8d | %-8d\n', ...
        config_table{d, 1}, config_table{d, 3}, config_table{d, 4}, ...
        config_table{d, 5}, config_table{d, 6}, config_table{d, 7}, ...
        config_table{d, 8}, config_table{d, 9}, config_table{d, 10}, ...
        config_table{d, 11}, config_table{d, 12});
end

%% 导出 JSON (部署用)
deployment_config = struct();
deployment_config.description = 'Depth-Aware Receiver Pre-Configuration for UASN Nodes';
deployment_config.generated = datestr(now);
deployment_config.parameters.sampling_rate_hz = fs;
deployment_config.parameters.symbol_rate_sps = symbol_rate;
deployment_config.parameters.modulation = 'QPSK';
deployment_config.parameters.equalizer_type = 'RLS-DFE';

nodes = [];
for d = 1:4
    node = struct();
    node.depth_m = config_table{d+1, 1};
    node.label = config_table{d+1, 2};
    node.channel_stats = struct(...
        'tau_rms_ms', config_table{d+1, 3}, ...
        'tau_mean_ms', config_table{d+1, 4}, ...
        'coherence_time_ms', config_table{d+1, 5}, ...
        'coherence_bandwidth_hz', config_table{d+1, 6}, ...
        'doppler_rms_hz', config_table{d+1, 7});
    node.receiver_config = struct(...
        'forward_taps', config_table{d+1, 8}, ...
        'feedback_taps', config_table{d+1, 9}, ...
        'forgetting_factor', config_table{d+1, 10});
    node.complexity = struct(...
        'rls_mult_per_symbol', config_table{d+1, 11}, ...
        'lms_mult_per_symbol', config_table{d+1, 12}, ...
        'lms_vs_rls_percent', config_table{d+1, 13});
    nodes = [nodes; node];
end
deployment_config.nodes = nodes;

json_str = jsonencode(deployment_config, 'PrettyPrint', true);
fid = fopen('deployment_config.json', 'w');
fprintf(fid, '%s', json_str);
fclose(fid);
fprintf('\nJSON 配置已导出到 deployment_config.json\n');

%% 导出 CSV
csv_table = cell2table(config_table(2:end, :), ...
    'VariableNames', config_table(1, :));
writetable(csv_table, 'deployment_config.csv');
fprintf('CSV 配置已导出到 deployment_config.csv\n');

%% 生成 LaTeX 表格 (可直接插入论文)
fprintf('\n========== LaTeX 表格代码 ==========\n');
fprintf('%% 论文 Table: Deployment Configuration for UASN Nodes\n');
fprintf('\\begin{table}[htbp]\n');
fprintf('  \\centering\n');
fprintf('  \\caption{Depth-Aware Receiver Pre-Configuration for Fixed Sensor Nodes}\n');
fprintf('  \\label{tab:deployment}\n');
fprintf('  \\begin{tabular}{ccccccccc}\n');
fprintf('    \\hline\n');
fprintf('    Depth (m) & $\\tau_\\mathrm{rms}$ (ms) & $T_c$ (ms) & $B_c$ (Hz) & $f_d^\\mathrm{rms}$ (Hz) & $N_\\mathrm{fb}$ & $\\lambda$ & RLS ops/sym & LMS ops/sym \\\\\n');
fprintf('    \\hline\n');
for d = 2:5
    fprintf('    %d & %.2f & %.1f & %.1f & %.4f & %d & %.4f & %d & %d \\\\\n', ...
        config_table{d, 1}, config_table{d, 3}, config_table{d, 5}, ...
        config_table{d, 6}, config_table{d, 7}, config_table{d, 9}, ...
        config_table{d, 10}, config_table{d, 11}, config_table{d, 12});
end
fprintf('    \\hline\n');
fprintf('  \\end{tabular}\n');
fprintf('\\end{table}\n');

save('deployment_config.mat', 'config_table', 'deployment_config', 'nodes');
disp('全部配置已保存。');
