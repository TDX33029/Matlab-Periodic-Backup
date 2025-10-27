close all;
%Define
TOAEnable = 1;
AOAEnable = 1;
TDOAEnable = 1;

syms x;
%创建数据变量
data_input_toa = [];
data_input_aoa = [];
data_input_tdoa = [];
data_original = [];
data_output = [];

%读取配置
list = [];
disp('Reading Config');
baseInfo = csvread('basementInfo.csv');
NumberOfBasement = baseInfo(1,1);
basex = csvread('basementInfo.csv',1,0,[1,0,NumberOfBasement,0]);
basey = csvread('basementInfo.csv',1,1,[1,1,NumberOfBasement,1]);
NumberOfBasement = baseInfo(1,1);
if NumberOfBasement == 2 || NumberOfBasement == 1
    disp('Tow Few Basements !');
    return;
end

%导入数据
data_input_toa  = csvread('TOA_RES.csv');
data_input_aoa  = csvread('AOA_RES.csv');
data_input_tdoa = csvread('TDOA_RES.csv');
data_original = csvread('Original_RES.csv');
data1 = data_input_toa;
data2 = data_input_aoa;
data3 = data_input_tdoa;

%旧版
%基站二维列表
for NumberOfCirculate1 = 1:NumberOfBasement
    if NumberOfCirculate1 ~= NumberOfBasement-1 && NumberOfCirculate1 ~= NumberOfBasement
        list = [list;NumberOfCirculate1,NumberOfCirculate1+1,NumberOfCirculate1+2];
    elseif NumberOfCirculate1 ~= NumberOfBasement
        list = [list;NumberOfCirculate1,NumberOfCirculate1+1,1];
    else
        list = [list;NumberOfCirculate1,1,2];
    end
end

% disp(list);
figure(1);
grid on;
hold on;
xlabel('Displacement X (m)');
ylabel('Displacement Y (m)');
ylim([0,2000]);
axis equal;
set(gca,'ylim',[0 200]);
set(gca,'xlim',[-100 100]);

for cir0 = 1:NumberOfBasement
    label = sprintf('RU%d',cir0);
    text(basex(cir0)+2, basey(cir0)+2, label, 'FontSize', 8, 'Color', 'red');
end
scatter(basex,basey,30,'k');


dt = 0.1;
[trajectory, results] = processDataVariableBS(data_input_toa, data_input_aoa, data_input_tdoa, data_original, basex, basey, dt);

fprintf('Done\nEnd Position: (%.2f, %.2f)\n', trajectory(end,1), trajectory(end,2));


function ekf = initializeEKFVariableBS(basex, basey, init_state, dt)

    ekf = struct();

    ekf.bs_positions = [basex(:), basey(:)];
    ekf.num_bs = length(basex);
    
    if length(basex) ~= length(basey)
        error('RU error');
    end
    
    if ekf.num_bs < 3
        warning('RU lack');
    end
    
    ekf.x = init_state(:);
    ekf.dim_x = length(init_state);
    
    ekf.P = eye(ekf.dim_x) * 10;    %开始的不确定度
    %ekf.P = eye(ekf.dim_x) * 1;
    
    ekf.F = [1, 0, dt, 0;
            0, 1, 0, dt;
            0, 0, 1, 0;
            0, 0, 0, 1];

    %0.1
    process_noise = 0.1;
    process_noise_pos = 1;   %位置过程噪声
    process_noise_vel = 1;   %速度过程噪声

    ekf.Q = [dt^4/4 , 0     , dt^3/2, 0     ;
             0      , dt^4/4, 0     , dt^3/2;
             dt^3/2 , 0     , dt^2  , 0     ;
             0      , dt^3/2, 0     , dt^2  ] * process_noise_pos;

    ekf.Q(3,3) = dt^2 * process_noise_vel;
    ekf.Q(4,4) = dt^2 * process_noise_vel;

    %观测噪声协方差矩阵              标准差
    ekf.R_toa = eye(ekf.num_bs)     * 2.0^2;   
    ekf.R_tdoa = eye(ekf.num_bs-1)  * 2.0^2;
    ekf.R_aoa = eye(ekf.num_bs)     * 0.08^2;  

    %老版
    % ekf.R_toa = eye(ekf.num_bs) * 0.1^2;
    % ekf.R_tdoa = eye(ekf.num_bs-1) * 0.05^2;
    % ekf.R_aoa = eye(ekf.num_bs) * 0.001^2;

    
    ekf.R = blkdiag(ekf.R_toa, ekf.R_tdoa, ekf.R_aoa);
    
    ekf.c = 3e8;
    
    fprintf('Done: %d RU\n', ekf.num_bs);
end

function ekf = predictStep(ekf)
    ekf.x = ekf.F * ekf.x;
    ekf.P = ekf.F * ekf.P * ekf.F' + ekf.Q;
end

function [ekf, innovation] = updateStep(ekf, z_toa, z_tdoa, z_aoa)
    gain_limit = 0.8;
    z_meas = [z_toa(:); z_tdoa(:); z_aoa(:)];
    z_pred = computeExpectedMeasurements(ekf);
    H = computeJacobian(ekf);
    
    S = H * ekf.P * H' + ekf.R;
    K = ekf.P * H' / S;

    K = min(K, gain_limit);
    K = max(K, -gain_limit);
    
    innovation = z_meas - z_pred;
    ekf.x = ekf.x + K * innovation;
    ekf.P = (eye(ekf.dim_x) - K * H) * ekf.P;
end

function z_pred = computeExpectedMeasurements(ekf)  

    x = ekf.x(1); y = ekf.x(2);
    bs = ekf.bs_positions;
    
    % TOA
    pred_toa = zeros(ekf.num_bs, 1);
    for i = 1:ekf.num_bs
        dist = sqrt((x - bs(i,1))^2 + (y - bs(i,2))^2);
        pred_toa(i) = dist;
    end
    
    % TDOA
    pred_tdoa = zeros(ekf.num_bs-1, 1);
    dist_ref = sqrt((x - bs(1,1))^2 + (y - bs(1,2))^2);
    for i = 2:ekf.num_bs
        dist_i = sqrt((x - bs(i,1))^2 + (y - bs(i,2))^2);
        pred_tdoa(i-1) = dist_i - dist_ref;
    end
    
    % AOA
    pred_aoa = zeros(ekf.num_bs, 1);
    for i = 1:ekf.num_bs
        angle = atan2(y - bs(i,2), x - bs(i,1));
        pred_aoa(i) = angle;
    end
    
    z_pred = [pred_toa; pred_tdoa; pred_aoa];   %向量合成结果
end

function H = computeJacobian(ekf)
    x = ekf.x(1); y = ekf.x(2);
    bs = ekf.bs_positions;
    
    total_measurements = ekf.num_bs * 3 - 1;
    H = zeros(total_measurements, ekf.dim_x);
    
    % TOA
    for i = 1:ekf.num_bs
        dx = x - bs(i,1);
        dy = y - bs(i,2);
        dist = sqrt(dx^2 + dy^2);
        
        if dist > 1e-6
            H(i, 1) = dx / dist;
            H(i, 2) = dy / dist;
        end
    end
    
    % TDOA
    dist_ref = sqrt((x - bs(1,1))^2 + (y - bs(1,2))^2);
    for i = 2:ekf.num_bs
        idx = ekf.num_bs + (i-1);
        
        dx_i = x - bs(i,1);
        dy_i = y - bs(i,2);
        dist_i = sqrt(dx_i^2 + dy_i^2);
        
        dx_ref = x - bs(1,1);
        dy_ref = y - bs(1,2);
        
        if dist_i > 1e-6 && dist_ref > 1e-6
            H(idx, 1) = dx_i/dist_i - dx_ref/dist_ref;
            H(idx, 2) = dy_i/dist_i - dy_ref/dist_ref;
        end
    end
    
    % AOA
    aoa_start_idx = ekf.num_bs + (ekf.num_bs-1);
    for i = 1:ekf.num_bs
        idx = aoa_start_idx + i;
        
        dx = x - bs(i,1);
        dy = y - bs(i,2);
        dist_sq = dx^2 + dy^2;
        
        if dist_sq > 1e-6
            H(idx, 1) = -dy / dist_sq;
            H(idx, 2) = dx / dist_sq;
        end
    end
end

function pos = getPosition(ekf)
    pos = ekf.x(1:2)';
end

function [estimated_trajectory, results] = processDataVariableBS(data_toa, data_aoa, data_tdoa, data_original, basex, basey, dt)

    %初始位置跳过
    delay_points = 1;
    %
    
    num_points = size(data_toa, 1);
    num_bs = length(basex);   

    data_toa_m = data_toa * 3e8;
    data_tdoa_m = data_tdoa * 3e8;

    
    data_toa_processed = data_toa_m(:, 1:num_bs);
    data_aoa_processed = data_aoa(:, 1:num_bs);
    data_tdoa_processed = data_tdoa_m(:, 1:(num_bs-1));
    
    init_positions = zeros(delay_points, 2);

    for k = 1:delay_points
        init_pos = improvedInitialPosition(data_toa_processed(k,:), data_aoa_processed(k,:), basex, basey);
        init_positions(k, :) = init_pos;
    end

    avg_init_pos = mean(init_positions, 1);
    init_state = [avg_init_pos(1); avg_init_pos(2); 0; 0];
    
    % Initialize EKF
    ekf = initializeEKFVariableBS(basex, basey, init_state, dt);
    
    estimated_trajectory = zeros(num_points, 2);
    total_measurements = num_bs + (num_bs-1) + num_bs;  % TOA + TDOA + AOA
    innovations = zeros(num_points, total_measurements);
    covariance_trace = zeros(num_points, 1);
    
    
    for k = 1:num_points
        if k <= delay_points
            simple_pos = estimateInitialPositionVariableBS(data_toa_processed(k,:), data_aoa_processed(k,:), basex, basey);
            estimated_trajectory(k, :) = simple_pos;

            % 直接继承结果
            ekf.x(1:2) = simple_pos';
            ekf.x(3:4) = [0; 0];
        else
            z_toa = data_toa_processed(k, :)';
            z_aoa = data_aoa_processed(k, :)';
            z_tdoa = data_tdoa_processed(k, :)';
            
            % 正常的EKF
            ekf = predictStep(ekf);
            [ekf, innovation] = updateStepWithGeometryAwareness(ekf, z_toa, z_tdoa, z_aoa);
            %[ekf, innovation] = updateStepWithGeometryAwareness(ekf, z_toa, z_tdoa, z_aoa);
            %[ekf, innovation] = updateStep(ekf, z_toa, z_tdoa, z_aoa);
            
            estimated_trajectory(k, :) = getPosition(ekf);
            innovations(k, :) = innovation';
        end

        covariance_trace(k) = trace(ekf.P);
        
    end
    
    results.estimated_trajectory = estimated_trajectory;
    results.innovations = innovations;
    results.covariance_trace = covariance_trace;
    results.time_vector = (0:num_points-1)' * dt;
    results.num_bs = num_bs;
    results.bs_positions = [basex(:), basey(:)];
    
    disp('Done');
    
    % plot
    plotResultsVariableBS(results,data_original);
end

function plotResultsVariableBS(results,data_original)
    plot(results.estimated_trajectory(:,1), results.estimated_trajectory(:,2), 'r-', 'LineWidth', 2);
    plot(data_original(:,1),data_original(:,2),'b-', 'LineWidth', 2);
    hold on;
    xlabel('Displacement X (m)');
    ylabel('Displacement Y (m)');
    grid on;

end

function init_pos = estimateInitialPositionVariableBS(toa, aoa, basex, basey)
    num_bs = length(basex);
    candidate_positions = zeros(num_bs, 2);
    
    for i = 1:num_bs
        dist = toa(i);
        angle = aoa(i);
        candidate_positions(i, 1) = basex(i) + dist * cos(angle);
        candidate_positions(i, 2) = basey(i) + dist * sin(angle);
    end
    
    init_pos = mean(candidate_positions, 1);
end

function init_pos = improvedInitialPosition(toa, aoa, basex, basey)

    num_bs = length(basex);
    
    %最小二乘法处理初始点
    if num_bs >= 3
        init_pos = leastSquaresPosition(toa, aoa, basex, basey);
        return;

    end
  
    for i = 1:num_bs
        candidate_positions(i, 1) = basex(i) + toa(i) * cos(aoa(i));%简单的三角定位
        candidate_positions(i, 2) = basey(i) + toa(i) * sin(aoa(i));
        
        % 权重
        weights(i) = 1 / (toa(i) + 0.001);  %防止分子无效
    end
    
    weights = weights / sum(weights);  % 归一化
    init_pos = sum(candidate_positions .* weights, 1);
    
end

function pos = leastSquaresPosition(toa, aoa, basex, basey)
   
    num_bs = length(basex);
    
    % 线性系统: A * pos = b
    A = zeros(2*num_bs, 2);
    b = zeros(2*num_bs, 1);
    
    for i = 1:num_bs
        % TOA
        idx = 2*i - 1;
        A(idx, 1) = 2 * basex(i);
        A(idx, 2) = 2 * basey(i);
        b(idx) = basex(i)^2 + basey(i)^2 - toa(i)^2;
        
        % AOA
        idx = 2*i;
        A(idx, 1) = -sin(aoa(i));
        A(idx, 2) = cos(aoa(i));
        b(idx) = -basex(i)*sin(aoa(i)) + basey(i)*cos(aoa(i));
    end
    
    % 最小二乘
    pos = (A' * A) \ (A' * b);
    pos = pos';
end

function ekf = buildAdaptiveNoiseMatrix(ekf, geometry_quality)  %根据基站距离决定AOA权重
    
    if geometry_quality < 0.3
        R_aoa_adapted = ekf.R_aoa * 100;
    else
        R_aoa_adapted = ekf.R_aoa;
    end
    %再次应用
    ekf.R = blkdiag(ekf.R_toa, ekf.R_tdoa, R_aoa_adapted);
end

function [ekf, innovation] = updateStepWithGeometryAwareness(ekf, z_toa, z_tdoa, z_aoa)
    
    % 步骤1: 检测几何质量
    geometry_quality = checkGeometryQuality(ekf, 3.0);
    
    % 步骤2: 根据几何质量调整噪声矩阵 ← 在这里调用！
    ekf = buildAdaptiveNoiseMatrix(ekf, geometry_quality);
    
    % 步骤3: 正常EKF更新
    z_meas = [z_toa(:); z_tdoa(:); z_aoa(:)];
    z_pred = computeExpectedMeasurements(ekf);
    H = computeJacobian(ekf);
    
    innovation = z_meas - z_pred;
    S = H * ekf.P * H' + ekf.R;
    K = ekf.P * H' / S;
    
    ekf.x = ekf.x + K * innovation;
    ekf.P = (eye(ekf.dim_x) - K * H) * ekf.P;
end

function geometry_quality = checkGeometryQuality(ekf, min_safe_distance)
    
    if nargin < 2
        min_safe_distance = 10.0;  % 小于10m加大AOA观测噪声
    end
    
    x = ekf.x(1); y = ekf.x(2);
    bs_positions = ekf.bs_positions;
    
    geometry_quality = 1.0;  % 初始质量系数
    
    min_distance = inf;
    for i = 1:size(bs_positions, 1)
        distance = norm([x, y] - bs_positions(i, :));
        min_distance = min(min_distance, distance);
        
        if distance < min_safe_distance     %应用距离检测
            penalty = exp((min_safe_distance - distance) / 2.0);
            geometry_quality = geometry_quality / penalty;
        end
    end
    %几何精度和角度多样化
    H = computeJacobian(ekf);
    if rank(H) < size(H, 2)
        geometry_quality = geometry_quality * 0.1;  % 效果差亏损，有秩亏出现
    else
        cond_number = cond(H);
        if cond_number > 1e8
            gdop_penalty = 1.0 / log10(cond_number);
            geometry_quality = geometry_quality * gdop_penalty;
        end
    end
    angle_diversity = checkAngleDiversity(ekf);
    geometry_quality = geometry_quality * angle_diversity;
    
    % 限幅
    geometry_quality = max(0.01, min(1.0, geometry_quality));
end

function diversity_score = checkAngleDiversity(ekf)
    
    x = ekf.x(1); y = ekf.x(2);
    bs_positions = ekf.bs_positions;
    num_bs = size(bs_positions, 1);
    
    if num_bs < 3
        diversity_score = 0.5;  % RU Lack 
        return;
    end
    
    % 计算基站角度
    angles = zeros(num_bs, 1);
    for i = 1:num_bs
        dx = bs_positions(i, 1) - x;
        dy = bs_positions(i, 2) - y;
        angles(i) = atan2(dy, dx);
    end
    
    angles_sorted = sort(angles);
    angles_sorted = [angles_sorted; angles_sorted(1) + 2*pi];  % 循环
    
    angle_gaps = diff(angles_sorted);
    max_gap = max(angle_gaps);
    
    ideal_gap = 2 * pi / num_bs;
    
    if max_gap > pi                 % 方向差
        diversity_score = 0.3;
    elseif max_gap > ideal_gap * 2  %一般
        diversity_score = 0.6;
    else
        diversity_score = 0.9;      %好
    end
    
end