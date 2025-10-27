close all;
%Define
TOAEnable = 1;
AOAEnable = 1;
TDOAEnable = 1;

syms x;
%Create System Value
data_input_toa = [];
data_input_aoa = [];
data_input_tdoa = [];
data_output = [];

%Get Config
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

%Data input
data_input_toa  = csvread('TOA_RES.csv');
data_input_aoa  = csvread('AOA_RES.csv');
data_input_tdoa = csvread('TDOA_RES.csv');
data1 = data_input_toa;
data2 = data_input_aoa;
data3 = data_input_tdoa;
%

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
[trajectory, results] = processDataVariableBS(data_input_toa, data_input_aoa, data_input_tdoa, basex, basey, dt);

fprintf('Done ! Trcking Point Numbers: %d\n', size(trajectory,1));
fprintf('Finally Position: (%.2f, %.2f)\n', trajectory(end,1), trajectory(end,2));


function ekf = initializeEKFVariableBS(basex, basey, init_state, dt)

    ekf = struct();

    ekf.bs_positions = [basex(:), basey(:)];
    ekf.num_bs = length(basex);
    
    if length(basex) ~= length(basey)
        error('basex?basey??????');
    end
    
    if ekf.num_bs < 3
        warning('??????3???????????');
    end
    
    ekf.x = init_state(:);
    ekf.dim_x = length(init_state);
    
    ekf.P = eye(ekf.dim_x) * 10;
    
    ekf.F = [1, 0, dt, 0;
            0, 1, 0, dt;
            0, 0, 1, 0;
            0, 0, 0, 1];
    
    process_noise = 0.1;
    ekf.Q = [dt^4/4, 0, dt^3/2, 0;
            0, dt^4/4, 0, dt^3/2;
            dt^3/2, 0, dt^2, 0;
            0, dt^3/2, 0, dt^2] * process_noise;

    ekf.R_toa = eye(ekf.num_bs) * 2.0^2;           % TOA??
    ekf.R_tdoa = eye(ekf.num_bs-1) * 1.0^2;        % TDOA??
    ekf.R_aoa = eye(ekf.num_bs) * 0.05^2;          % AOA??
    
    ekf.R = blkdiag(ekf.R_toa, ekf.R_tdoa, ekf.R_aoa);
    
    ekf.c = 3e8;
    
    fprintf('?????: %d???\n', ekf.num_bs);
end

function ekf = predictStep(ekf)
    % ????
    ekf.x = ekf.F * ekf.x;
    ekf.P = ekf.F * ekf.P * ekf.F' + ekf.Q;
end

function [ekf, innovation] = updateStep(ekf, z_toa, z_tdoa, z_aoa)
    % ????
    z_meas = [z_toa(:); z_tdoa(:); z_aoa(:)];
    z_pred = computeExpectedMeasurements(ekf);
    H = computeJacobian(ekf);
    
    S = H * ekf.P * H' + ekf.R;
    K = ekf.P * H' / S;
    
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
    
    z_pred = [pred_toa; pred_tdoa; pred_aoa];
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

function [estimated_trajectory, results] = processDataVariableBS(data_toa, data_aoa, data_tdoa, basex, basey, dt)
    % ?????????????????
    % ??:
    %   data_toa: TOA??, N×M??, M?????
    %   data_aoa: AOA??, N×M??
    %   data_tdoa: TDOA??, N×(M-1)??
    %   basex: ??x????
    %   basey: ??y????
    %   dt: ????
    
    % ??????????
    num_points = size(data_toa, 1);
    num_bs = length(basex);
    
    fprintf('=== ???? ===\n');
    fprintf('????: %d\n', num_points);
    fprintf('????: %d\n', num_bs);
    fprintf('TOA????: %dx%d\n', size(data_toa));
    fprintf('AOA????: %dx%d\n', size(data_aoa));
    fprintf('TDOA????: %dx%d\n', size(data_tdoa));
    
    % ??????????????
    if size(data_toa, 2) < num_bs
        error('TOA??????????');
    end
    if size(data_aoa, 2) < num_bs
        error('AOA??????????');
    end
    if size(data_tdoa, 2) < (num_bs-1)
        error('TDOA??????(????-1)');
    end
    
    % ??????
    if max(data_toa(:)) < 0.01
        data_toa_m = data_toa * 3e8;  % ?????
        data_tdoa_m = data_tdoa * 3e8;
        fprintf('??????????????\n');
    else
        data_toa_m = data_toa;
        data_tdoa_m = data_tdoa;
        fprintf('??????\n');
    end
    
    % ???????????
    data_toa_processed = data_toa_m(:, 1:num_bs);
    data_aoa_processed = data_aoa(:, 1:num_bs);
    data_tdoa_processed = data_tdoa_m(:, 1:(num_bs-1));
    
    % ??????
    init_pos = estimateInitialPositionVariableBS(data_toa_processed(1,:), data_aoa_processed(1,:), basex, basey);
    init_state = [init_pos(1); init_pos(2); 0; 0];
    
    fprintf('??????: (%.2f, %.2f)\n', init_pos(1), init_pos(2));
    
    % ???EKF
    ekf = initializeEKFVariableBS(basex, basey, init_state, dt);
    
    % ????
    estimated_trajectory = zeros(num_points, 2);
    total_measurements = num_bs + (num_bs-1) + num_bs;  % TOA + TDOA + AOA
    innovations = zeros(num_points, total_measurements);
    covariance_trace = zeros(num_points, 1);
    
    fprintf('??????...\n');
    
    for k = 1:num_points
        % ??????????
        z_toa = data_toa_processed(k, :)';
        z_aoa = data_aoa_processed(k, :)';
        z_tdoa = data_tdoa_processed(k, :)';
        
        % EKF??
        ekf = predictStep(ekf);
        [ekf, innovation] = updateStep(ekf, z_toa, z_tdoa, z_aoa);
        
        % ????
        estimated_trajectory(k, :) = getPosition(ekf);
        innovations(k, :) = innovation';
        covariance_trace(k) = trace(ekf.P);
        
        % ????
        if mod(k, 50) == 0
            fprintf('??? %d/%d ??\n', k, num_points);
        end
    end
    
    % ????
    results.estimated_trajectory = estimated_trajectory;
    results.innovations = innovations;
    results.covariance_trace = covariance_trace;
    results.time_vector = (0:num_points-1)' * dt;
    results.num_bs = num_bs;
    results.bs_positions = [basex(:), basey(:)];
    
    fprintf('???????\n');
    
    % ??
    plotResultsVariableBS(results);
end

function plotResultsVariableBS(results)
    plot(results.estimated_trajectory(:,1), results.estimated_trajectory(:,2), 'b-', 'LineWidth', 2);
    hold on;
    
%     % ??????
%     bs_positions = results.bs_positions;
%     plot(bs_positions(:,1), bs_positions(:,2), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
%     
    % ??????
%     for i = 1:size(bs_positions,1)
%         text(bs_positions(i,1), bs_positions(i,2), sprintf('BS%d', i), ...
%              'FontSize', 10, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
%     end
    
%     plot(results.estimated_trajectory(1,1), results.estimated_trajectory(1,2), 'go', ...
%          'MarkerSize', 8, 'MarkerFaceColor', 'g');
%     plot(results.estimated_trajectory(end,1), results.estimated_trajectory(end,2), 'ms', ...
%          'MarkerSize', 8, 'MarkerFaceColor', 'm');
    
    xlabel('Displacement X (m)');
    ylabel('Displacement Y (m)');
    %title('Tracking');
    %legend('RU', 'Estimate Tracking', 'Start Point', 'End Point', 'Location', 'best');
    %legend('RU', 'Estimate Tracking', 'Location', 'best');
    grid on;
    
    figure(2);
    plot(results.time_vector, results.covariance_trace, 'r-', 'LineWidth', 1.5);
    xlabel('Time (s)');
    ylabel('????');
    title('??????');
    grid on;
    
    figure(3);
    colors = lines(size(results.innovations,2));
    for i = 1:size(results.innovations,2)
        plot(results.time_vector, results.innovations(:,i), 'Color', colors(i,:), ...
             'DisplayName', sprintf('??%d', i));
        hold on;
    end
    xlabel('Time (s)');
    ylabel('???');
    title(sprintf('???? (%d????)', size(results.innovations,2)));
    legend('show', 'Location', 'best');
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
    
    % ??: ?????????????????
    % weights = 1 ./ (toa .^ 2);  % ????????
    % weights = weights / sum(weights);
    % init_pos = sum(candidate_positions .* weights', 1);
end
