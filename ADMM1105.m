clear; clc;

%% 参数与初始化
max_iter = 50;
tol = 1e-4;          % 收敛阈值
nbeta = 36;            % 故障网中松弛开关数量示例
T = 2;               % 时间步长
rho_s = 5000;
rho_t = 3000;
eta_beta = 0.3;        % 残差/变量平滑因子（指数移动平均），越小平滑越强
lambda_beta = 200;     % 使 beta 向 0/1 集中之凸二次
Ess_Num1=5;Ess_Num=3;

%% 初始化 松弛二进制变量及其平滑版本
beta_prev = 0.5*ones(nbeta,1);    % 上一轮 beta（平滑用）
alphaIEEE_prev = 0.5;                 % 上一轮 alpha（示例单个）
beta_smoothed = beta_prev;
alpha_smoothed = alphaIEEE_prev;

% 全局-边界共识变量（24节点出口/21节点入口）
zP = zeros(1, T); zQ = zeros(1, T); zU = zeros(1, T);
% 局部子网边界变量
local_P1 = zeros(1, T); local_Q1 = zeros(1, T); local_U1 = ones(1, T);  % 正常网
local_P2 = zeros(1, T); local_Q2 = zeros(1, T); local_U2 = ones(1, T);  % 故障网

% ADMM乘子初始化
uP = zeros(1, T); uQ = zeros(1, T); uU = zeros(1, T);
lamdaP = zeros(1, T); lamdaQ = zeros(1, T); lamdaU = zeros(1, T);
% 记录残差与目标
r_hist = []; s_hist = []; obj_hist = [];
beta2_hist = [];
fprintf('ADMM主循环（边界联络线共识协调）\n');

%% SOC相关一致性变量及对偶
% z_soc = 0.6 * ones(max(Ess_Num1,Ess_Num), T+1);     % 全局SOC一致性变量
lam1 = zeros(Ess_Num1, T+1);   % 正常网对偶
lam2 = zeros(Ess_Num, T+1);   % 故障网对偶
% 保存收敛
% prim_soc_hist = [];
% dual_soc_hist = [];

for k = 1:max_iter
    fprintf('---------- ADMM 第 %d 轮 ----------\n', k);
    %% 子网一：正常网（24节点出口边界），带惩罚项与乘子变量
    [~, f1, local_P1, local_Q1, local_U1,Ess1] = subnormal1( ...
        rho_s,rho_t, zP, zQ, zU, uP, uQ, uU, ...    % ADMM参数
        lam1 ...                        % 边界初值
    );
% fprintf('  subnormal1 返回：f1=%.4f, local_P1(1) = %.4f\n', f1, local_P1(1));
    %% 子网二：故障网（21节点入口边界），带惩罚项与乘子变量
    [~, f2, local_P2, local_Q2, local_U2,beta2,alphaIEEE_12,Ess2] = subguzhang( ...
        rho_s, rho_t, zP, zQ, zU, lamdaP, lamdaQ, lamdaU, ...    % ADMM参数
         lam2  ...                        % 边界初值
    );
% fprintf('  subguzhang 返回：f2=%.4f, local_P2(1) = %.4f\n', f2, local_P2(1));

% -------- Step 2 SOC一致性主变量更新 --------
    % z_new = zeros(size(z_soc));
    % count = zeros(size(z_soc));
    % z_new(1:Ess_Num1,:) = Ess1 + lam1/rho_t; count(1:Ess_Num1,:) = count(1:Ess_Num1,:) + 1;
    % z_new(1:Ess_Num,:) = z_new(1:Ess_Num,:) + Ess2 + lam2/rho_t; count(1:Ess_Num,:) = count(1:Ess_Num,:) + 1;
    % z_soc = z_new ./ max(count,1);
    % % 边界条件
    % z_soc(:,1) = 0.6;
    % z_soc(:,end) = min(max(z_soc(:,end),0.4),0.6);
    % z_soc = min(max(z_soc,0.2),0.9);
%%SOC对偶变量更新 --------
    % lam1 = lam1 + rho_t*(Ess1 - z_soc(1:Ess_Num1,:));
    % lam2 = lam2 + rho_t*(Ess2 - z_soc(1:Ess_Num,:));
   %% 对故障网的松弛二进制变量做指数移动平均平滑
    beta_smoothed = eta_beta * beta2 + (1-eta_beta) * beta_prev;
    alpha_smoothed = eta_beta * alphaIEEE_12 + (1-eta_beta) * alphaIEEE_prev;
    %% 局部变量（边界值）汇总
    p1 = local_P1(:).'; q1 = local_Q1(:).'; u1 = local_U1(:).';
    p2 = local_P2(:).'; q2 = local_Q2(:).'; u2 = local_U2(:).';

    %% 共识变量更新（简单加权平均/可自定义）
    zP_prev = zP; zQ_prev = zQ; zU_prev = zU;
    zP = 0.5*(p1 + p2);
    zQ = 0.5*(q1 + q2);
    zU = 0.5*(u1 + u2);
%% 故障网松弛变量保存
    beta2_hist{k} = beta2;
    alphaIEEE_12_hist{k} = alphaIEEE_12;
    %% 原始残差r_k: 两个子网边界变量差异
    r_k = [p1 - p2; q1 - q2; u1 - u2];
    r_norm = norm(r_k);
    r_hist = [r_hist, r_norm];

    %% 对偶残差s_k: 共识变量变化
    s_k = [zP - zP_prev; zQ - zQ_prev; zU - zU_prev];
    s_norm = rho_s * norm(s_k);
    s_hist = [s_hist, s_norm];

    obj_hist = [obj_hist, f1 + f2];

    %% 乘子项更新（ADMM标准式）
    uP = uP + p1 - zP;
    uQ = uQ + q1 - zQ;
    uU = uU + u1 - zU;

    lamdaP = lamdaP + p2 - zP;
    lamdaQ = lamdaQ + q2 - zQ;
    lamdaU = lamdaU + u2 - zU;
    %更新数值
    beta_prev = beta_smoothed;
    alphaIEEE_prev = alpha_smoothed;

    %% 输出迭代信息
    % fprintf('  外层iter%d：原始残差 ||r|| = %.4e, 对偶残差 ||s|| = %.4e, 总目标 = %.4e\n', r_norm, s_norm, f1+f2);
 %%SOC一致性收敛监控 ---------
    % prim_soc = norm(Ess1 - z_soc(1:Ess_Num1,:), 'fro') + norm(Ess2 - z_soc(1:Ess_Num,:), 'fro');
    % dual_soc = norm(z_soc - z_new./max(count,1), 'fro');
    % prim_soc_hist = [prim_soc_hist, prim_soc];
    % dual_soc_hist = [dual_soc_hist, dual_soc];


    %% 收敛判据
    if r_norm < tol && s_norm < tol
        fprintf('外层ADMM收敛于第 %d 轮\n', k);
        break;
    end
end
%% 收敛后：对松弛二进制变量取阈值（0.5）并可选地做可行化重求解
beta_final = double(beta_smoothed > 0.5);
alpha_final = double(alpha_smoothed > 0.5);
obj=value(obj_hist);
%% 结果展示
fprintf('\n---- ADMM结束 ----\n');
fprintf('最终残差: 原始 %.4e, 对偶 %.4e\n', r_norm, s_norm);
fprintf('边界共识（功率/电压）：P=%.4f Q=%.4f U=%.4f\n', zP(1), zQ(1), zU(1));
fprintf('正常网边界（P Q U）：%.4f %.4f %.4f\n', p1(1), q1(1), u1(1));
fprintf('故障网边界（P Q U）：%.4f %.4f %.4f\n', p2(1), q2(1), u2(1));
figure;
subplot(3,1,1); plot(r_hist,'-o'); ylabel('原始残差||r||');
subplot(3,1,2); plot(s_hist,'-o'); ylabel('对偶残差||s||');
subplot(3,1,3); plot(obj,'-o'); ylabel('总目标');
% subplot(4,1,4); plot(prim_soc_hist,'-o'); ylabel('SOC一致性残差');
xlabel('ADMM次数');
