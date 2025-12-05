clear; clc;

%% 参数与初始化
max_iter = 60;
tol = 1e-2;          % 收敛阈值
nbeta = 36;            % 故障网中松弛开关数量示例
T = 2;               % 时间步长
rho_s = 0.1;
rho_t = 0.1;
eta_beta = 0.1;        % 残差/变量平滑因子（指数移动平均），越小平滑越强
lambda_beta = 300;     % 使 beta 向 0/1 集中之凸二次

% 自适应规则参数
mu_adapt = 10;        % 比例阈值
tau_inc = 2.5;          % 增大因子
tau_dec = 2.5;          % 减小因子
%% 初始化 松弛二进制变量及其平滑版本
% beta_prev = 0.5*ones(nbeta,1);    % 上一轮 beta（平滑用）
beta0=[0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,0,0,0,0];
beta_init = repmat(beta0(:), 1, T);
beta_prev = beta_init; 
alphaIEEE_prev = 0.5*ones(1,T);               % 上一轮 alpha（示例单个）
beta_smoothed = beta_prev;
alpha_smoothed = alphaIEEE_prev;

% 全局-边界共识变量（24节点出口/21节点入口）
zP = zeros(1, T); zQ = zeros(1, T); zU = zeros(1, T);
% 局部子网边界变量
local_P1 = zeros(1, T); local_Q1 = zeros(1, T); local_U1 = ones(1, T);  % 正常网
local_P2 = zeros(1, T); local_Q2 = zeros(1, T); local_U2 = ones(1, T);  % 故障网

% ADMM乘子初始化(边界)
uP = zeros(1, T); uQ = zeros(1, T); uU = zeros(1, T);
lamdaP = zeros(1, T); lamdaQ = zeros(1, T); lamdaU = zeros(1, T);
%内层时间耦合
Ess_Num1=3; Ess_Num=3;
z_soc1=0.6*ones(Ess_Num1,T+1);lam1=zeros(Ess_Num1,T+1);
z_soc =0.6*ones(Ess_Num,T+1); lam2=zeros(Ess_Num,T+1);
% 记录残差与目标
r_hist = []; s_hist = []; obj_hist = [];
beta2_hist = [];f1_hist = [];f2_hist = [];
 Unode1_hist =[];Unode2_hist =[];
 % prim_res_hist=[];dual_res_hist=[];
fprintf('ADMM主循环空间解耦（边界联络线共识协调）\n');

for k = 1:max_iter
    fprintf('---------- ADMM 第 %d 轮 ----------\n', k);

    %% 子网一：正常网（24节点出口边界），带惩罚项与乘子变量
    [Unode1, f1_val, local_P1, local_Q1, local_U1,Pdis1,z_soc1_hist] = subnormal1125(rho_s,rho_t,zP, zQ, zU, uP, uQ, uU,z_soc1, lam1);
                                     
   

    %% 子网二：故障网（21节点入口边界），带惩罚项与乘子变量
    [Unode2, f2_val, local_P2, local_Q2, local_U2,beta2,alphaIEEE_12,z2,P_disf,z_soc_hist] = subguzhang1125( rho_s,rho_t,zP, zQ, zU, lamdaP, lamdaQ, lamdaU,z_soc, lam2,lambda_beta,beta_init);  


   %% 对故障网的松弛二进制变量做指数移动平均平滑
    beta_smoothed = eta_beta * beta2 + (1-eta_beta) * beta_prev;
    alpha_smoothed = eta_beta * alphaIEEE_12 + (1-eta_beta) * alphaIEEE_prev;
    %% 局部变量（边界值）汇总
    p1 = local_P1(:).'; q1 = local_Q1(:).'; u1 = local_U1(:).';
    p2 = local_P2(:).'; q2 = local_Q2(:).'; u2 = local_U2(:).';

    %% 共识变量更新（简单加权平均/可自定义）
    zP_prev = zP; zQ_prev = zQ; zU_prev = zU;
    zP = 0.5*(p1 + p2 +((uP + lamdaP)/rho_s));
    zQ = 0.5*(q1 + q2 +((uQ + lamdaQ)/rho_s));
    zU = 0.5*(u1 + u2 +((uU + lamdaU)/rho_s));

    % zP = 0.5*(p1 + p2 + uP + lamdaP);
    % zQ = 0.5*(q1 + q2 + uQ + lamdaQ);
    % zU = 0.5*(u1 + u2 + uU + lamdaU);
%% 故障网松弛变量保存
    beta2_hist{k} = beta2;
    alphaIEEE_12_hist{k} = alphaIEEE_12;
    %% 原始残差r_k: 两个子网边界变量差异
    r_k = [p1 - zP; p2 - zP;q1 - zQ; q2 - zQ;u1 - zU;u2 - zU];
    r_norm = norm(r_k);
    r_hist = [r_hist, r_norm];

    %% 对偶残差s_k: 共识变量变化
    s_k = [zP - zP_prev; zQ - zQ_prev; zU - zU_prev];
    s_norm = rho_s * norm(s_k);
    s_hist = [s_hist, s_norm];
    obj_hist = [obj_hist, f1_val + f2_val]; 
    f1_hist = [f1_hist,f1_val];f2_hist = [f2_hist,f2_val];
     Unode1_hist =[Unode1_hist,Unode1];Unode2_hist =[Unode2_hist,Unode2];

    %% 乘子项更新（ADMM标准式）
    uP = uP + (p1 - zP);
    uQ = uQ + (q1 - zQ);
    uU = uU + (u1 - zU);

    lamdaP = lamdaP + (p2 - zP);
    lamdaQ = lamdaQ + (q2 - zQ);
    lamdaU = lamdaU + (u2 - zU);
    %更新数值
    beta_prev = beta_smoothed;
    alphaIEEE_prev = alpha_smoothed;
 %% 自适应惩罚参数更新
    if r_norm > mu_adapt * s_norm
        rho_old_s = rho_s; rho_old_t = rho_t;
        rho_s = rho_s * tau_inc;
        rho_t = rho_t * tau_inc;
        % 缩放乘子以匹配 rho 变化
        uP = uP / tau_inc; uQ = uQ / tau_inc; uU = uU / tau_inc;
        lamdaP = lamdaP / tau_inc; lamdaQ = lamdaQ / tau_inc; lamdaU = lamdaU / tau_inc;
        fprintf('  Adaptive rho: 增大 rho_s %g -> %g, rho_t %g -> %g\n', rho_old_s, rho_s, rho_old_t, rho_t);
    elseif s_norm > mu_adapt * r_norm
        rho_old_s = rho_s; rho_old_t = rho_t;
        rho_s = max(1, rho_s / tau_dec);
        rho_t = max(1, rho_t / tau_dec);
        uP = uP * tau_dec; uQ = uQ * tau_dec; uU = uU * tau_dec;
        lamdaP = lamdaP * tau_dec; lamdaQ = lamdaQ * tau_dec; lamdaU = lamdaU * tau_dec;
        fprintf('  Adaptive rho: 减小 rho_s %g -> %g, rho_t %g -> %g\n', rho_old_s, rho_s, rho_old_t, rho_t);
    end
    %% 输出迭代信息
     % fprintf('  原始残差 ||r|| = %.4e, 对偶残差 ||s|| = %.4e, 目标 = %.4e\n', r_norm, s_norm, obj_hist(end));

    %% 收敛判据
    if r_norm < tol && s_norm < tol
        fprintf('ADMM收敛于第 %d 轮\n', k);
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
% subplot(5,1,4); plot(prim_res1_hist,'-o'); ylabel('内层原始残差||r||');
% subplot(5,1,5); plot(prim_res1_hist,'-o'); ylabel('内层对偶残差||s||');
xlabel('ADMM次数');
