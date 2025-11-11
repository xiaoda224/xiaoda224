function [Unode2, f2,local_P2,local_Q2,local_U2,P_chf,P_disf]=subguzhang(rho_s,rho_t,zP,zQ,zU,uP,uQ,uU,z_soc, lam2,PIEEE_1f,QIEEE_1f,UjIEEE_1f)
max_iter = 50;
max_iiter = 50;
T = 2;               % 时间步长
rho_s = 5000;
rho_t = 3000;
tol=1e-4;

bus_Num=32;
Ess_Num=3;PV_Num=5;
Sijmax=10;
Lijmax=10;
Umax=1.05^2;
Umin=0.85^2;
Ubase=1.03^2;
%% 故障网光伏参数
Spv=[0.02;0.03;0.05;0.02;0.03];    %光伏变化T=2
f_pv_bian=[1,1,1,1,1,1];  
PVf=Spv*f_pv_bian;  %光伏出力
%% 节点负荷实时数据
%配电网
Pnet(1:bus_Num,T)=[0.1, 0.09, 0.12, 0.06, 0.06, 0.2, 0.2, 0.06, 0.06, 0.045, 0.06, 0.06, 0.12, 0.06, 0.06, 0.06, 0.09, 0.09, 0.09, 0.09, 0.09, 0.09, 0.42, 0.42, 0.06, 0.06, 0.06, 0.12, 0.2, 0.15, 0.21, 0.06]/10;
Qnet(1:bus_Num,T)=[0.06, 0.04, 0.08, 0.03, 0.02, 0.1, 0.1, 0.02, 0.02, 0.03, 0.035, 0.035, 0.08, 0.01, 0.01, 0.02, 0.04, 0.04, 0.04, 0.04, 0.04, 0.05, 0.2, 0.2, 0.025, 0.025, 0.02, 0.07, 0.6, 0.07, 0.1, 0.04]/10;
Pj(1:bus_Num,T)=Pnet(1:bus_Num,T);
Qj(1:bus_Num,T)=Qnet(1:bus_Num,T);
%% 线路阻抗参数
%配电网
L=32;L2=36;
Ldata=[0.0922,0.047;  0.493,0.2511;  0.366,0.1864;  0.3811,0.1941;  0.0819,0.707;  0.1872,0.6188;  0.7114,0.2351;  1.03,0.74;       1.044,0.74;     0.1966,0.065;  0.3744,0.1238;
                 1.468,1.155;  0.5416,0.7129; 0.591,0.526;   0.7463,0.545;   1.289,1.721;   0.372,0.574;    0.164,0.1565;   1.5042,1.3554;   0.4095,0.4784;  0.7089,0.9373; 0.4512, 0.3083;
                 0.898,0.7091; 0.896,0.7011;  0.203,0.1034;  0.2842,0.1447;  1.059,0.9337;  0.8042,0.7006;  0.5075,0.2585;  0.9744,0.963;    0.3105, 0.3619; 0.341, 0.5362;    0.5, 0.5;     0.5,0.5;
                 0.5,0.5;     0.5,0.5];

Rij(1:L,1)=Ldata(1:L,1)/16.02756;
Xij(1:L,2)=Ldata(1:L,2)/16.02756;
Rijf(1:L2,1)=Ldata(1:L2,1)/16.02756;%故障网内部联络线开关4个
Xijf(1:L2,2)=Ldata(1:L2,2)/16.02756;
RIEEE_1=0.5/16.02756;XIEEE_1=0.5/16.02756;
%故障配电网
PIEEE_1=sdpvar(1,T);QIEEE_1=sdpvar(1,T);LIEEE_1=sdpvar(1,T);%正常网24节点向外发出的给故障网的
PIEEE_1f=sdpvar(1,T);QIEEE_1f=sdpvar(1,T);%故障网接受的参数
alphaIEEE_1=sdpvar(1,T);%故障网与正常网连的边界值，故障网节点21
UjIEEE_1f=sdpvar(1,T);%故障网与正常网相连的边界电压
M=1000;
Pijf=sdpvar(L2,T);%故障配电网内部线路
Qijf=sdpvar(L2,T);
Ujf=sdpvar(bus_Num,T);UTf=sdpvar(1,T);
Lijf=sdpvar(L2,T);
Pdecf=sdpvar(5,T);%光伏削减量
Qdecf=sdpvar(5,T);
Pchf=sdpvar(Ess_Num,T);uchf=sdpvar(Ess_Num,T);
Pdisf=sdpvar(Ess_Num,T);udisf=sdpvar(Ess_Num,T);
Ess=sdpvar(Ess_Num,T+1);
Pre=sdpvar(32,T);  %表示负荷是否拾取
Qre=sdpvar(32,T);  %表示负荷是否拾取
beta=sdpvar(36,T);%网架结构，开关状态，0断 1接
z=sdpvar(32,T);
%储能充放电参数
P_essmax=2*1.5/10; % 充放最大功率  基值是10MW
%% 开关初始状态
beta0=zeros(L2,1);
beta0(1:L2)=[0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,0,0,0,0];

Ess_var = 0.6*ones(Ess_Num, T+1);
lam2 = zeros(Ess_Num, T+1);
Ess_hist = zeros(Ess_Num, T+1, max_iiter);
for inner_iter = 1:max_iiter
    [Unode2, f2, local_P2, local_Q2, local_U2, P_chf, P_disf] = ...
        subguzhang(Ess_var,zP,zQ,zU,uP,uQ,uU, lam2,PIEEE_1f,QIEEE_1f,UjIEEE_1f);
     r = zeros(Ess_Num, T);
    for t = 1:T
        r(:,t) = Ess_var(:,t+1) - (Ess_var(:,t) + 0.9*P_chf(:,t) - 1.11*P_disf(:,t));
    end
    Ess_var_old = Ess_var;
    for t = 1:T
        Ess_var(:, t+1) = Ess_var(:,t) + 0.9*P_chf(:,t) - 1.11*P_disf(:,t) - lam2(:,t)/rho_t;
    end
    Ess_var = min(max(Ess_var, 0.2), 0.9);
    Ess_var(:,1) = 0.6;
    Ess_var(:,end) = min(max(Ess_var(:,end), 0.4), 0.6);

    for t = 1:T
        lam2(:,t) = lam2(:,t) + rho_t*r(:,t);
    end

    Ess_hist(:,:,inner_iter) = Ess_var;
%% 优化目标
for t=1:T
  Ess(:,t)=z_soc(:,t); Ess(:,t+1)=z_soc(:,t+1);
 Cl=5000;Cp=4000;Cpv=7000;Cb=0.5;%有功缩减成本系数
 Ce=4000;%储能系统放电成本系数

  Rf=repmat(Rijf,1,T);
  lostf=Lijf.*Rf;
  loadr=z.*Pj;
  f_loadr=Cl*sum(sum(loadr(:,1:T)));
  f_l=Cp*sum(sum(lostf(:,1:T)));
  f_pv=Cpv*(sum(Pdecf(:,1:T)));
  f_ess=Ce*sum(Pchf(:)) + Ce*sum(Pdisf(:));
  f_b=0.5*rho_s*(PIEEE_1f(1,t)+zP(t)+uP(t))^2+0.5*rho_s*(QIEEE_1f(1,t)+zQ(t)+uQ(t))^2+0.5*rho_s*(UjIEEE_1f(1,t)+zU(t)+uU(t))^2;
  % f_d=0.5*rho_t* (Ess - z_soc +lam2).^2 + 0.5*rho_t *(Ess - z_soc + lam2).^2;
  % f_d=0.5*rho_t* (r(:,t+1) - r(:,t) +lam2).^2 + 0.5*rho_t *(r(:,t+1) - r(:,t) + lam2).^2;
  f2 = f_loadr+ f_l+ f_pv + f_ess + f_b  ;%+ f_d;

end
 
  for iteri = 1:max_iter
%%故障网约束
 F=[];
 for t=1:T
     %储能约束
      F = [F;
          udisf(:,t)+uchf(:,t)<=1;       udisf(:,t)+uchf(:,t)>=0;%表示充电，放电，不充不放三种状态
          Pdisf(:,t)>=0;Pdisf(:,t)<=udisf(:,t)*P_essmax;%储能放电功率约束
          Pchf(:,t)>=0;Pchf(:,t)<=uchf(:,t)*P_essmax;%储能放电功率约束
          % Ess1(:,1)==0.6;     
          % Ess1(:,t)>=0.2; Ess1(:,t)<=0.9;
          % Ess1(:,T+1)>=0.4; Ess1(:,T+1)<=0.6; 
          
          Ess(:,t+1)==Ess(:,t)+0.9*Pchf(:,t)-1.11*Pdisf(:,t); 
          
          ];
F=[F;
        0<=beta(:);beta(:)<=1;0<=alphaIEEE_1(:);alphaIEEE_1(:)<=1;0<=z(:);z(:)<=1;
        % UjIEEE_1f(1,t)==Uj1(24,t);
        0<=Lijf(:,t);     Lijf(:,t)<=Lijmax;  %线路电流约束
        Umin<=Ujf(:,t);    Ujf(:,t)<=Umax;   %节点电压约束
        0<=Pdecf(:,t);     Pdecf(:,t)<=PVf(:,t);
        Pre(:,t)==z(:,t).*Pj(:,t);Qre(:,t)==z(:,t).*Qj(:,t);
        Pijf(:,t)>=-beta(:,t)*Sijmax;Pijf(:,t)<=beta(:,t)*Sijmax;
        Qijf(:,t)>=-beta(:,t)*Sijmax;Qijf(:,t)<=beta(:,t)*Sijmax;
        Lijf(:,t)>=-beta(:,t)*Lijmax;Lijf(:,t)<=beta(:,t)*Lijmax;    

    -Sijmax*(alphaIEEE_1(1,t))<=PIEEE_1(1,t);  PIEEE_1(1,t)<=Sijmax*(alphaIEEE_1(1,t));
    -Sijmax*(alphaIEEE_1(1,t))<=QIEEE_1(1,t);  QIEEE_1(1,t)<=Sijmax*(alphaIEEE_1(1,t));
    -Lijmax*(alphaIEEE_1(1,t))<=LIEEE_1(1,t);  LIEEE_1(1,t)<=Lijmax*(alphaIEEE_1(1,t));
    
        Pijf(1,t)==0; Qijf(1,t)==0;Lijf(1,t)==0;%因为故障了
        
        Pijf(1,t)-Lijf(1,t)*Rijf(1)-Pj(1,t)+Pre(1,t)==Pijf(2,t)+Pijf(18,t);
        Pijf(2,t)-Lijf(2,t)*Rijf(2)-Pj(2,t)+Pre(2,t)==Pijf(3,t)+Pijf(22,t);
        Pijf(3,t)-Lijf(3,t)*Rijf(3)-Pj(3,t)+Pre(3,t)==Pijf(4,t);
        Pijf(4,t)-Lijf(4,t)*Rijf(4)-Pj(4,t)+Pre(4,t)==Pijf(5,t);
        Pijf(5,t)-Lijf(5,t)*Rijf(5)-Pj(5,t)+Pre(5,t)==Pijf(6,t)+Pijf(25,t);
        Pijf(6,t)-Lijf(6,t)*Rijf(6)-Pj(6,t)+Pre(6,t)==Pijf(7,t);
        % Pij(7,t)-Lij(7,t)*Rijf(7)+Pij(34,t)-Lij(34,t)*Rijf(34)-Pj(7,t)+Pre(7,t)==Pij(8,t);
        Pijf(7,t)-Lijf(7,t)*Rijf(7)+Pijf(34,t)-Lijf(34,t)*Rijf(34)-Pj(7,t)+Pdisf(1,t)-Pchf(1,t)+Pre(7,t)==Pijf(8,t);%7节点处有储能
%       Ess(1,t+1)==Ess(1,t)+0.9*P_ch(1,t)-1.11*P_dis(1,t);
        Pijf(8,t)-Lijf(8,t)*Rijf(8)-Pj(8,t)+PVf(1,t)-Pdecf(1,t)+Pre(8,t)==Pijf(9,t);  
        Pijf(9,t)-Lijf(9,t)*Rijf(9)-Pj(9,t)+Pre(9,t)==Pijf(10,t);
        Pijf(10,t)-Lijf(10,t)*Rijf(10)-Pj(10,t)+Pre(10,t)==Pijf(11,t);
        Pijf(11,t)-Lijf(11,t)*Rijf(11)+Pijf(33,t)-Lijf(33,t)*Rijf(33)-Pj(11,t)+Pre(11,t)==Pijf(12,t);
        Pijf(12,t)-Lijf(12,t)*Rijf(12)-Pj(12,t)+PVf(2,t)-Pdecf(2,t)+Pre(12,t)==Pijf(13,t);%12节点含有光伏
        Pijf(13,t)-Lijf(13,t)*Rijf(13)-Pj(13,t)+Pre(13,t)==Pijf(14,t);
        Pijf(14,t)-Lijf(14,t)*Rijf(14)-Pj(14,t)+Pre(14,t)==Pijf(15,t);
        Pijf(15,t)-Lijf(15,t)*Rijf(15)-Pj(15,t)+Pdisf(2,t)-Pchf(2,t)+Pre(15,t)==Pijf(16,t);%ess
        Pijf(16,t)-Lijf(16,t)*Rijf(16)-Pj(16,t)+Pre(16,t)==Pijf(17,t);
        Pijf(17,t)-Lijf(17,t)*Rijf(17)+Pijf(36,t)-Lijf(36,t)*Rijf(36)-Pj(17,t)+Pre(17,t)==0;
        Pijf(18,t)-Lijf(18,t)*Rijf(18)-Pj(18,t)+Pre(18,t)==Pijf(19,t);
        Pijf(19,t)-Lijf(19,t)*Rijf(19)-Pj(19,t)+PVf(3,t)-Pdecf(3,t)+Pre(19,t)==Pijf(20,t);%19节点含有光伏
        Pijf(20,t)-Lijf(20,t)*Rijf(20)-Pj(20,t)+Pre(20,t)==Pijf(21,t)+Pijf(34,t);
        PIEEE_1f(1,t)+Pijf(21,t)-Lijf(21,t)*Rijf(21)-Pj(21,t)+Pre(21,t)==Pijf(33,t);%21节点连接正常配网中的24节点
        Pijf(22,t)-Lijf(22,t)*Rijf(22)-Pj(22,t)+Pre(22,t)==Pijf(23,t);
        Pijf(23,t)-Lijf(23,t)*Rijf(23)-Pj(23,t)+PVf(4,t)-Pdecf(4,t)+Pre(23,t)==Pijf(24,t);
        Pijf(24,t)-Lijf(24,t)*Rijf(24)-Pj(24,t)+Pre(24,t)==Pijf(35,t);
        Pijf(25,t)-Lijf(25,t)*Rijf(25)-Pj(25,t)+Pre(25,t)==Pijf(26,t);
        Pijf(26,t)-Lijf(26,t)*Rijf(26)-Pj(26,t)+Pre(26,t)==Pijf(27,t);
        Pijf(27,t)-Lijf(27,t)*Rijf(27)-Pj(27,t)+Pre(27,t)==Pijf(28,t);
        Pijf(28,t)-Lijf(28,t)*Rijf(28)+Pijf(35,t)-Lijf(35,t)*Rijf(35)-Pj(28,t)+Pdisf(3,t)-Pchf(3,t)+Pre(28,t)==Pijf(29,t);%ess
        Pijf(29,t)-Lijf(29,t)*Rijf(29)-Pj(29,t)+Pre(29,t)==Pijf(30,t);
        Pijf(30,t)-Lijf(30,t)*Rijf(30)-Pj(30,t)+Pre(30,t)==Pijf(31,t);
        Pijf(31,t)-Lijf(31,t)*Rijf(31)-Pj(31,t)+PVf(5,t)-Pdecf(5,t)+Pre(31,t)==Pijf(32,t);
        Pijf(32,t)-Lijf(32,t)*Rijf(32)-Pj(32,t)+Pre(32,t)==Pijf(36,t);
        
        Qijf(1,t)-Lijf(1,t)*Xijf(1)-Qj(1,t)+Qre(1,t)==Qijf(2,t)+Qijf(18,t);
        Qijf(2,t)-Lijf(2,t)*Xijf(2)-Qj(2,t)+Qre(2,t)==Qijf(3,t)+Qijf(22,t);
        Qijf(3,t)-Lijf(3,t)*Xijf(3)-Qj(3,t)+Qre(3,t)==Qijf(4,t);
        Qijf(4,t)-Lijf(4,t)*Xijf(4)-Qj(4,t)+Qre(4,t)==Qijf(5,t);
        Qijf(5,t)-Lijf(5,t)*Xijf(5)-Qj(5,t)+Qre(5,t)==Qijf(6,t)+Qijf(25,t);
        Qijf(6,t)-Lijf(6,t)*Xijf(6)-Qj(6,t)+Qre(6,t)==Qijf(7,t);
        Qijf(7,t)-Lijf(7,t)*Xijf(7)+Qijf(34,t)-Lijf(34,t)*Xijf(34)-Qj(7,t)+Qre(7,t)==Qijf(8,t);
        Qijf(8,t)-Lijf(8,t)*Xijf(8)-Qj(8,t)+Qdecf(1,t)+Qre(8,t)==Qijf(9,t);  
        Qdecf(1,t)<=(PVf(1,t)-Pdecf(1,t))*0.32868421;  %PV功率因数约束，后面要考虑容量约束
        Qdecf(1,t)>=(Pdecf(1,t)-PVf(1,t))*0.32868421;
        Qijf(9,t)-Lijf(9,t)*Xijf(9)-Qj(9,t)+Qre(9,t)==Qijf(10,t);
        Qijf(10,t)-Lijf(10,t)*Xijf(10)-Qj(10,t)+Qre(10,t)==Qijf(11,t);
        Qijf(11,t)-Lijf(11,t)*Xijf(11)+Qijf(33,t)-Lijf(33,t)*Xijf(33)-Qj(11,t)+Qre(11,t)==Qijf(12,t);
        Qijf(12,t)-Lijf(12,t)*Xijf(12)-Qj(12,t)+Qdecf(2,t)+Qre(12,t)==Qijf(13,t);
        Qdecf(2,t)<=(PVf(2,t)-Pdecf(2,t))*0.32868421;  %PV功率因数约束，后面要考虑容量约束
        Qdecf(2,t)>=(Pdecf(2,t)-PVf(2,t))*0.32868421;

        Qijf(13,t)-Lijf(13,t)*Xijf(13)-Qj(13,t)+Qre(13,t)==Qijf(14,t);
        Qijf(14,t)-Lijf(14,t)*Xijf(14)-Qj(14,t)+Qre(14,t)==Qijf(15,t);
        Qijf(15,t)-Lijf(15,t)*Xijf(15)-Qj(15,t)+Qre(15,t)==Qijf(16,t);
        Qijf(16,t)-Lijf(16,t)*Xijf(16)-Qj(16,t)+Qre(16,t)==Qijf(17,t);
        Qijf(17,t)-Lijf(17,t)*Xijf(17)+Qijf(36,t)-Lijf(36,t)*Xijf(36)-Qj(17,t)+Qre(17,t)==0;
        Qijf(18,t)-Lijf(18,t)*Xijf(18)-Qj(18,t)+Qre(18,t)==Qijf(19,t);
        Qijf(19,t)-Lijf(19,t)*Xijf(19)-Qj(19,t)+Qdecf(3,t)+Qre(19,t)==Qijf(20,t);
        Qdecf(3,t)<=(PVf(3,t)-Pdecf(3,t))*0.32868421;  %PV功率因数约束，后面要考虑容量约束
        Qdecf(3,t)>=(Pdecf(3,t)-PVf(3,t))*0.32868421;

        Qijf(20,t)-Lijf(20,t)*Xijf(20)-Qj(20,t)+Qre(20,t)==Qijf(21,t)+Qijf(34,t);
        QIEEE_1f(1,t)+Qijf(21,t)-Lijf(21,t)*Xijf(21)-Qj(21,t)+Qre(21,t)==Qijf(33,t);
        Qijf(22,t)-Lijf(22,t)*Xijf(22)-Qj(22,t)+Qre(22,t)==Qijf(23,t);
        Qijf(23,t)-Lijf(23,t)*Xijf(23)-Qj(23,t)+Qdecf(4,t)+Qre(23,t)==Qijf(24,t);
        Qdecf(4,t)<=(PVf(4,t)-Pdecf(4,t))*0.32868421;  %PV功率因数约束，后面要考虑容量约束
        Qdecf(4,t)>=(Pdecf(4,t)-PVf(4,t))*0.32868421;

        Qijf(24,t)-Lijf(24,t)*Xijf(24)-Qj(24,t)+Qre(24,t)==Qijf(35,t);
        Qijf(25,t)-Lijf(25,t)*Xijf(25)-Qj(25,t)+Qre(25,t)==Qijf(26,t);
        Qijf(26,t)-Lijf(26,t)*Xijf(26)-Qj(26,t)+Qre(26,t)==Qijf(27,t);
        Qijf(27,t)-Lijf(27,t)*Xijf(27)-Qj(27,t)+Qre(27,t)==Qijf(28,t);
        Qijf(28,t)-Lijf(28,t)*Xijf(28)+Qijf(35,t)-Lijf(35,t)*Xijf(35)-Qj(28,t)+Qre(28,t)==Qijf(29,t);
        Qijf(29,t)-Lijf(29,t)*Xijf(29)-Qj(29,t)+Qre(29,t)==Qijf(30,t);
        Qijf(30,t)-Lijf(30,t)*Xijf(30)-Qj(30,t)+Qre(30,t)==Qijf(31,t);
        Qijf(31,t)-Lijf(31,t)*Xijf(31)-Qj(31,t)+Qdecf(5,t)+Qre(31,t)==Qijf(32,t);
        Qdecf(5,t)<=(PVf(5,t)-Pdecf(5,t))*0.32868421;  %PV功率因数约束，后面要考虑容量约束
        Qdecf(5,t)>=(Pdecf(5,t)-PVf(5,t))*0.32868421;
        Qijf(32,t)-Lijf(32,t)*Xijf(32)-Qj(32,t)+Qre(32,t)==Qijf(36,t);
        
        UTf(1,t)==0;%主变故障变为0
        Ujf(1,t)>=-M*(1-beta(1,t))+UTf(1,t)-2*(Pijf(1,t)*Rijf(1)+Qijf(1,t)*Xijf(1))+Lijf(1,t)*(Rijf(1)^2+Xijf(1)^2);%%%UTf变0
        Ujf(1,t)<=M*(1-beta(1,t)) +UTf(1,t)-2*(Pijf(1,t)*Rijf(1)+Qijf(1,t)*Xijf(1))+Lijf(1,t)*(Rijf(1)^2+Xijf(1)^2);
        Ujf(2,t)>=-M*(1-beta(2,t))+Ujf(1,t)-2*(Pijf(2,t)*Rijf(2)+Qijf(2,t)*Xijf(2))+Lijf(2,t)*(Rijf(2)^2+Xijf(2)^2);
        Ujf(2,t)<=M*(1-beta(2,t))+Ujf(1,t)-2*(Pijf(2,t)*Rijf(2)+Qijf(2,t)*Xijf(2))+Lijf(2,t)*(Rijf(2)^2+Xijf(2)^2);
        Ujf(3,t)>=-M*(1-beta(3,t))+Ujf(2,t)-2*(Pijf(3,t)*Rijf(3)+Qijf(3,t)*Xijf(3))+Lijf(3,t)*(Rijf(3)^2+Xijf(3)^2);
        Ujf(3,t)<=M*(1-beta(3,t))+Ujf(2,t)-2*(Pijf(3,t)*Rijf(3)+Qijf(3,t)*Xijf(3))+Lijf(3,t)*(Rijf(3)^2+Xijf(3)^2);
        Ujf(4,t)>=-M*(1-beta(4,t))+Ujf(3,t)-2*(Pijf(4,t)*Rijf(4)+Qijf(4,t)*Xijf(4))+Lijf(4,t)*(Rijf(4)^2+Xijf(4)^2);
        Ujf(4,t)<=M*(1-beta(4,t))+Ujf(3,t)-2*(Pijf(4,t)*Rijf(4)+Qijf(4,t)*Xijf(4))+Lijf(4,t)*(Rijf(4)^2+Xijf(4)^2);
        Ujf(5,t)>=-M*(1-beta(5,t))+Ujf(4,t)-2*(Pijf(5,t)*Rijf(5)+Qijf(5,t)*Xijf(5))+Lijf(5,t)*(Rijf(5)^2+Xijf(5)^2);
        Ujf(5,t)<=M*(1-beta(5,t))+Ujf(4,t)-2*(Pijf(5,t)*Rijf(5)+Qijf(5,t)*Xijf(5))+Lijf(5,t)*(Rijf(5)^2+Xijf(5)^2);
        Ujf(6,t)>=-M*(1-beta(6,t))+Ujf(5,t)-2*(Pijf(6,t)*Rijf(6)+Qijf(6,t)*Xijf(6))+Lijf(6,t)*(Rijf(6)^2+Xijf(6)^2);
        Ujf(6,t)<=M*(1-beta(6,t))+Ujf(5,t)-2*(Pijf(6,t)*Rijf(6)+Qijf(6,t)*Xijf(6))+Lijf(6,t)*(Rijf(6)^2+Xijf(6)^2);
        Ujf(7,t)>=-M*(1-beta(7,t))+Ujf(6,t)-2*(Pijf(7,t)*Rijf(7)+Qijf(7,t)*Xijf(7))+Lijf(7,t)*(Rijf(7)^2+Xijf(7)^2);
        Ujf(7,t)<=M*(1-beta(7,t))+Ujf(6,t)-2*(Pijf(7,t)*Rijf(7)+Qijf(7,t)*Xijf(7))+Lijf(7,t)*(Rijf(7)^2+Xijf(7)^2);
        Ujf(8,t)>=-M*(1-beta(8,t))+Ujf(7,t)-2*(Pijf(8,t)*Rijf(8)+Qijf(8,t)*Xijf(8))+Lijf(8,t)*(Rijf(8)^2+Xijf(8)^2);
        Ujf(8,t)<=M*(1-beta(8,t))+Ujf(7,t)-2*(Pijf(8,t)*Rijf(8)+Qijf(8,t)*Xijf(8))+Lijf(8,t)*(Rijf(8)^2+Xijf(8)^2);
        Ujf(9,t)>=-M*(1-beta(9,t))+Ujf(8,t)-2*(Pijf(9,t)*Rijf(9)+Qijf(9,t)*Xijf(9))+Lijf(9,t)*(Rijf(9)^2+Xijf(9)^2);
        Ujf(9,t)<=M*(1-beta(9,t))+Ujf(8,t)-2*(Pijf(9,t)*Rijf(9)+Qijf(9,t)*Xijf(9))+Lijf(9,t)*(Rijf(9)^2+Xijf(9)^2);
        Ujf(10,t)>=-M*(1-beta(10,t))+Ujf(9,t)-2*(Pijf(10,t)*Rijf(10)+Qijf(10,t)*Xijf(10))+Lijf(10,t)*(Rijf(10)^2+Xijf(10)^2);
        Ujf(10,t)<=M*(1-beta(10,t))+Ujf(9,t)-2*(Pijf(10,t)*Rijf(10)+Qijf(10,t)*Xijf(10))+Lijf(10,t)*(Rijf(10)^2+Xijf(10)^2);
        Ujf(11,t)>=-M*(1-beta(11,t))+Ujf(10,t)-2*(Pijf(11,t)*Rijf(11)+Qijf(11,t)*Xijf(11))+Lijf(11,t)*(Rijf(11)^2+Xijf(11)^2);
        Ujf(11,t)<=M*(1-beta(11,t))+Ujf(10,t)-2*(Pijf(11,t)*Rijf(11)+Qijf(11,t)*Xijf(11))+Lijf(11,t)*(Rijf(11)^2+Xijf(11)^2);
        Ujf(12,t)>=-M*(1-beta(12,t))+Ujf(11,t)-2*(Pijf(12,t)*Rijf(12)+Qijf(12,t)*Xijf(12))+Lijf(12,t)*(Rijf(12)^2+Xijf(12)^2);
        Ujf(12,t)<=M*(1-beta(12,t))+Ujf(11,t)-2*(Pijf(12,t)*Rijf(12)+Qijf(12,t)*Xijf(12))+Lijf(12,t)*(Rijf(12)^2+Xijf(12)^2);
        Ujf(13,t)>=-M*(1-beta(13,t))+Ujf(12,t)-2*(Pijf(13,t)*Rijf(13)+Qijf(13,t)*Xijf(13))+Lijf(13,t)*(Rijf(13)^2+Xijf(13)^2);
        Ujf(13,t)<=M*(1-beta(13,t))+Ujf(12,t)-2*(Pijf(13,t)*Rijf(13)+Qijf(13,t)*Xijf(13))+Lijf(13,t)*(Rijf(13)^2+Xijf(13)^2);
        Ujf(14,t)>=-M*(1-beta(14,t))+Ujf(13,t)-2*(Pijf(14,t)*Rijf(14)+Qijf(14,t)*Xijf(14))+Lijf(14,t)*(Rijf(14)^2+Xijf(14)^2);
        Ujf(14,t)<=M*(1-beta(14,t))+Ujf(13,t)-2*(Pijf(14,t)*Rijf(14)+Qijf(14,t)*Xijf(14))+Lijf(14,t)*(Rijf(14)^2+Xijf(14)^2);
        Ujf(15,t)>=-M*(1-beta(15,t))+Ujf(14,t)-2*(Pijf(15,t)*Rijf(15)+Qijf(15,t)*Xijf(15))+Lijf(15,t)*(Rijf(15)^2+Xijf(15)^2);
        Ujf(15,t)<=M*(1-beta(15,t))+Ujf(14,t)-2*(Pijf(15,t)*Rijf(15)+Qijf(15,t)*Xijf(15))+Lijf(15,t)*(Rijf(15)^2+Xijf(15)^2);
        Ujf(16,t)>=-M*(1-beta(16,t))+Ujf(15,t)-2*(Pijf(16,t)*Rijf(16)+Qijf(16,t)*Xijf(16))+Lijf(16,t)*(Rijf(16)^2+Xijf(16)^2);
        Ujf(16,t)<=M*(1-beta(16,t))+Ujf(15,t)-2*(Pijf(16,t)*Rijf(16)+Qijf(16,t)*Xijf(16))+Lijf(16,t)*(Rijf(16)^2+Xijf(16)^2);
        Ujf(17,t)>=-M*(1-beta(17,t))+Ujf(16,t)-2*(Pijf(17,t)*Rijf(17)+Qijf(17,t)*Xijf(17))+Lijf(17,t)*(Rijf(17)^2+Xijf(17)^2);
        Ujf(17,t)<=M*(1-beta(17,t))+Ujf(16,t)-2*(Pijf(17,t)*Rijf(17)+Qijf(17,t)*Xijf(17))+Lijf(17,t)*(Rijf(17)^2+Xijf(17)^2);
        Ujf(18,t)>=-M*(1-beta(18,t))+Ujf(1,t)-2*(Pijf(18,t)*Rijf(18)+Qijf(18,t)*Xijf(18))+Lijf(18,t)*(Rijf(18)^2+Xijf(18)^2);
        Ujf(18,t)<=M*(1-beta(18,t))+Ujf(1,t)-2*(Pijf(18,t)*Rijf(18)+Qijf(18,t)*Xijf(18))+Lijf(18,t)*(Rijf(18)^2+Xijf(18)^2);
        Ujf(19,t)>=-M*(1-beta(19,t))+Ujf(18,t)-2*(Pijf(19,t)*Rijf(19)+Qijf(19,t)*Xijf(19))+Lijf(19,t)*(Rijf(19)^2+Xijf(19)^2);
        Ujf(19,t)<=M*(1-beta(19,t))+Ujf(18,t)-2*(Pijf(19,t)*Rijf(19)+Qijf(19,t)*Xijf(19))+Lijf(19,t)*(Rijf(19)^2+Xijf(19)^2);
        Ujf(20,t)>=-M*(1-beta(20,t))+Ujf(19,t)-2*(Pijf(20,t)*Rijf(20)+Qijf(20,t)*Xijf(20))+Lijf(20,t)*(Rijf(20)^2+Xijf(20)^2);
        Ujf(20,t)<=M*(1-beta(20,t))+Ujf(19,t)-2*(Pijf(20,t)*Rijf(20)+Qijf(20,t)*Xijf(20))+Lijf(20,t)*(Rijf(20)^2+Xijf(20)^2);
        Ujf(21,t)>=-M*(1-beta(21,t))+Ujf(20,t)-2*(Pijf(21,t)*Rijf(21)+Qijf(21,t)*Xijf(21))+Lijf(21,t)*(Rijf(21)^2+Xijf(21)^2);
        Ujf(21,t)<=M*(1-beta(21,t))+Ujf(20,t)-2*(Pijf(21,t)*Rijf(21)+Qijf(21,t)*Xijf(21))+Lijf(21,t)*(Rijf(21)^2+Xijf(21)^2);
        Ujf(22,t)>=-M*(1-beta(22,t))+Ujf(2,t)-2*(Pijf(22,t)*Rijf(22)+Qijf(22,t)*Xijf(22))+Lijf(22,t)*(Rijf(22)^2+Xijf(22)^2);
        Ujf(22,t)<=M*(1-beta(22,t))+Ujf(2,t)-2*(Pijf(22,t)*Rijf(22)+Qijf(22,t)*Xijf(22))+Lijf(22,t)*(Rijf(22)^2+Xijf(22)^2);
        Ujf(23,t)>=-M*(1-beta(23,t))+Ujf(22,t)-2*(Pijf(23,t)*Rijf(23)+Qijf(23,t)*Xijf(23))+Lijf(23,t)*(Rijf(23)^2+Xijf(23)^2);
        Ujf(23,t)<=M*(1-beta(23,t))+Ujf(22,t)-2*(Pijf(23,t)*Rijf(23)+Qijf(23,t)*Xijf(23))+Lijf(23,t)*(Rijf(23)^2+Xijf(23)^2);
        Ujf(24,t)>=-M*(1-beta(24,t))+Ujf(23,t)-2*(Pijf(24,t)*Rijf(24)+Qijf(24,t)*Xijf(24))+Lijf(24,t)*(Rijf(24)^2+Xijf(24)^2);
        Ujf(24,t)<=M*(1-beta(24,t))+Ujf(23,t)-2*(Pijf(24,t)*Rijf(24)+Qijf(24,t)*Xijf(24))+Lijf(24,t)*(Rijf(24)^2+Xijf(24)^2);
        Ujf(25,t)>=-M*(1-beta(25,t))+Ujf(5,t)-2*(Pijf(25,t)*Rijf(25)+Qijf(25,t)*Xijf(25))+Lijf(25,t)*(Rijf(25)^2+Xijf(25)^2);
        Ujf(25,t)<=M*(1-beta(25,t))+Ujf(5,t)-2*(Pijf(25,t)*Rijf(25)+Qijf(25,t)*Xijf(25))+Lijf(25,t)*(Rijf(25)^2+Xijf(25)^2);
        Ujf(26,t)>=-M*(1-beta(26,t))+Ujf(25,t)-2*(Pijf(26,t)*Rijf(26)+Qijf(26,t)*Xijf(26))+Lijf(26,t)*(Rijf(26)^2+Xijf(26)^2);
        Ujf(26,t)<=M*(1-beta(26,t))+Ujf(25,t)-2*(Pijf(26,t)*Rijf(26)+Qijf(26,t)*Xijf(26))+Lijf(26,t)*(Rijf(26)^2+Xijf(26)^2);
        Ujf(27,t)>=-M*(1-beta(27,t))+Ujf(26,t)-2*(Pijf(27,t)*Rijf(27)+Qijf(27,t)*Xijf(27))+Lijf(27,t)*(Rijf(27)^2+Xijf(27)^2);
        Ujf(27,t)<=M*(1-beta(27,t))+Ujf(26,t)-2*(Pijf(27,t)*Rijf(27)+Qijf(27,t)*Xijf(27))+Lijf(27,t)*(Rijf(27)^2+Xijf(27)^2);
        Ujf(28,t)>=-M*(1-beta(28,t))+Ujf(27,t)-2*(Pijf(28,t)*Rijf(28)+Qijf(28,t)*Xijf(28))+Lijf(28,t)*(Rijf(28)^2+Xijf(28)^2);
        Ujf(28,t)<=M*(1-beta(28,t))+Ujf(27,t)-2*(Pijf(28,t)*Rijf(28)+Qijf(28,t)*Xijf(28))+Lijf(28,t)*(Rijf(28)^2+Xijf(28)^2);
        Ujf(29,t)>=-M*(1-beta(29,t))+Ujf(28,t)-2*(Pijf(29,t)*Rijf(29)+Qijf(29,t)*Xijf(29))+Lijf(29,t)*(Rijf(29)^2+Xijf(29)^2);
        Ujf(29,t)<=M*(1-beta(29,t))+Ujf(28,t)-2*(Pijf(29,t)*Rijf(29)+Qijf(29,t)*Xijf(29))+Lijf(29,t)*(Rijf(29)^2+Xijf(29)^2);
        Ujf(30,t)>=-M*(1-beta(30,t))+Ujf(29,t)-2*(Pijf(30,t)*Rijf(30)+Qijf(30,t)*Xijf(30))+Lijf(30,t)*(Rijf(30)^2+Xijf(30)^2);
        Ujf(30,t)<=M*(1-beta(30,t))+Ujf(29,t)-2*(Pijf(30,t)*Rijf(30)+Qijf(30,t)*Xijf(30))+Lijf(30,t)*(Rijf(30)^2+Xijf(30)^2);
        Ujf(31,t)>=-M*(1-beta(31,t))+Ujf(30,t)-2*(Pijf(31,t)*Rijf(31)+Qijf(31,t)*Xijf(31))+Lijf(31,t)*(Rijf(31)^2+Xijf(31)^2);
        Ujf(31,t)<=M*(1-beta(31,t))+Ujf(30,t)-2*(Pijf(31,t)*Rijf(31)+Qijf(31,t)*Xijf(31))+Lijf(31,t)*(Rijf(31)^2+Xijf(31)^2);
        Ujf(32,t)>=-M*(1-beta(32,t))+Ujf(31,t)-2*(Pijf(32,t)*Rijf(32)+Qijf(32,t)*Xijf(32))+Lijf(32,t)*(Rijf(32)^2+Xijf(32)^2);
        Ujf(32,t)<=M*(1-beta(32,t))+Ujf(31,t)-2*(Pijf(32,t)*Rijf(32)+Qijf(32,t)*Xijf(32))+Lijf(32,t)*(Rijf(32)^2+Xijf(32)^2);
        
        Ujf(11,t)>=-M*(1-beta(33,t))+Ujf(21,t)-2*(Pijf(33,t)*Rijf(33)+Qijf(33,t)*Xijf(33))+Lijf(33,t)*(Rijf(33)^2+Xijf(33)^2);%
        Ujf(11,t)<=M*(1-beta(33,t))+Ujf(21,t)-2*(Pijf(33,t)*Rijf(33)+Qijf(33,t)*Xijf(33))+Lijf(33,t)*(Rijf(33)^2+Xijf(33)^2);%
        Ujf(7,t)>=-M*(1-beta(34,t))+Ujf(20,t)-2*(Pijf(34,t)*Rijf(34)+Qijf(34,t)*Xijf(34))+Lijf(34,t)*(Rijf(34)^2+Xijf(34)^2);%
        Ujf(7,t)<=M*(1-beta(34,t))+Ujf(20,t)-2*(Pijf(34,t)*Rijf(34)+Qijf(34,t)*Xijf(34))+Lijf(34,t)*(Rijf(34)^2+Xijf(34)^2);%
        Ujf(28,t)>=-M*(1-beta(35,t))+Ujf(24,t)-2*(Pijf(35,t)*Rijf(35)+Qijf(35,t)*Xijf(35))+Lijf(35,t)*(Rijf(35)^2+Xijf(35)^2);%
        Ujf(28,t)<=M*(1-beta(35,t))+Ujf(24,t)-2*(Pijf(35,t)*Rijf(35)+Qijf(35,t)*Xijf(35))+Lijf(35,t)*(Rijf(35)^2+Xijf(35)^2);%
        Ujf(17,t)>=-M*(1-beta(36,t))+Ujf(32,t)-2*(Pijf(36,t)*Rijf(36)+Qijf(36,t)*Xijf(36))+Lijf(36,t)*(Rijf(36)^2+Xijf(36)^2);%
        Ujf(17,t)<=M*(1-beta(36,t))+Ujf(32,t)-2*(Pijf(36,t)*Rijf(36)+Qijf(36,t)*Xijf(36))+Lijf(36,t)*(Rijf(36)^2+Xijf(36)^2);%
%1，2网向故障网送电的电压约束       
        Ujf(21,t)<=M*(1-(alphaIEEE_1(1,t)))+UjIEEE_1f(1,t)-2*(PIEEE_1(1,t)*RIEEE_1+QIEEE_1(1,t)*XIEEE_1)+LIEEE_1(1,t)*(RIEEE_1^2+XIEEE_1^2); %M后边乘的数是为支路ij的开关
        Ujf(21,t)>=-M*(1-(alphaIEEE_1(1,t)))+UjIEEE_1f(1,t)-2*(PIEEE_1(1,t)*RIEEE_1+QIEEE_1(1,t)*XIEEE_1)+LIEEE_1(1,t)*(RIEEE_1^2+XIEEE_1^2);  
%       Ujf(24)<=M*(1-betaIEEE_1)+UjIEEE_2-2*(PIEEE_2*RIEEE_2+QIEEE_2*XIEEE_2)+LIEEE_2*(RIEEE_2^2+XIEEE_2^2); 
%       Ujf(24)>=-M*(1-betaIEEE_1)+UjIEEE_2-2*(PIEEE_2*RIEEE_2+QIEEE_2*XIEEE_2)+LIEEE_2*(RIEEE_2^2+XIEEE_2^2); 

 %%对潮流加二阶锥松弛，利用公式S=UI推导至潮流公式以二阶锥（二范数）形式，并进行松弛，利用matlab中的二范数计算函数norm（A，2），若A里面是向量，2代表是二范数 二范数表示向量A的模，也就是A里元素的平方和开方   
        cone([2*Pijf(1,t);   2*Qijf(1,t);  Lijf(1,t)-UTf(1,t)],Lijf(1,t)+UTf(1,t)) ;%因为故障UTf(1,t)=0
        cone([2*Pijf(2,t);   2*Qijf(2,t);  Lijf(2,t)-Ujf(1,t)], Lijf(2,t)+Ujf(1,t)) ;
        cone([2*Pijf(3,t);   2*Qijf(3,t);  Lijf(3,t)-Ujf(2,t)], Lijf(3,t)+Ujf(2,t)) ;
        cone([2*Pijf(4,t);   2*Qijf(4,t);  Lijf(4,t)-Ujf(3,t)], Lijf(4,t)+Ujf(3,t)) ;
        cone([2*Pijf(5,t);   2*Qijf(5,t);  Lijf(5,t)-Ujf(4,t)], Lijf(5,t)+Ujf(4,t)) ;
        cone([2*Pijf(6,t);   2*Qijf(6,t);  Lijf(6,t)-Ujf(5,t)], Lijf(6,t)+Ujf(5,t)) ;
        cone([2*Pijf(7,t);   2*Qijf(7,t);  Lijf(7,t)-Ujf(6,t)], Lijf(7,t)+Ujf(6,t)) ;
        cone([2*Pijf(8,t);   2*Qijf(8,t);  Lijf(8,t)-Ujf(7,t)], Lijf(8,t)+Ujf(7,t)) ;
        cone([2*Pijf(9,t);   2*Qijf(9,t);  Lijf(9,t)-Ujf(8,t)], Lijf(9,t)+Ujf(8,t)) ;
        cone([2*Pijf(10,t);  2*Qijf(10,t);  Lijf(10,t)-Ujf(9,t)], Lijf(10,t)+Ujf(9,t)) ;
        cone([2*Pijf(11,t);  2*Qijf(11,t);  Lijf(11,t)-Ujf(10,t)],Lijf(11,t)+Ujf(10,t)) ;
        cone([2*Pijf(12,t);  2*Qijf(12,t);  Lijf(12,t)-Ujf(11,t)],Lijf(12,t)+Ujf(11,t)) ;
        cone([2*Pijf(13,t);  2*Qijf(13,t);  Lijf(13,t)-Ujf(12,t)],Lijf(13,t)+Ujf(12,t)) ;
        cone([2*Pijf(14,t);  2*Qijf(14,t);  Lijf(14,t)-Ujf(13,t)],Lijf(14,t)+Ujf(13,t)) ;
        cone([2*Pijf(15,t);  2*Qijf(15,t);  Lijf(15,t)-Ujf(14,t)],Lijf(15,t)+Ujf(14,t)) ;
        cone([2*Pijf(16,t);  2*Qijf(16,t);  Lijf(16,t)-Ujf(15,t)],Lijf(16,t)+Ujf(15,t)) ;
        cone([2*Pijf(17,t);  2*Qijf(17,t);  Lijf(17,t)-Ujf(16,t)],Lijf(17,t)+Ujf(16,t)) ;
        cone([2*Pijf(18,t);  2*Qijf(18,t);  Lijf(18,t)-Ujf(1,t)], Lijf(18,t)+Ujf(1,t)) ;
        cone([2*Pijf(19,t);  2*Qijf(19,t);  Lijf(19,t)-Ujf(18,t)],Lijf(19,t)+Ujf(18,t)) ;
        cone([2*Pijf(20,t);  2*Qijf(20,t);  Lijf(20,t)-Ujf(19,t)],Lijf(20,t)+Ujf(19,t)) ;
        cone([2*Pijf(21,t);  2*Qijf(21,t);  Lijf(21,t)-Ujf(20,t)],Lijf(21,t)+Ujf(20,t)) ;
        cone([2*Pijf(22,t);  2*Qijf(22,t);  Lijf(22,t)-Ujf(2,t)], Lijf(22,t)+Ujf(2,t)) ;
        cone([2*Pijf(23,t);  2*Qijf(23,t);  Lijf(23,t)-Ujf(22,t)],Lijf(23,t)+Ujf(22,t)) ;
        cone([2*Pijf(24,t);  2*Qijf(24,t);  Lijf(24,t)-Ujf(23,t)],Lijf(24,t)+Ujf(23,t)) ;
        cone([2*Pijf(25,t);  2*Qijf(25,t);  Lijf(25,t)-Ujf(5,t)], Lijf(25,t)+Ujf(5,t)) ;
        cone([2*Pijf(26,t);  2*Qijf(26,t);  Lijf(26,t)-Ujf(25,t)],Lijf(26,t)+Ujf(25,t)) ;
        cone([2*Pijf(27,t);  2*Qijf(27,t);  Lijf(27,t)-Ujf(26,t)],Lijf(27,t)+Ujf(26,t)) ;
        cone([2*Pijf(28,t);  2*Qijf(28,t);  Lijf(28,t)-Ujf(27,t)],Lijf(28,t)+Ujf(27,t)) ;
        cone([2*Pijf(29,t);  2*Qijf(29,t);  Lijf(29,t)-Ujf(28,t)],Lijf(29,t)+Ujf(28,t)) ;
        cone([2*Pijf(30,t);  2*Qijf(30,t);  Lijf(30,t)-Ujf(29,t)],Lijf(30,t)+Ujf(29,t)) ;
        cone([2*Pijf(31,t);  2*Qijf(31,t);  Lijf(31,t)-Ujf(30,t)],Lijf(31,t)+Ujf(30,t)) ;
        cone([2*Pijf(32,t);  2*Qijf(32,t);  Lijf(32,t)-Ujf(31,t)],Lijf(32,t)+Ujf(31,t)) ;
        
        cone([2*Pijf(33,t);  2*Qijf(33,t);  Lijf(33,t)-Ujf(21,t)],Lijf(33,t)+Ujf(21,t)) ;%
        cone([2*Pijf(34,t);  2*Qijf(34,t);  Lijf(34,t)-Ujf(20,t)],Lijf(34,t)+Ujf(20,t)) ;%
        cone([2*Pijf(35,t);  2*Qijf(35,t);  Lijf(35,t)-Ujf(24,t)],Lijf(35,t)+Ujf(24,t)) ;%
        cone([2*Pijf(36,t);  2*Qijf(36,t);  Lijf(36,t)-Ujf(32,t)],Lijf(36,t)+Ujf(32,t)) ;%
        %联络线的二阶锥松弛
        cone([2*PIEEE_1(1,t); 2*QIEEE_1(1,t); LIEEE_1(1,t)-UjIEEE_1f(1,t)],LIEEE_1(1,t)+UjIEEE_1f(1,t));

];

%联络线约束，数环保证配电网辐射状结构
F=[F;
    sum(beta(8:11,t))+beta(21,t)+sum(beta(33:34,t))<=6;%一个供电区域内允许闭合的开关总数
    sum(beta(2:7,t))+sum(beta(18:20,t))+beta(34,t)<=9;
    sum(beta(3:5,t))+sum(beta(22:28,t))+beta(35,t)<=10;
    sum(beta(6:17,t))+sum(beta(25:32,t))+beta(36,t)<=20;
    
    sum(beta(2:11,t))+sum(beta(18:21,t))+beta(33,t)<=14;
    sum(beta(6:7,t))+sum(beta(12:17,t))+beta(21,t)+sum(beta(25:34,t))+beta(36,t)<=19;
    beta(2,t)+sum(beta(6:7,t))+sum(beta(18:20,t))+sum(beta(22:28,t))+sum(beta(34:35,t))<=14;
    sum(beta(2:5,t))+sum(beta(8:20,t))+sum(beta(25:32,t))+beta(34,t)+beta(36,t)<=26;
    sum(beta(3:17,t))+sum(beta(22:24,t))+sum(beta(29:32,t))+beta(35,t)+beta(36,t)<=23;
    
    beta(2,t)+sum(beta(6:11,t))+sum(beta(18:28,t))+beta(33,t)+beta(35,t)<=19;
    sum(beta(2:5,t))+sum(beta(12:21,t))+sum(beta(25:33,t))+beta(36,t)<=23;
    sum(beta(3:7,t))+sum(beta(12:17,t))+sum(beta(21:24,t))+sum(beta(29:36,t))<=22;
    beta(2,t)+sum(beta(8:20,t))+sum(beta(22:24,t))+sum(beta(29:32,t))+sum(beta(34:36,t))<=23;
    
    beta(2,t)+sum(beta(12:24,t))+sum(beta(29:33,t))+beta(35,t)+beta(36,t)<=20;
    
    ];    
 end
  end
 %% 调用求解器计算潮流
        

        options=sdpsettings('solver','gurobi');
        
        optimize(F,  f2,  options);
        
        TF=strcmp(ans.info,'Successfully solved (GUROBI)');
        
        % 检查收敛性
    prim_res = norm(Ess_var - Ess_var_old, 'fro');
    fprintf('ADMM iter %d: primal_residual = %.3e\n', inner_iter, prim_res);
    if prim_res < tol
        disp('ADMM Converged');
      
    end
end
   Unode2=sqrt(value(Ujf)); %f2=value(ff);
  Pijnode=value(Pijf); Qijnode=value(Qijf);  
  local_P2=value(PIEEE_1f);local_Q2=value(QIEEE_1f);local_U2=value(UjIEEE_1f);
  Ess2 = value(Ess); P_chf = value(Pchf); P_disf = value(Pdisf);