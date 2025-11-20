function [Unode1, f1_val,local_P1,local_Q1,local_U1,Pch1,Pdis1,z_soc1_hist]=subnormal1(zP,zQ,zU,uP,uQ,uU,z_soc1,lam1)
max_iter = 100;
max_iiter = 50;
T = 2;               % 时间步长
rho_s = 0.01;
rho_t = 0.001;
tol =  1e-3;
%变压器参数设置 
UTmax=1.05^2;  %电压运行上限
UTmin=0.95^2;    
Tran_Num=1; %有载调压变数量
k0=1;  %档位为零时，OLTC高低压侧的变比
Klimit1=4;  %档位上下限
Kt1=sdpvar(1,5);  %有载调压变压器的档位
Kt_1=sdpvar(1,6);  %OLTC上一时刻档位数
bk_1=binvar(2*Klimit1+1,5); 
Uc_1=sdpvar(2*Klimit1+1,5);
Ktp1=sdpvar(1,5);  Ktn1=sdpvar(1,5);
deltk=0.025;  %OLTC档位变比的单步调节步长

%% 多时段变量定义
%正常配电网
bus_Num=32;
L=32;
Sijmax=10;
Lijmax=10;
Umax=1.05^2;
Umin=0.85^2;
Ubase=1.03^2;
PV_Num1=5; 
Ess_Num1=5;
%% 光伏参数
Spv=[0.01;0.02;0.03;0.02;0.015]; 
f_pv_bian=[1,1,1,1,1,1];   %光伏变化T=6
PV1=Spv*f_pv_bian;  %光伏出力


%% 节点负荷实时数据
%配电网
Pnet= [0.12 ;0.108 ;0.144 ;0.072 ;0.072 ;0.24 ;0.24 ;0.072 ;0.072 ;0.054 ;0.072 ;0.072 ;0.144 ;0.072 ;0.072 ;0.072 ;0.108 ;0.108 ;
0.108 ;0.108 ;0.108 ;0.108 ;0.504 ;0.504 ;0.072 ;0.072 ;0.072 ;0.144 ;0.24 ;0.18 ;0.252 ;0.072 ;]/10;
Qnet= [0.072 ;0.048 ;0.096 ;0.036 ;0.024 ;0.12 ;0.12 ;0.024 ;0.024 ;0.036 ;0.042 ;0.042 ;0.096 ;0.012 ;0.024 ;0.024 ;0.048 ;
0.048 ;0.048 ;0.048 ;0.048 ;0.06 ;0.24 ;0.24 ;0.03 ;0.03 ;0.024 ;0.084 ;0.72 ;0.084 ;0.12 ;0.048 ;]/10;

for t=1:T
    for n=1:32
      P(n,t)=Pnet(n);
      Q(n,t)=Qnet(n);
     end
end
Pj=P;  %单位10MW
Qj=Q;
%% 线路阻抗参数
%配电网
Ldata=[0.0922,0.047;  0.493,0.2511;  0.366,0.1864;  0.3811,0.1941;  0.0819,0.707;  0.1872,0.6188;  0.7114,0.2351;  1.03,0.74;       1.044,0.74;     0.1966,0.065;  0.3744,0.1238;
                 1.468,1.155;  0.5416,0.7129; 0.591,0.526;   0.7463,0.545;   1.289,1.721;   0.372,0.574;    0.164,0.1565;   1.5042,1.3554;   0.4095,0.4784;  0.7089,0.9373; 0.4512, 0.3083;
                 0.898,0.7091; 0.896,0.7011;  0.203,0.1034;  0.2842,0.1447;  1.059,0.9337;  0.8042,0.7006;  0.5075,0.2585;  0.9744,0.963;    0.3105, 0.3619; 0.341, 0.5362;       
                   ];

Rij=Ldata(1:L,1)/16.02756;
Xij=Ldata(1:L,2)/16.02756;
%% 定义优化变量及其范围
% 基本线路节点参数
%正常配电网1
Pij1=sdpvar(32,T);
Qij1=sdpvar(32,T);
Uj1=sdpvar(32,T);
UT1=sdpvar(1,T);%变压器二次侧电压
Lij1=sdpvar(32,T); 
Pdec1=sdpvar(5,T);
Qdec1=sdpvar(5,T);
subp1=sdpvar(1,T);
subn1=sdpvar(1,T);
P_dis1=sdpvar(5,T);
P_ch1=sdpvar(5,T);
u_ch1=sdpvar(Ess_Num1,T);%充电状态
u_dis1=sdpvar(Ess_Num1,T); %放电状态
Ess1=sdpvar(Ess_Num1,T+1);
PIEEE_1=sdpvar(1,T);
QIEEE_1=sdpvar(1,T);
UjIEEE_1=sdpvar(1,T);

%储能充放电参数
P_essmax=2*1.5/10; % 充放最大功率  基值是10MW
z_soc1 = 0.6*ones(Ess_Num1, T+1);%全局soc变量
Ess1_local = 0.6*ones(Ess_Num1, T+1);%局部soc变量
lam1 = zeros(Ess_Num1, T+1); %soc拉格朗日乘子
z_soc1_hist = zeros(Ess_Num1, T+1, max_iiter);
for k=1:max_iiter
    % fprintf('===ADMM step %d===\n',k);

%% 优化目标
for t=1:T
   % Ess1(:,t)=z_soc1(:,t); Ess1(:,t+1)=z_soc1(:,t+1);
Cl=5000;Cp=4000;Cpv=7000;Cb=0.5;%有功缩减成本
Ce=4000;
  R=repmat(Rij,1,T);
  lost1=Lij1.*R;
 

   f_l1=Cp*(sum(lost1(:,t)));
  f_pv1=Cpv*(sum(Pdec1(:,t)));
   f_ess1=Ce*sum(P_ch1(:,t)) + Ce*sum(P_dis1(:,t));
   f_a= 0.5*rho_s*(PIEEE_1(1,t)-zP(t)+uP(t))^2+0.5*rho_s*(QIEEE_1(1,t)-zQ(t)+uQ(t))^2+0.5*rho_s*(UjIEEE_1(1,t)-zU(t)+uU(t))^2;
    % f_c=0.5*rho_t *(Ess - z_soc +lam1).^2 + 0.5*rho_t *(Ess - z_soc + lam1).^2;
    % f_c=0.5*rho_t *(r(:,t+1) - r(:,t) +lam1).^2 + 0.5*rho_t *(r(:,t+1) - r(:,t) + lam1).^2;
    vec= Ess1(:,t:t+1) - z_soc1(:,t:t+1) +lam1(:,t:t+1); %Ess_Num1×2
    f_c=0.5*rho_t* sum(vec(:).^2);
%   f_l=Cp*(sum(sum(lost(:,1:T))));
   f1 =  f_l1+ f_pv1 + f_ess1 + f_a + f_c;

  % f1 = f_a + f_c;

end



%% 约束条件
%正常配电网
F=[];
%光伏约束
for t=1:T
      F=[F;     
          
        0<=Pdec1(:,t);     Pdec1(:,t)<=PV1(:,t);
%         PV1(k,t) == PVmax(k,t) - Pdec1(k,t);% PV实际出力
        -tan(acos(0.9)) * PV1(:,t) <= Qdec1(:,t); Qdec1(:,t)<= tan(acos(0.9)) * PV1(:,t);
        ];
     
     
 
end
%储能约束
      F = [F;
          
          u_dis1(:)+u_ch1(:)<=1;       u_dis1(:)+u_ch1(:)>=0;%表示充电，放电，不充不放三种状态
          P_dis1(:,t)>=0;P_dis1(:,t)<=u_dis1(:,t)*P_essmax;%储能放电功率约束
          P_ch1(:,t)>=0;P_ch1(:,t)<=u_ch1(:,t)*P_essmax;%储能放电功率约束
          Ess1(:,1)==0.6;     
          Ess1(:,t)>=0.2; Ess1(:,t)<=0.9;
          Ess1(:,T+1)>=0.4; Ess1(:,T+1)<=0.6; 
          % Ess1(:,t) == Ess_guess(:,t);  % fix SOC 
          % Ess1(:,t+1) == Ess_guess(:,t+1);
           Ess1_local(:,t+1)==z_soc1(:,t)+0.9*P_ch1(:,t)-1.11*P_dis1(:,t); 


    ];  
%多时段程序
for t=1:T

F=[F;
        0<=Lij1(:,t);      Lij1(:,t)<=Lijmax;  %线路电流约束
        Umin<=Uj1(:,t);    Uj1(:,t)<=Umax;   %节点电压约束   
        UTmin<=UT1(:,t);    UT1(:,t)<=UTmax;
%       subp(1,t)-subn(1,t)==Pij1(1,t)/12
        subp1(1,t)>=0;
        subn1(1,t)>=0;   
        Uj1(24,t) == UjIEEE_1(1,t);
%正常配电网多时间尺度各段线路的潮流等式方程
        Pij1(1,t)- Lij1(1,t)*Rij(1)-Pj(1,t)==Pij1(2,t)+Pij1(18,t);
        Pij1(2,t)- Lij1(2,t)*Rij(2)-Pj(2,t)==Pij1(3,t)+Pij1(22,t);
        Pij1(3,t)- Lij1(3,t)*Rij(3)-Pj(3,t)==Pij1(4,t);
        Pij1(4,t)- Lij1(4,t)*Rij(4)-Pj(4,t)==Pij1(5,t);
        Pij1(5,t)- Lij1(5,t)*Rij(5)-Pj(5,t)==Pij1(6,t)+Pij1(25,t);
        Pij1(6,t)- Lij1(6,t)*Rij(6)-Pj(6,t)==Pij1(7,t);
        Pij1(7,t)- Lij1(7,t)*Rij(7)-Pj(7,t)+P_dis1(1,t)-P_ch1(1,t)==Pij1(8,t);%7节点有储能
        % Ess1(1,t+1)==Ess1(1,t)+0.9*P_ch1(1,t)-1.11*P_dis1(1,t); 
        
        Pij1(8,t)- Lij1(8,t)*Rij(8)-Pj(8,t)+PV1(1,t)-Pdec1(1,t)==Pij1(9,t);  %8节点有光伏
        Pij1(9,t)- Lij1(9,t)*Rij(9)-Pj(9,t)==Pij1(10,t);
        Pij1(10,t)-Lij1(10,t)*Rij(10)-Pj(10,t)+P_dis1(2,t)-P_ch1(2,t)==Pij1(11,t);%10有储能
        % Ess1(2,t+1)==Ess1(2,t)+0.9*P_ch1(2,t)-1.11*P_dis1(2,t); 

        Pij1(11,t)- Lij1(11,t)*Rij(11)-Pj(11,t)==Pij1(12,t);
        Pij1(12,t)- Lij1(12,t)*Rij(12)-Pj(12,t)+PV1(2,t)-Pdec1(2,t)==Pij1(13,t);%12节点有光伏
        Pij1(13,t)- Lij1(13,t)*Rij(13)-Pj(13,t)==Pij1(14,t);
        Pij1(14,t)- Lij1(14,t)*Rij(14)-Pj(14,t)==Pij1(15,t);
        Pij1(15,t)- Lij1(15,t)*Rij(15)-Pj(15,t)==Pij1(16,t);
        Pij1(16,t)- Lij1(16,t)*Rij(16)-Pj(16,t)==Pij1(17,t);
        Pij1(17,t)- Lij1(17,t)*Rij(17)-Pj(17,t)==0;
        
        Pij1(18,t)- Lij1(18,t)*Rij(18)-Pj(18,t)==Pij1(19,t);
        Pij1(19,t)- Lij1(19,t)*Rij(19)-Pj(19,t)+PV1(3,t)-Pdec1(3,t)==Pij1(20,t);%19节点有光伏
        Pij1(20,t)- Lij1(20,t)*Rij(20)-Pj(20,t)+P_dis1(3,t)-P_ch1(3,t)==Pij1(21,t);%20有储能
        % Ess1(3,t+1)==Ess1(3,t)+0.9*P_ch1(3,t)-1.11*P_dis1(3,t); 

        Pij1(21,t)- Lij1(21,t)*Rij(21)-Pj(21,t)==0;
        Pij1(22,t)- Lij1(22,t)*Rij(22)-Pj(22,t)==Pij1(23,t);
        Pij1(23,t)- Lij1(23,t)*Rij(23)-Pj(23,t)+PV1(4,t)-Pdec1(4,t)==Pij1(24,t);%23节点有光伏
        Pij1(24,t)- Lij1(24,t)*Rij(24)-Pj(24,t)==PIEEE_1(1,t);
       
        
        Pij1(25,t)- Lij1(25,t)*Rij(25)-Pj(25,t)==Pij1(26,t);
        Pij1(26,t)- Lij1(26,t)*Rij(26)-Pj(26,t)==Pij1(27,t);
        Pij1(27,t)- Lij1(27,t)*Rij(27)-Pj(27,t)==Pij1(28,t); 
        Pij1(28,t)- Lij1(28,t)*Rij(28)-Pj(28,t)+P_dis1(4,t)-P_ch1(4,t)==Pij1(29,t);%28节点有储能
        % Ess1(4,t+1)==Ess1(4,t)+0.9*P_ch1(4,t)-1.11*P_dis1(4,t); 

        Pij1(29,t)- Lij1(29,t)*Rij(29)-Pj(29,t)==Pij1(30,t);
        Pij1(30,t)- Lij1(30,t)*Rij(30)-Pj(30,t)==Pij1(31,t);
        Pij1(31,t)- Lij1(31,t)*Rij(31)-Pj(31,t)+PV1(5,t)-Pdec1(5,t)+P_dis1(5,t)-P_ch1(5,t)==Pij1(32,t);%31节点有光伏、储能
        % Ess1(5,t+1)==Ess1(5,t)+0.9*P_ch1(5,t)-1.11*P_dis1(5,t); 
        Pij1(32,t)- Lij1(32,t)*Rij(32)-Pj(32,t)==0;
        
        Qij1(1,t) - Lij1(1,t)*Xij(1)-  Qj(1,t) == Qij1(2,t) + Qij1(18,t);  
        Qij1(2,t) - Lij1(2,t)*Xij(2) - Qj(2,t) == Qij1(3,t) + Qij1(22,t); 
        Qij1(3,t) - Lij1(3,t)*Xij(3) - Qj(3,t) == Qij1(4,t);            
        Qij1(4,t) - Lij1(4,t)*Xij(4) - Qj(4,t) == Qij1(5,t);             
        Qij1(5,t) - Lij1(5,t)*Xij(5) - Qj(5,t) == Qij1(6,t) + Qij1(25,t);
        Qij1(6,t) - Lij1(6,t)*Xij(6) - Qj(6,t) == Qij1(7,t);             
        Qij1(7,t) - Lij1(7,t)*Xij(7) - Qj(7,t) == Qij1(8,t); % 节点7(ESS1)
        
        Qij1(8,t) - Lij1(8,t)*Xij(8) - Qj(8,t) + Qdec1(1,t) == Qij1(9,t); % 节点8(PV1)
        Qdec1(1,t) <= (PV1(1,t) - Pdec1(1,t)) * 0.32868421; % 功率因数约束
        Qdec1(1,t) >= ( Pdec1(1,t)-PV1(1,t) ) * 0.32868421;
        
        Qij1(9,t) - Lij1(9,t)*Xij(9)   - Qj(9,t) == Qij1(10,t);            
        Qij1(10,t)- Lij1(10,t)*Xij(10) - Qj(10,t) == Qij1(11,t);          
        Qij1(11,t)- Lij1(11,t)*Xij(11) - Qj(11,t) == Qij1(12,t);        
        Qij1(12,t)- Lij1(12,t)*Xij(12) - Qj(12,t) + Qdec1(2,t) == Qij1(13,t); % 节点12(PV2)
        Qdec1(2,t) <= (PV1(2,t) - Pdec1(2,t)) * 0.32868421;
        Qdec1(2,t) >= (Pdec1(2,t)- PV1(2,t) ) * 0.32868421;

        Qij1(13,t) - Lij1(13,t)*Xij(13) - Qj(13,t) == Qij1(14,t);       
        Qij1(14,t) - Lij1(14,t)*Xij(14) - Qj(14,t) == Qij1(15,t);        
        Qij1(15,t) - Lij1(15,t)*Xij(15) - Qj(15,t) == Qij1(16,t);        
        Qij1(16,t) - Lij1(16,t)*Xij(16) - Qj(16,t) == Qij1(17,t);       
        Qij1(17,t) - Lij1(17,t)*Xij(17) - Qj(17,t) == 0;  
        
        Qij1(18,t) - Lij1(18,t)*Xij(18) - Qj(18,t) == Qij1(19,t);       
        Qij1(19,t) - Lij1(19,t)*Xij(19) - Qj(19,t) + Qdec1(3,t) == Qij1(20,t); % 节点19(PV3)
        Qdec1(3,t) <= (PV1(3,t) - Pdec1(3,t)) * 0.32868421;
        Qdec1(3,t) >= ( Pdec1(3,t)- PV1(3,t)) * 0.32868421;

        Qij1(20,t) - Lij1(20,t)*Xij(20) - Qj(20,t) == Qij1(21,t);     
        Qij1(21,t) - Lij1(21,t)*Xij(21) - Qj(21,t) == 0;                 
        
        Qij1(22,t) - Lij1(22,t)*Xij(22) - Qj(22,t) == Qij1(23,t);        
        Qij1(23,t) - Lij1(23,t)*Xij(23) - Qj(23,t) + Qdec1(4,t) == Qij1(24,t); % 节点23(PV4)
        Qdec1(4,t) <= (PV1(4,t) - Pdec1(4,t)) * 0.32868421;
        Qdec1(4,t) >= (Pdec1(4,t)- PV1(4,t) ) * 0.32868421;

        Qij1(24,t) - Lij1(24,t)*Xij(24) - Qj(24,t) == QIEEE_1(1,t);
        Qij1(25,t) - Lij1(25,t)*Xij(25) - Qj(25,t) == Qij1(26,t);    
        Qij1(26,t) - Lij1(26,t)*Xij(26) - Qj(26,t) == Qij1(27,t);       
        Qij1(27,t) - Lij1(27,t)*Xij(27) - Qj(27,t) == Qij1(28,t);        
        Qij1(28,t) - Lij1(28,t)*Xij(28) - Qj(28,t) == Qij1(29,t); % 节点28(ESS2)
        Qij1(29,t) - Lij1(29,t)*Xij(29) - Qj(29,t) == Qij1(30,t);       
        Qij1(30,t) - Lij1(30,t)*Xij(30) - Qj(30,t) == Qij1(31,t);        
        Qij1(31,t) - Lij1(31,t)*Xij(31) - Qj(31,t) + Qdec1(5,t) == Qij1(32,t); % 节点31(PV5)
        Qdec1(5,t) <= (PV1(5,t) - Pdec1(5,t)) * 0.32868421;
        Qdec1(5,t) >= (Pdec1(5,t)- PV1(5,t)) * 0.32868421;
        Qij1(32,t)-Lij1(32,t)*Xij(32)-Qj(32,t)==0;
        
        Uj1(1,t) == UT1(1,t) - 2*(Pij1(1,t)*Rij(1) + Qij1(1,t)*Xij(1)) + Lij1(1,t)*(Rij(1)^2 + Xij(1)^2);
        Uj1(2,t) == Uj1(1,t) - 2*(Pij1(2,t)*Rij(2) + Qij1(2,t)*Xij(2)) + Lij1(2,t)*(Rij(2)^2 + Xij(2)^2);
        Uj1(3,t) == Uj1(2,t) - 2*(Pij1(3,t)*Rij(3) + Qij1(3,t)*Xij(3)) + Lij1(3,t)*(Rij(3)^2 + Xij(3)^2);
        Uj1(4,t) == Uj1(3,t) - 2*(Pij1(4,t)*Rij(4) + Qij1(4,t)*Xij(4)) + Lij1(4,t)*(Rij(4)^2 + Xij(4)^2);
        Uj1(5,t) == Uj1(4,t) - 2*(Pij1(5,t)*Rij(5) + Qij1(5,t)*Xij(5)) + Lij1(5,t)*(Rij(5)^2 + Xij(5)^2);
        Uj1(6,t) == Uj1(5,t) - 2*(Pij1(6,t)*Rij(6) + Qij1(6,t)*Xij(6)) + Lij1(6,t)*(Rij(6)^2 + Xij(6)^2);
        Uj1(7,t) == Uj1(6,t) - 2*(Pij1(7,t)*Rij(7) + Qij1(7,t)*Xij(7)) + Lij1(7,t)*(Rij(7)^2 + Xij(7)^2);
        Uj1(8,t) == Uj1(7,t) - 2*(Pij1(8,t)*Rij(8) + Qij1(8,t)*Xij(8)) + Lij1(8,t)*(Rij(8)^2 + Xij(8)^2);
        Uj1(9,t) == Uj1(8,t) - 2*(Pij1(9,t)*Rij(9) + Qij1(9,t)*Xij(9)) + Lij1(9,t)*(Rij(9)^2 + Xij(9)^2);
        Uj1(10,t) == Uj1(9,t) - 2*(Pij1(10,t)*Rij(10) + Qij1(10,t)*Xij(10)) + Lij1(10,t)*(Rij(10)^2 + Xij(10)^2);
        Uj1(11,t) == Uj1(10,t) - 2*(Pij1(11,t)*Rij(11) + Qij1(11,t)*Xij(11)) + Lij1(11,t)*(Rij(11)^2 + Xij(11)^2);
        Uj1(12,t) == Uj1(11,t) - 2*(Pij1(12,t)*Rij(12) + Qij1(12,t)*Xij(12)) + Lij1(12,t)*(Rij(12)^2 + Xij(12)^2);
        Uj1(13,t) == Uj1(12,t) - 2*(Pij1(13,t)*Rij(13) + Qij1(13,t)*Xij(13)) + Lij1(13,t)*(Rij(13)^2 + Xij(13)^2);
        Uj1(14,t) == Uj1(13,t) - 2*(Pij1(14,t)*Rij(14) + Qij1(14,t)*Xij(14)) + Lij1(14,t)*(Rij(14)^2 + Xij(14)^2);
        Uj1(15,t) == Uj1(14,t) - 2*(Pij1(15,t)*Rij(15) + Qij1(15,t)*Xij(15)) + Lij1(15,t)*(Rij(15)^2 + Xij(15)^2);
        Uj1(16,t) == Uj1(15,t) - 2*(Pij1(16,t)*Rij(16) + Qij1(16,t)*Xij(16)) + Lij1(16,t)*(Rij(16)^2 + Xij(16)^2);
        Uj1(17,t) == Uj1(16,t) - 2*(Pij1(17,t)*Rij(17) + Qij1(17,t)*Xij(17)) + Lij1(17,t)*(Rij(17)^2 + Xij(17)^2);
        Uj1(18,t) == Uj1(1,t)  - 2*(Pij1(18,t)*Rij(18) + Qij1(18,t)*Xij(18)) + Lij1(18,t)*(Rij(18)^2 + Xij(18)^2);
        Uj1(19,t) == Uj1(18,t) - 2*(Pij1(19,t)*Rij(19) + Qij1(19,t)*Xij(19)) + Lij1(19,t)*(Rij(19)^2 + Xij(19)^2);
        Uj1(20,t) == Uj1(19,t) - 2*(Pij1(20,t)*Rij(20) + Qij1(20,t)*Xij(20)) + Lij1(20,t)*(Rij(20)^2 + Xij(20)^2);
        Uj1(21,t) == Uj1(20,t) - 2*(Pij1(21,t)*Rij(21) + Qij1(21,t)*Xij(21)) + Lij1(21,t)*(Rij(21)^2 + Xij(21)^2);
        Uj1(22,t) == Uj1(2,t)  - 2*(Pij1(22,t)*Rij(22) + Qij1(22,t)*Xij(22)) + Lij1(22,t)*(Rij(22)^2 + Xij(22)^2);
        Uj1(23,t) == Uj1(22,t) - 2*(Pij1(23,t)*Rij(23) + Qij1(23,t)*Xij(23)) + Lij1(23,t)*(Rij(23)^2 + Xij(23)^2);
        Uj1(24,t) == Uj1(23,t) - 2*(Pij1(24,t)*Rij(24) + Qij1(24,t)*Xij(24)) + Lij1(24,t)*(Rij(24)^2 + Xij(24)^2);
        Uj1(25,t) == Uj1(5,t)  - 2*(Pij1(25,t)*Rij(25) + Qij1(25,t)*Xij(25)) + Lij1(25,t)*(Rij(25)^2 + Xij(25)^2);
        Uj1(26,t) == Uj1(25,t) - 2*(Pij1(26,t)*Rij(26) + Qij1(26,t)*Xij(26)) + Lij1(26,t)*(Rij(26)^2 + Xij(26)^2);
        Uj1(27,t) == Uj1(26,t) - 2*(Pij1(27,t)*Rij(27) + Qij1(27,t)*Xij(27)) + Lij1(27,t)*(Rij(27)^2 + Xij(27)^2);
        Uj1(28,t) == Uj1(27,t) - 2*(Pij1(28,t)*Rij(28) + Qij1(28,t)*Xij(28)) + Lij1(28,t)*(Rij(28)^2 + Xij(28)^2);
        Uj1(29,t) == Uj1(28,t) - 2*(Pij1(29,t)*Rij(29) + Qij1(29,t)*Xij(29)) + Lij1(29,t)*(Rij(29)^2 + Xij(29)^2);
        Uj1(30,t) == Uj1(29,t) - 2*(Pij1(30,t)*Rij(30) + Qij1(30,t)*Xij(30)) + Lij1(30,t)*(Rij(30)^2 + Xij(30)^2);
        Uj1(31,t) == Uj1(30,t) - 2*(Pij1(31,t)*Rij(31) + Qij1(31,t)*Xij(31)) + Lij1(31,t)*(Rij(31)^2 + Xij(31)^2);
        Uj1(32,t) == Uj1(31,t) - 2*(Pij1(32,t)*Rij(32) + Qij1(32,t)*Xij(32)) + Lij1(32,t)*(Rij(32)^2 + Xij(32)^2);  
        
%对潮流加二阶锥松弛，利用公式S=UI推导至潮流公式以二阶锥（二范数）形式，并进行松弛，利用matlab中的二范数计算函数norm（A，2），若A里面是向量，2代表是二范数 二范数表示向量A的模，也就是A里元素的平方和开方  
    
        cone([2*Pij1(1,t);   2*Qij1(1,t);  Lij1(1,t)-UT1(1,t)],Lij1(1,t)+UT1(1,t));
        cone([2*Pij1(2,t);   2*Qij1(2,t);  Lij1(2,t)-Uj1(1,t)],Lij1(2,t)+Uj1(1,t)) ;
        cone([2*Pij1(3,t);   2*Qij1(3,t);  Lij1(3,t)-Uj1(2,t)],Lij1(3,t)+Uj1(2,t)) ;
        cone([2*Pij1(4,t);   2*Qij1(4,t);  Lij1(4,t)-Uj1(3,t)],Lij1(4,t)+Uj1(3,t)) ;
        cone([2*Pij1(5,t);   2*Qij1(5,t);  Lij1(5,t)-Uj1(4,t)],Lij1(5,t)+Uj1(4,t)) ;
        cone([2*Pij1(6,t);   2*Qij1(6,t);  Lij1(6,t)-Uj1(5,t)],Lij1(6,t)+Uj1(5,t)) ;
        cone([2*Pij1(7,t);   2*Qij1(7,t);  Lij1(7,t)-Uj1(6,t)],Lij1(7,t)+Uj1(6,t)) ;
        cone([2*Pij1(8,t);   2*Qij1(8,t);  Lij1(8,t)-Uj1(7,t)],Lij1(8,t)+Uj1(7,t)) ;
        cone([2*Pij1(9,t);   2*Qij1(9,t);  Lij1(9,t)-Uj1(8,t)],Lij1(9,t)+Uj1(8,t)) ;
        cone([2*Pij1(10,t);  2*Qij1(10,t);  Lij1(10,t)-Uj1(9,t)], Lij1(10,t)+Uj1(9,t)) ;
        cone([2*Pij1(11,t);  2*Qij1(11,t);  Lij1(11,t)-Uj1(10,t)],Lij1(11,t)+Uj1(10,t));
        cone([2*Pij1(12,t);  2*Qij1(12,t);  Lij1(12,t)-Uj1(11,t)],Lij1(12,t)+Uj1(11,t));
        cone([2*Pij1(13,t);  2*Qij1(13,t);  Lij1(13,t)-Uj1(12,t)],Lij1(13,t)+Uj1(12,t)) ;
        cone([2*Pij1(14,t);  2*Qij1(14,t);  Lij1(14,t)-Uj1(13,t)],Lij1(14,t)+Uj1(13,t)) ;
        cone([2*Pij1(15,t);  2*Qij1(15,t);  Lij1(15,t)-Uj1(14,t)],Lij1(15,t)+Uj1(14,t)) ;
        cone([2*Pij1(16,t);  2*Qij1(16,t);  Lij1(16,t)-Uj1(15,t)],Lij1(16,t)+Uj1(15,t)) ;
        cone([2*Pij1(17,t);  2*Qij1(17,t);  Lij1(17,t)-Uj1(16,t)],Lij1(17,t)+Uj1(16,t)) ;
        cone([2*Pij1(18,t);  2*Qij1(18,t);  Lij1(18,t)-Uj1(1,t)], Lij1(18,t)+Uj1(1,t)) ;
        cone([2*Pij1(19,t);  2*Qij1(19,t);  Lij1(19,t)-Uj1(18,t)],Lij1(19,t)+Uj1(18,t));
        cone([2*Pij1(20,t);  2*Qij1(20,t);  Lij1(20,t)-Uj1(19,t)],Lij1(20,t)+Uj1(19,t));
        cone([2*Pij1(21,t);  2*Qij1(21,t);  Lij1(21,t)-Uj1(20,t)],Lij1(21,t)+Uj1(20,t)) ;
        cone([2*Pij1(22,t);  2*Qij1(22,t);  Lij1(22,t)-Uj1(2,t)], Lij1(22,t)+Uj1(2,t)) ;
        cone([2*Pij1(23,t);  2*Qij1(23,t);  Lij1(23,t)-Uj1(22,t)],Lij1(23,t)+Uj1(22,t)) ;
        cone([2*Pij1(24,t);  2*Qij1(24,t);  Lij1(24,t)-Uj1(23,t)],Lij1(24,t)+Uj1(23,t));
        cone([2*Pij1(25,t);  2*Qij1(25,t);  Lij1(25,t)-Uj1(5,t)], Lij1(25,t)+Uj1(5,t)) ;
        cone([2*Pij1(26,t);  2*Qij1(26,t);  Lij1(26,t)-Uj1(25,t)],Lij1(26,t)+Uj1(25,t)) ;
        cone([2*Pij1(27,t);  2*Qij1(27,t);  Lij1(27,t)-Uj1(26,t)],Lij1(27,t)+Uj1(26,t)) ;
        cone([2*Pij1(28,t);  2*Qij1(28,t);  Lij1(28,t)-Uj1(27,t)],Lij1(28,t)+Uj1(27,t)) ;
        cone([2*Pij1(29,t);  2*Qij1(29,t);  Lij1(29,t)-Uj1(28,t)],Lij1(29,t)+Uj1(28,t)) ;
        cone([2*Pij1(30,t);  2*Qij1(30,t);  Lij1(30,t)-Uj1(29,t)],Lij1(30,t)+Uj1(29,t)) ;
        cone([2*Pij1(31,t);  2*Qij1(31,t);  Lij1(31,t)-Uj1(30,t)],Lij1(31,t)+Uj1(30,t)) ;
        cone([2*Pij1(32,t);  2*Qij1(32,t);  Lij1(32,t)-Uj1(31,t)],Lij1(32,t)+Uj1(31,t));

        
        
        ];
end


% 主变压器运行约束
%   ADN1网变压器
cc=1;
for t=1:T
 F=[ F;
    sum(bk_1(:,cc))==1; %变压器只能工作在唯一档位
     % 档位变化量约束 
    Kt_1(1,1)==0;  Kt_1(1,cc+1)==Kt1(1,cc);
    Ktp1(:,cc)-Ktn1(:,cc)==Kt1(:,cc)-Kt_1(:,cc);  Ktp1(:,cc)>=0; Ktn1(:,cc)>=0;
    Kt1(1,cc)==0*bk_1(1,cc)+1*bk_1(2,cc)+2*bk_1(3,cc)+3*bk_1(4,cc)+4*bk_1(5,cc)+5*bk_1(6,cc)+6*bk_1(7,cc)+7*bk_1(8,cc)+8*bk_1(9,cc)-Klimit1; %档位与实际档位的映射
%变压器输出电压满足U0=Uin*(k0+Kt*deltk)
   1.0609==(k0+(0-Klimit1)*deltk)^2*Uc_1(1,cc)+(k0+(1-Klimit1)*deltk)^2*Uc_1(2,cc)+(k0+(2-Klimit1)*deltk)^2*Uc_1(3,cc)+(k0+(3-Klimit1)*deltk)^2*Uc_1(4,cc)+(k0+(4-Klimit1)*deltk)^2*Uc_1(5,cc)+(k0+(5-Klimit1)*deltk)^2*Uc_1(6,cc)+(k0+(6-Klimit1)*deltk)^2*Uc_1(7,cc)+(k0+(7-Klimit1)*deltk)^2*Uc_1(8,cc)+(k0+(8-Klimit1)*deltk)^2*Uc_1(9,cc);
    UTmin*bk_1(:,cc)<=Uc_1(:,cc);  Uc_1(:,cc)<=UTmax*bk_1(:,cc); %对于选中档位的约束
    UTmin*(1-bk_1(:,cc))<=UT1(1,cc)-Uc_1(:,cc);  UT1(1,cc)-Uc_1(:,cc)<=UTmax*(1-bk_1(:,cc));%对于未选中档位的约束

];
end

%% 调用求解器计算潮流
        

       options=sdpsettings('solver','gurobi') ;      
       optimize(F,  f1,  options);  
       TF=strcmp(ans.info,'Successfully solved (GUROBI)');

   for t = 1:T
        z_soc1_new(:,t+1) = Ess1_local(:,t+1)+lam1(:,t+1)/rho_t;
        % z_soc_new = Ess_local+lam2/rho_t;
        z_soc1_new = min(max(z_soc1_new(:,t+1),0.2),0.9) ;  
        z_soc1_new(:,1) = 0.6;   %初始soc
        z_soc1_new(:,end) = min(max(z_soc1_new(:,end),0.4),0.6) ;%终端约束
  %乘子更新
        lam1 = lam1+ rho_t*(Ess1_local-z_soc1_new);
 end
%检查收敛性

    prim_res = norm( Ess1_local - z_soc1_new , 'fro');
    dual_res = norm( z_soc1_new -z_soc1, 'fro');
    fprintf('ADMM iter %d: primal_res = %.3e,dual_res = %.3e\n',prim_res,dual_res );
    
    if (prim_res < tol)&&(dual_res < tol)
        disp('ADMM 收敛');
      
    end
end
      Unode1=sqrt(max(0,value(Uj1)));
    f1_val=value(f1);
    local_P1=value(PIEEE_1);local_Q1=value(QIEEE_1);local_U1=value(UjIEEE_1);  
    Pch1= value(P_ch1);Ess0 =value(Ess1);
    Pdis1= value(P_dis1);

    z_soc1_hist=z_soc1_new;
    z_soc1 = z_soc1_new;%次轮主变量
end
