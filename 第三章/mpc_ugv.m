clc;
clear;

%%
%参考轨迹生成
N = 100;       %参考轨迹点数量
T = 0.05;      %采样周期
Xout = zeros(N,3);
Tout = zeros(N,1);
for k = 1:1:N
    Xout(k,1) = k*T;        %参考x坐标
    Xout(k,2) = 2;          %参考y坐标
    Xout(k,3) = 0;          %参考航向角
    Tout(k,1) = (k-1)*T;    %第k个参考点时刻
end
%%
%控制系统基本情况介绍
Nx = 3;        %状态量个数
Nu = 2;        %控制量个数
[Nr,Nc] = size(Xout);    %Nr=100,Nc=3
Tsim = 20;     %仿真时间
X0 = [0 0 pi/3]; %车辆初始状态
L = 1;         %车辆轴距
vd1 = 1;       %参考系统的纵向速度
vd2 = 0;       %参考系统前轮偏角
%%
%根据控制系统的维度信息，提前定义好相关矩阵并赋值
x_real = zeros(Nr,Nc);   %状态量实测信息
x_piao = zeros(Nr,Nc);   %状态量预测信息
u_real = zeros(Nr,Nu);   %控制量实际信息
u_piao = zeros(Nr,Nu);   %控制量预测信息
x_real(1,:) = X0;           
x_piao(1,:) = x_real(1,:)-Xout(1,:);
X_PIAO = zeros(Nr,Nx*Tsim);
XXX = zeros(Nr,Nx*Tsim); %用于保存每个时刻预测的所有状态值
q = [1 0 0;0 1 0;0 0 0.5];
Q_cell = cell(Tsim,Tsim);
for i = 1:1:Tsim
    for j = 1:1:Tsim
        if i == j
            Q_cell{i,j} = q;
        else
            Q_cell{i,j} = zeros(Nx,Nx);
        end
    end
end
Q = cell2mat(Q_cell);           %目标函数中Q矩阵
R = 0.1 * eye(Nu*Tsim,Nu*Tsim); %目标函数中R矩阵
%%
%模型预测控制的主体，在每一周期开始，获取系统当前状态量，更新状态空间方程。
%求解一个标准二次规划问题。将所得第一个解施加到系统上。
for i = 1:1:Nr
    t_d = Xout(i,3);
    a = [1 0 -vd1*sin(t_d)*T;
         0 1 vd1*cos(t_d)*T;
         0 0 1;];               %状态空间中的矩阵A
    b = [cos(t_d)*T 0;
         sin(t_d)*T 0;
         vd2*T/L vd1*T/(cos(vd2)^2);]; %状态空间中的矩阵B
    A_cell = cell(Tsim,1);
    B_cell = cell(Tsim,Tsim);
    for j = 1:1:Tsim
        A_cell{j,1} = a^j;
        for k = 1:1:Tsim
            if k <= j
                B_cell{j,k}=(a^(j-k))*b;
            else
                B_cell{j,k}=zeros(Nx,Nu);
            end
        end
    end
    A = cell2mat(A_cell);   %预测模型中A矩阵
    B = cell2mat(B_cell);   %预测模型中B矩阵
    
    H = 2 * (B' * Q * B + R);
    f = 2 * B' * Q * A * x_piao(i,:)';
    A_cons = [];
    b_cons = [];
    lb = [-2.3; -0.64];       %控制量下界
    ub = [0.3 ; 0.64];        %控制量上界
    [X,fval(i,1), exitflag(i,1),output(i,1)] = quadprog(H,f,A_cons,b_cons,[],[],lb,ub);
    X_PIAO(i,:) = (A*x_piao(i,:)'+ B * X)';
    if i+j < Nr
        for j = 1:1:Tsim
            XXX(i,1 + 3 * (j-1))=X_PIAO(i,1+3*(j-1))+Xout(i+j,1);
            XXX(i,2 + 3 * (j-1))=X_PIAO(i,2+3*(j-1))+Xout(i+j,2);
            XXX(i,3 + 3 * (j-1))=X_PIAO(i,3+3*(j-1))+Xout(i+j,3);
        end
    else
        for j = 1:1:Tsim
            XXX(i,1 + 3 * (j-1))=X_PIAO(i,1+3*(j-1))+Xout(Nr,1);
            XXX(i,2 + 3 * (j-1))=X_PIAO(i,2+3*(j-1))+Xout(Nr,2);
            XXX(i,3 + 3 * (j-1))=X_PIAO(i,3+3*(j-1))+Xout(Nr,3);
        end
    end
    u_piao(i,1) = X(1,1);   %将控制序列的第一个元素作为预测控制量
    u_piao(i,2) = X(2,1);
    Tvec = [0 :0.05:4];
    X00 = x_real(i,:);
    vd11 = vd1 + u_piao(i,1);    %
    vd22 = vd2 + u_piao(i,2);
    XOUT = dsolve('Dx - vd11 * cos(z)=0','Dy - vd11 * sin(z) = 0','Dz - vd22 = 0', 'x(0) = X00(1)','y(0) = X00(2)','z(0) = X00(3)');
    t = T;
    x_real(i+1,1) = eval(XOUT.x);
    x_real(i+1,2) = eval(XOUT.y);
    x_real(i+1,3) = eval(XOUT.z);
    if(i<Nr)
        x_piao(i+1,:) = x_real(i+1,:)-Xout(i+1,:);
    end
    u_real(i,1) = vd1 + u_piao(i,1);
    u_real(i,2) = vd2 + u_piao(i,2);

    figure(1)
    plot(Xout(1:Nr,1),Xout(1:Nr,2));       %绘制参考轨迹
    hold on;
    plot(x_real(i,1),x_real(i,2),'r*');    %绘制真实轨迹
    xlabel('横向位置X[m]');
%     axis([-15 -13]);
    ylabel('纵向位置Y[m]');
    hold on;
    for k= 1:1:Tsim
        X(i,k+1) = XXX(i,1+3*(k-1));
        Y(i,k+1) = XXX(i,2+3*(k-1));
    end
    X(i,1) = x_real(i,1);
    Y(i,1) = x_real(i,2);
    plot(X(i,:),Y(i,:),'yo');
    hold on
end
 %%
 %系统状态量随时间变化
 figure(2)
 subplot(3,1,1)
 plot(Tout(1:Nr),Xout(1:Nr,1));
 hold on
 plot(Tout(1:Nr),x_real(1:Nr,1),'r');
 xlabel('时间T');
 ylabel('横向位置X');
 subplot(3,1,2);
 plot(Tout(1:Nr),Xout(1:Nr,2));
 hold on
 plot(Tout(1:Nr),x_real(1:Nr,2),'r');
 xlabel('时间T');
 ylabel('纵向位置X');
  subplot(3,1,3);
 plot(Tout(1:Nr),Xout(1:Nr,3));
 hold on
 plot(Tout(1:Nr),x_real(1:Nr,3),'r');
 xlabel('时间T');
 ylabel('\theta');

 %控制量随时间变化
 figure(3)
 subplot(2,1,1);
 plot(Tout(1:Nr),u_real(1:Nr,1),'r');
 xlabel('时间T');
 ylabel('纵向速度');
 subplot(2,1,2);
 plot(Tout(1:Nr),u_real(1:Nr,2),'r');
 xlabel('时间T');
 ylabel('角速度');
 
 %状态量偏差随时间变化
 figure(4)
 subplot(3,1,1);
 plot(Tout(1:Nr),x_piao(1:Nr,1));
 xlabel('时间T');
 ylabel('e(x)');
 subplot(3,1,2);
 plot(Tout(1:Nr),x_piao(1:Nr,2));
 xlabel('时间T');
 ylabel('e(y)');
  subplot(3,1,3);
 plot(Tout(1:Nr),x_piao(1:Nr,3));
 xlabel('时间T');
 ylabel('e(\theta)');