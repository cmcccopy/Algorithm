%% 利用惯导测加速度，预测小车的位置
clc
clear all;
close  all;
%% 导入数据
%惯导数据用这个
load('temp_acc.mat');   %惯导采的真实数据

%放真数据用这个
% temp_acc = [[0:0.1:2],[1.9:-0.1:-2],[-1.9:0.1:0]];   %仿真数据
% temp_acc = temp_acc';
% step=length(temp_acc);
% Voice=normrnd(0,0.1,1,step);
% Voice = Voice';
% temp_acc =temp_acc + Voice;

%% 正文开始
step=length(temp_acc);
dtt = 1/250;   %惯导数据dtt=1/250    仿真数据dtt=0.1
A=[1,dtt,0;0,1,dtt;0,0,1]; %忽略高次项
H=[0,0,1];  %观测阵

position = zeros(1,step+1);   %记录数据
speed = zeros(1,step+1);
acc = zeros(1,step+1);

Q=[0.1,0,0;0,0.1,0;0,0,0.1];  %系统噪声
R=[0.1];   %量测噪声
Pk=[1,0,0;0,1,0;0,0,1];  
x = [0 0 0]';

position(1) = x(1);
speed(1) = x(2);
acc(1) = x(3);

%% 卡尔曼滤波
for i=1:step
    x_ = A * x ; %①
    Pk_ = A * Pk * A' + Q ; %②
    
    Kk = (Pk_ * H')/(H * Pk_ * H' + R);  %③
    x = x_ + Kk * (temp_acc(i)-H*x_);  %④
    Pk = (diag([1 1 1 ]) - Kk * H) * Pk_;  %⑤
    position(i+1) = x(1);
    speed(i+1) = x(2);
    acc(i+1) = x(3);
end

%% 绘图
figure(1)
plot([1:step+1]/250,position)    %估计位置
title('位置估计')

figure(2)
plot([1:step+1]/250,speed)  %估计速度
title('速度估计')

figure(3)
plot([1:step+1]/250,[0;temp_acc]) %实际加速度
hold on;
plot([1:step+1]/250,acc) %估计加速度
legend('实际加速度','估计加速度')
title('加速度估计')