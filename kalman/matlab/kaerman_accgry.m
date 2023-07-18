%% ���ùߵ�����ٶȣ�Ԥ��С����λ��
clc
clear all;
close  all;
%% ��������
%�ߵ����������
load('temp_acc.mat');   %�ߵ��ɵ���ʵ����

%�������������
% temp_acc = [[0:0.1:2],[1.9:-0.1:-2],[-1.9:0.1:0]];   %��������
% temp_acc = temp_acc';
% step=length(temp_acc);
% Voice=normrnd(0,0.1,1,step);
% Voice = Voice';
% temp_acc =temp_acc + Voice;

%% ���Ŀ�ʼ
step=length(temp_acc);
dtt = 1/250;   %�ߵ�����dtt=1/250    ��������dtt=0.1
A=[1,dtt,0;0,1,dtt;0,0,1]; %���Ըߴ���
H=[0,0,1];  %�۲���

position = zeros(1,step+1);   %��¼����
speed = zeros(1,step+1);
acc = zeros(1,step+1);

Q=[0.1,0,0;0,0.1,0;0,0,0.1];  %ϵͳ����
R=[0.1];   %��������
Pk=[1,0,0;0,1,0;0,0,1];  
x = [0 0 0]';

position(1) = x(1);
speed(1) = x(2);
acc(1) = x(3);

%% �������˲�
for i=1:step
    x_ = A * x ; %��
    Pk_ = A * Pk * A' + Q ; %��
    
    Kk = (Pk_ * H')/(H * Pk_ * H' + R);  %��
    x = x_ + Kk * (temp_acc(i)-H*x_);  %��
    Pk = (diag([1 1 1 ]) - Kk * H) * Pk_;  %��
    position(i+1) = x(1);
    speed(i+1) = x(2);
    acc(i+1) = x(3);
end

%% ��ͼ
figure(1)
plot([1:step+1]/250,position)    %����λ��
title('λ�ù���')

figure(2)
plot([1:step+1]/250,speed)  %�����ٶ�
title('�ٶȹ���')

figure(3)
plot([1:step+1]/250,[0;temp_acc]) %ʵ�ʼ��ٶ�
hold on;
plot([1:step+1]/250,acc) %���Ƽ��ٶ�
legend('ʵ�ʼ��ٶ�','���Ƽ��ٶ�')
title('���ٶȹ���')