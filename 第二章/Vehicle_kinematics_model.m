function [vx,vy,wv] = Vehicle_kinematics_model(v,theta,l)
% v:�ٶ�
% theta:ǰ��ƫ��
% l:���
% vx:x�����ٶ�
% vy:y�����ٶ�
% wv:��ڽ��ٶ�
degree2rad=pi/180;
w_rad = v*tan(theta*degree2rad/l);
vx = v*cos(w_rad);
vy = v*sin(w_rad);
wv = w_rad;
end

