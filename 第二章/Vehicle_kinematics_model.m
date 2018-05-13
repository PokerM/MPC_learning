function [vx,vy,wv] = Vehicle_kinematics_model(v,theta,l)
% v:速度
% theta:前轮偏角
% l:轴距
% vx:x方向速度
% vy:y方向速度
% wv:横摆角速度
degree2rad=pi/180;
w_rad = v*tan(theta*degree2rad/l);
vx = v*cos(w_rad);
vy = v*sin(w_rad);
wv = w_rad;
end

