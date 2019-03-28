function [vx,vy,vz,dvx,dvy,dvz] = HaloOdeGPU(f,x,y,z,vx,vy,vz,mu,e)
%   根据 HaloOde.m 修改使其可以在 GPU 上积分
%
%   created by PH at 2016-03-06:2047
%

%%
%Xdot = HaloOde(f,X,mu,e)
%   ERTBP 动力学方程
%   Xdot = HaloOde(f,X,mu)
%   Halo 轨道 ODE ，默认e=0
%
%   注意：勿轻易修改，小的改动会导致小的数值误差，先前计算的 ME-Halo 轨道库的精度就会受到影响
%
%   ！！！ 2014-08-12:1246 后对 RTBPOde.m 进行修改
%
% last modified by PH at 2013-10-12:1046
% last modified by PH at 2013-11-18:精简了最后的计算，引入了中间变量r1,r2,kappa
% last modified by PH at 2013-12-06:修改了注释
% last modified by PH at 2014-03-05 再次检验了模型。


%% 由symbolic.m推导来的ODE

r2cubic = ((mu + x - 1)^2 + y^2 + z^2)^(3/2);
r1cubic = ((mu + x)^2 + y^2 + z^2)^(3/2);
kappa = (e*cos(f) + 1);

dvx =  2*vy + ( x + ((2*mu + 2*x)*(mu - 1))/(2*r1cubic) - (mu*(2*mu + 2*x - 2))/(2*r2cubic) ) / kappa;
dvy = -2*vx + ( y - (mu*y)/r2cubic + (y*(mu - 1))/r1cubic ) / kappa;
dvz = -   z + ( z - (mu*z)/r2cubic + (z*(mu - 1))/r1cubic ) / kappa;

end