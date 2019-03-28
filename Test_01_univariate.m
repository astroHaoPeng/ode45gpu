%% 单变量改写ode45的测试
%
%   created by PH
%   updated by PH at 2019/03/27 22:39:00

%%
clear all
close all
fun = @(t,x)x;
tspan = [0,1]; nn = 100;
tspanvector = linspace(tspan(1),tspan(end),20);
y0 = 0.5;
AbsTol = 1e-6;
RelTol = 1e-6;

%% use ode45
odeOptions = odeset('AbsTol',AbsTol,'RelTol',RelTol);
tic;
[T,Y] = ode45(@UnivariateODE,tspanvector,y0,odeOptions);
toc;

%% use modifed ode45_GPU_univariate
tic;
y0GPU = gpuArray(y0);
[TGPU,YGPU] = arrayfun(@ode45_GPU_univariate,tspan(1),tspanvector,y0GPU,AbsTol,RelTol,0.1*abs(tspan(end)-tspan(1)));
% [TGPU,YGPU] = arrayfun(@ode45_GPU_univariate,tspan(1),tspan(end),y0GPU,AbsTol,RelTol,0.1*abs(tspan(end)-tspan(1)));
% [TGPU,YGPU] = feval(@ode45_GPU_univariate,tspan(1),tspan(2),y0GPU,AbsTol,RelTol,0.1*abs(tspan(end)-tspan(1)));
toc;

%% compare
figure(1);
set(gcf, 'Position', [50,100,900,500]);
clf;

subplot(121);
plot(T(:),Y(:),'b-o'); hold on;
plot(TGPU(:),YGPU(:),'r--*');
title('results: CPU vs. GPU');

subplot(122);
plot(T, Y(:)-YGPU(:), '-x')
title('difference: CPU - GPU');

