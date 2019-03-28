%% Circular/Elliptic Restricted Three-Body Problem
%   six variables
%
%   simple means no extra parameters are supported yet
%   unsupported features:
%       event
%       normcontrol
%       NonNegative
%
%   most basic codes are borrowed from ode45.m
%
%
%   created by PH at 2016-03-06:1515
%   updated by PH at 2019/03/27 22:11:00

function [t,y1,y2,y3,y4,y5,y6,flag_success] = ode45_GPU_RTBP(t0,tfinal,y10,y20,y30,y40,y50,y60,AbsTol,RelTol,MaxStep, mu, e)
% test version

% setting parameters
%   codes borrowed from odearguments.m
tdir = ( tfinal-t0 ) / abs( tfinal-t0 );
atol = AbsTol;
rtol = RelTol;
threshold = atol / rtol;
htspan = abs(tfinal - t0);
hmax = min(abs(tfinal-t0), MaxStep);
hmin = 16*eps(t0);
absh = min(hmax, htspan);
absh = max(absh, hmin);

% Initialize method parameters.
pow = 1/5;
% A = [1/5, 3/10, 4/5, 8/9, 1, 1];
A1 = 1/5; A2 = 3/10; A3 = 4/5; A4 = 8/9; A5 = 1; A6 = 1;
% B = [
%     1/5         3/40    44/45   19372/6561      9017/3168       35/384
%     0           9/40    -56/15  -25360/2187     -355/33         0
%     0           0       32/9    64448/6561      46732/5247      500/1113
%     0           0       0       -212/729        49/176          125/192
%     0           0       0       0               -5103/18656     -2187/6784
%     0           0       0       0               0               11/84
%     0           0       0       0               0               0
%     ]; % original ode45 code
% simply expand:
% B11 = 1/5; B12 = 3/40; B13 = 44/45; B14 = 19372/6561; B15 = 9017/3168; B16 = 35/384;
% B21 = 0; B22 = 9/40; B23 = -56/15; B24 = -25360/2187; B25 = -355/33; B26 = 0;
% B31 = 0; B32 = 0; B33 = 32/9; B34 = 64448/6561; B35 = 46732/5247; B36 = 500/1113;
% B41 = 0; B42 = 0; B43 = 0; B44 = -212/729; B45 = 49/176; B46 = 125/192;
% B51 = 0; B52 = 0; B53 = 0; B54 = 0; B55 = -5103/18656; B56 = -2187/6784;
% B61 = 0; B62 = 0; B63 = 0; B64 = 0; B65 = 0; B66 = 11/84;
% B71 = 0; B72 = 0; B73 = 0; B74 = 0; B75 = 0; B76 = 0;
% simplify: remove zeros and useless ones
B11 = 1/5; B12 = 3/40; B13 = 44/45; B14 = 19372/6561; B15 = 9017/3168; B16 = 35/384;
    B22 = 9/40; B23 = -56/15; B24 = -25360/2187; B25 = -355/33; B26 = 0;
        B33 = 32/9; B34 = 64448/6561; B35 = 46732/5247; B36 = 500/1113;
            B44 = -212/729; B45 = 49/176; B46 = 125/192;
                B55 = -5103/18656; B56 = -2187/6784;
                    B66 = 11/84;
% E = [71/57600; 0; -71/16695; 71/1920; -17253/339200; 22/525; -1/40]; % original ode45 code
E1 = 71/57600; E2 = 0; E3 = -71/16695; E4 = 71/1920; E5 = -17253/339200; E6 = 22/525; E7 = -1/40;
% f = zeros(neq,7); % original ode45 code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 修改此处来实现多变量输入
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RTBP
[f10,f20,f30,f40,f50,f60] = HaloOdeGPU(t0,y10,y20,y30,y40,y50,y60,mu,e); % specilized ODE function, only taking element-wise input
f11 = f10; f12 = 0; f13 = 0; f14 = 0; f15 = 0; f16 = 0; f17 = 0; % 行数为维数
f21 = f20; f22 = 0; f23 = 0; f24 = 0; f25 = 0; f26 = 0; f27 = 0; % 行数为维数
f31 = f30; f32 = 0; f33 = 0; f34 = 0; f35 = 0; f36 = 0; f37 = 0; % 行数为维数
f41 = f40; f42 = 0; f43 = 0; f44 = 0; f45 = 0; f46 = 0; f47 = 0; % 行数为维数
f51 = f50; f52 = 0; f53 = 0; f54 = 0; f55 = 0; f56 = 0; f57 = 0; % 行数为维数
f61 = f60; f62 = 0; f63 = 0; f64 = 0; f65 = 0; f66 = 0; f67 = 0; % 行数为维数

% THE MAIN LOOP

t = t0;
% sinble variable test
% y1 = y10;
% RTBP
y1 = y10; y2 = y20; y3 = y30; y4 = y40; y5 = y50; y6 = y60;

nfailed = 0; % 自动控制步长失败的次数
nsteps = 0; % 积分步长计数

err = 0;
tnew = 0;
% single variable test
% y1new = 0;
% RTBP
y1new = 0; y2new = 0; y3new = 0; y4new = 0; y5new = 0; y6new = 0;

done = false;
while ~done
    
    % By default, hmin is a small number such that t+hmin is only slightly
    % different than t.  It might be 0 if t is 0.
    hmin = 16*eps(t);
    absh = min(hmax, max(hmin, absh));    % couldn't limit absh until new hmin
    h = tdir * absh;
    
    % Stretch the step if within 10% of tfinal-t.
    % 如果距离最终点还差不到 1.1*h，则扩大步长直接到达终点
    if 1.1*absh >= abs(tfinal - t)
        h = tfinal - t;
        absh = abs(h);
        done = true;
    end
    
    % LOOP FOR ADVANCING ONE STEP.
    % 循环直到向前积分一步成功，有可能需要减小步长，所以嵌套在一个 while 循环中，当成功时用 break 跳出
    nofailed = true;                      % no failed attempts
    while true
%         hA = h * A;  % original ode45 code
%         hB = h * B;
%         f(:,2) = feval(ode,t+h*A1,y+f*hB(:,1));
%         f(:,3) = feval(ode,t+h*A2,y+f*hB(:,2));
%         f(:,4) = feval(ode,t+h*A3,y+f*hB(:,3));
%         f(:,5) = feval(ode,t+h*A4,y+f*hB(:,4));
%         f(:,6) = feval(ode,t+h*A5,y+f*hB(:,5));  % original ode45 code
        % 写成分量形式
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% 修改此处来实现多变量输入
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % RTBP
        [f12,f22,f32,f42,f52,f62] = HaloOdeGPU(t+h*A1,...
            y1+h*(f11*B11),...
            y2+h*(f21*B11),...
            y3+h*(f31*B11),...
            y4+h*(f41*B11),...
            y5+h*(f51*B11),...
            y6+h*(f61*B11), mu, e);
        [f13,f23,f33,f43,f53,f63] = HaloOdeGPU(t+h*A2,...
            y1+h*(f11*B12+f12*B22),...
            y2+h*(f21*B12+f22*B22),...
            y3+h*(f31*B12+f32*B22),...
            y4+h*(f41*B12+f42*B22),...
            y5+h*(f51*B12+f52*B22),...
            y6+h*(f61*B12+f62*B22), mu, e);
        [f14,f24,f34,f44,f54,f64] = HaloOdeGPU(t+h*A3,...
            y1+h*(f11*B13+f12*B23+f13*B33),...
            y2+h*(f21*B13+f22*B23+f23*B33),...
            y3+h*(f31*B13+f32*B23+f33*B33),...
            y4+h*(f41*B13+f42*B23+f43*B33),...
            y5+h*(f51*B13+f52*B23+f53*B33),...
            y6+h*(f61*B13+f62*B23+f63*B33), mu, e);        
        [f15,f25,f35,f45,f55,f65] = HaloOdeGPU(t+h*A4,...
            y1+h*(f11*B14+f12*B24+f13*B34+f14*B44),...
            y2+h*(f21*B14+f22*B24+f23*B34+f24*B44),...
            y3+h*(f31*B14+f32*B24+f33*B34+f34*B44),...
            y4+h*(f41*B14+f42*B24+f43*B34+f44*B44),...
            y5+h*(f51*B14+f52*B24+f53*B34+f54*B44),...
            y6+h*(f61*B14+f62*B24+f63*B34+f64*B44), mu, e);        
        [f16,f26,f36,f46,f56,f66] = HaloOdeGPU(t+h*A5,...
            y1+h*(f11*B15+f12*B25+f13*B35+f14*B45+f15*B55),...
            y2+h*(f21*B15+f22*B25+f23*B35+f24*B45+f25*B55),...
            y3+h*(f31*B15+f32*B25+f33*B35+f34*B45+f35*B55),...
            y4+h*(f41*B15+f42*B25+f43*B35+f44*B45+f45*B55),...
            y5+h*(f51*B15+f52*B25+f53*B35+f54*B45+f55*B55),...
            y6+h*(f61*B15+f62*B25+f63*B35+f64*B45+f65*B55), mu, e);        
        
        % 计算下一步的时间
        tnew = t + h*A6;
        if done
            tnew = tfinal;   % Hit end point exactly.
        end
        h = tnew - t;      % Purify h.
        
        % 计算下一步的积分状态
%         y1new = y + f*hB(:,6);  % original ode45 code
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% 修改此处来实现多变量输入
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % RTBP
        y1new = y1 + h*(f11*B16+f12*B26+f13*B36+f14*B46+f15*B56+f16*B66);
        y2new = y2 + h*(f21*B16+f22*B26+f23*B36+f24*B46+f25*B56+f26*B66);
        y3new = y3 + h*(f31*B16+f32*B26+f33*B36+f34*B46+f35*B56+f36*B66);
        y4new = y4 + h*(f41*B16+f42*B26+f43*B36+f44*B46+f45*B56+f46*B66);
        y5new = y5 + h*(f51*B16+f52*B26+f53*B36+f54*B46+f55*B56+f56*B66);
        y6new = y6 + h*(f61*B16+f62*B26+f63*B36+f64*B46+f65*B56+f66*B66);
        [f17,f27,f37,f47,f57,f67] = HaloOdeGPU(tnew,y1new,y2new,y3new,y4new,y5new,y6new,mu,e);
        
        % Estimate the error.
%         err = absh * norm((f * E) ./ max(max(abs(y),abs(y1new)),threshold),inf); % original ode45 code
        % RTBP
        err1 = abs( (f11*E1+f12*E2+f13*E3+f14*E4+f15*E5+f16*E6+f17*E7) ./ max(max(abs(y1),abs(y1new)),threshold) );
        err2 = abs( (f21*E1+f22*E2+f23*E3+f24*E4+f25*E5+f26*E6+f27*E7) ./ max(max(abs(y2),abs(y2new)),threshold) );
        err3 = abs( (f31*E1+f32*E2+f33*E3+f34*E4+f35*E5+f36*E6+f37*E7) ./ max(max(abs(y3),abs(y3new)),threshold) );
        err4 = abs( (f41*E1+f42*E2+f43*E3+f44*E4+f45*E5+f46*E6+f47*E7) ./ max(max(abs(y4),abs(y4new)),threshold) );
        err5 = abs( (f51*E1+f52*E2+f53*E3+f54*E4+f55*E5+f56*E6+f57*E7) ./ max(max(abs(y5),abs(y5new)),threshold) );
        err6 = abs( (f61*E1+f62*E2+f63*E3+f64*E4+f65*E5+f66*E6+f67*E7) ./ max(max(abs(y6),abs(y6new)),threshold) );
        err = absh * max(err6,max(err5,max(err4,max(err3,max(err2,max(err1,0))))));
        
        % Accept the solution only if the weighted error is no more than the
        % tolerance rtol.  Estimate an h that will yield an error of rtol on
        % the next step or the next try at taking this step, as the case may be,
        % and use 0.8 of this value to avoid failures.
        if err > rtol                       % Failed step
            nfailed = nfailed + 1;
            if absh <= hmin
                flag_success = 0;
                return;
            end
            if nofailed
                % 第一次减小步长
                nofailed = false;
                absh = max(hmin, absh * max(0.1, 0.8*(rtol/err)^pow));
            else
                % 再次减小步长
                absh = max(hmin, 0.5 * absh);
            end
            h = tdir * absh;
            done = false;
        else                                % Successful step
            break;
        end
    end
    nsteps = nsteps + 1;
    
    % If there were no failures compute a new h.
    if nofailed
        % Note that absh may shrink by 0.8, and that err may be 0.
        temp = 1.25*(err/rtol)^pow;
        if temp > 0.2
            absh = absh / temp;
        else
            absh = 5.0*absh;
        end
    end
    
    % Advance the integration one step.
    t = tnew;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% 修改此处来实现多变量输入
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CRTBP
    y1 = y1new;
    y2 = y2new;
    y3 = y3new;
    y4 = y4new;
    y5 = y5new;
    y6 = y6new;
    f11 = f17;  % Already have f(tnew,y1new)
    f21 = f27;  % Already have f(tnew,y1new)
    f31 = f37;  % Already have f(tnew,y1new)
    f41 = f47;  % Already have f(tnew,y1new)
    f51 = f57;  % Already have f(tnew,y1new)
    f61 = f67;  % Already have f(tnew,y1new)
    
end
    
    % generate outputs
    % t is the output
    % y1 is the output
    flag_success = 1;
    
end
