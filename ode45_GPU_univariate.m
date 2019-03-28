% % single variable test
%
%   simple means no extra parameters are supported yet
%   unsupported features:
%       event
%       normcontrol
%       NonNegative
%
%   most basic codes are borrowed from ode45.m
%
%   下一步工作：
%       整理代码，结合进多变量进行验证
%       
%
%
%   created by PH at 2016-03-06:1515
%

function [t,y1,flag_success] = ode45_GPU_univariate(t0,tfinal,y10,AbsTol,RelTol,MaxStep)
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
%     ];
% B11 = 1/5; B12 = 3/40; B13 = 44/45; B14 = 19372/6561; B15 = 9017/3168; B16 = 35/384;
% B21 = 0; B22 = 9/40; B23 = -56/15; B24 = -25360/2187; B25 = -355/33; B26 = 0;
% B31 = 0; B32 = 0; B33 = 32/9; B34 = 64448/6561; B35 = 46732/5247; B36 = 500/1113;
% B41 = 0; B42 = 0; B43 = 0; B44 = -212/729; B45 = 49/176; B46 = 125/192;
% B51 = 0; B52 = 0; B53 = 0; B54 = 0; B55 = -5103/18656; B56 = -2187/6784;
% B61 = 0; B62 = 0; B63 = 0; B64 = 0; B65 = 0; B66 = 11/84;
% B71 = 0; B72 = 0; B73 = 0; B74 = 0; B75 = 0; B76 = 0;
% simplified: remove zeros and useless
B11 = 1/5; B12 = 3/40; B13 = 44/45; B14 = 19372/6561; B15 = 9017/3168; B16 = 35/384;
    B22 = 9/40; B23 = -56/15; B24 = -25360/2187; B25 = -355/33; B26 = 0;
        B33 = 32/9; B34 = 64448/6561; B35 = 46732/5247; B36 = 500/1113;
            B44 = -212/729; B45 = 49/176; B46 = 125/192;
                B55 = -5103/18656; B56 = -2187/6784;
                    B66 = 11/84;
% E = [71/57600; 0; -71/16695; 71/1920; -17253/339200; 22/525; -1/40];
E1 = 71/57600; E2 = 0; E3 = -71/16695; E4 = 71/1920; E5 = -17253/339200; E6 = 22/525; E7 = -1/40;
% f = zeros(neq,7);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 修改此处来实现多变量输入
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% single variable test
f01 = UnivariateODE(t0,y10); % <--------------------------------------------------- 改动力学方程
f11 = f01; f12 = 0; f13 = 0; f14 = 0; f15 = 0; f16 = 0; f17 = 0; % 行数为维数

% THE MAIN LOOP

t = t0;
% sinble variable test
y1 = y10;

nfailed = 0; % 自动控制步长失败的次数
nsteps = 0; % 积分步长计数

err = 0;
tnew = 0;
% single variable test
y1new = 0;

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
%         hA = h * A;
%         hB = h * B;
%         f(:,2) = feval(ode,t+h*A1,y+f*hB(:,1));
%         f(:,3) = feval(ode,t+h*A2,y+f*hB(:,2));
%         f(:,4) = feval(ode,t+h*A3,y+f*hB(:,3));
%         f(:,5) = feval(ode,t+h*A4,y+f*hB(:,4));
%         f(:,6) = feval(ode,t+h*A5,y+f*hB(:,5));
        % 写成分量形式
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% 修改此处来实现多变量输入
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % single variable test
%         f12 = UnivariateODE(t+h*A1,y1+h*(f11*B11+f12*B21+f13*B31+f14*B41+f15*B51+f16*B61+f17*B71));
%         f13 = UnivariateODE(t+h*A2,y1+h*(f11*B12+f12*B22+f13*B32+f14*B42+f15*B52+f16*B62+f17*B72));
%         f14 = UnivariateODE(t+h*A3,y1+h*(f11*B13+f12*B23+f13*B33+f14*B43+f15*B53+f16*B63+f17*B73));
%         f15 = UnivariateODE(t+h*A4,y1+h*(f11*B14+f12*B24+f13*B34+f14*B44+f15*B54+f16*B64+f17*B74));
%         f16 = UnivariateODE(t+h*A5,y1+h*(f11*B15+f12*B25+f13*B35+f14*B45+f15*B55+f16*B65+f17*B75));
%         % simplified
        f12 = UnivariateODE(t+h*A1,y1+h*(f11*B11)); % <------------------------------------------------------------------------------- 改动力学方程
        f13 = UnivariateODE(t+h*A2,y1+h*(f11*B12+f12*B22)); % <------------------------------------------------------------------------ 改动力学方程
        f14 = UnivariateODE(t+h*A3,y1+h*(f11*B13+f12*B23+f13*B33)); % <----------------------------------------------------------------- 改动力学方程
        f15 = UnivariateODE(t+h*A4,y1+h*(f11*B14+f12*B24+f13*B34+f14*B44)); % <---------------------------------------------------------- 改动力学方程
        f16 = UnivariateODE(t+h*A5,y1+h*(f11*B15+f12*B25+f13*B35+f14*B45+f15*B55)); % <--------------------------------------------------- 改动力学方程

        % 计算下一步的时间
        tnew = t + h*A6;
        if done
            tnew = tfinal;   % Hit end point exactly.
        end
        h = tnew - t;      % Purify h.
        
        % 计算下一步的积分状态
%         y1new = y + f*hB(:,6);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% 修改此处来实现多变量输入
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % single variable test
        y1new = y1 + h*(f11*B16+f12*B26+f13*B36+f14*B46+f15*B56+f16*B66);
        f17 = UnivariateODE(tnew,y1new); % <--------------------------------------------------- 改动力学方程
        
        % Estimate the error.
%         err = absh * norm((f * E) ./ max(max(abs(y),abs(y1new)),threshold),inf);
        % single variable test
        err1 = abs( (f11*E1+f12*E2+f13*E3+f14*E4+f15*E5+f16*E6+f17*E7) ./ max(max(abs(y1),abs(y1new)),threshold) );
        err = absh * max(err1,0);
        
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
    % single variable test
    y1 = y1new;
    f11 = f17;  % Already have f(tnew,y1new)
    
end
    
    % generate outputs
    % t is output
    % y1 is output
    flag_success = 1;
    
end
