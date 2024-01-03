function [sys,x0,str,ts]=s_function(t,x,u,flag)
switch flag,
case 0,
[sys,x0,str,ts]=mdlInitializeSizes;
case 1,
sys=mdlDerivatives(t,x,u);
case 3,
sys=mdlOutputs(t,x,u);
case {2, 4, 9 }
sys = [];
otherwise
error(['Unhandled flag = ',num2str(flag)]);
end
end
function [sys,x0,str,ts]=mdlInitializeSizes

global  g m1 m2 a1 a2


sizes = simsizes;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
sizes.NumContStates = 4;
sizes.NumDiscStates = 0;
sizes.NumOutputs = 11;%
sizes.NumInputs =5;
sizes.DirFeedthrough = 0;
sizes.NumSampleTimes = 0;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
sys=simsizes(sizes);
x0=[0.1 0.03 0.03 0.03];
%x0=[0.3 0.5 0.5 0];
str=[];
ts=[];

g=9.8;
m1 = 2.6366;%5 ; 
m2 = 0.82,%2.6366+0.82 ;
a1 = 0.19,%0.38 ;
a2 = 0.0608,%0.0608+0.38 ;
global adk0 adk1 bdk1 adk2 bdk2 adk3 bdk3 adk4 bdk4 adk5 bdk5 wdk
global ada0 ada1 bda1 ada2 bda2 ada3 bda3 ada4 bda4 ada5 bda5 wda
global ag0 ag1 ag2 ag3 ag4 ag5 bg1 bg2 bg3 bg4 bg5 wg 
%ankle coeffcient (link 2)  
       ada0 =   -0.008241 ; %(-0.008938, -0.006626)
       ada1 =     -0.001481  ;%(-0.02222, -0.01857)
       bda1 =       0.1349  ;%(0.09676, 0.1)
       ada2 =      0.01367  ;%(0.008884, 0.01389)
       bda2 =    -0.1203  ;%(-0.1134, -0.1101)
       ada3 =    -0.05967  ;%(-0.06254, -0.05917)
       bda3 =    0.01831  ;%(0.008338, 0.01296)
       ada4 =     0.03788  ;%(0.03473, 0.03805)
       bda4 = -0.008675 ;%(0.0008274, 0.004929)
       ada5 =   0.0003328 ; %(0.004849, 0.008678)
       bda5 =    -0.02038; %(-0.02511, -0.02181)
       wda =          6.844   ;%(5.595, 5.605)
       
 %knee coeffcient (link 1)       
      adk0 =      0.4246  ;
       adk1 =     -0.1041 ;
       bdk1 =     -0.3332  ;
       adk2 =     -0.2508 ; 
       bdk2 =      0.1945  ;
       adk3 =     0.003604 ; 
       bdk3 =     0.08631 ; 
       adk4 =  -0.01391 ; 
       bdk4 =  0.01127;  
       adk5 =    0.00147 ; 
       bdk5 =     0.01428 ; 
       wdk =       6.845 ;
       
 %GrfCoeeficient
          ag0 =       12.86  ;%(12.77, 12.95)
       ag1 =      -1.719  ;%(-1.864, -1.573)
       bg1 =       1.707  ;%(1.569, 1.845)
       ag2 =       2.668  ;%(2.539, 2.797)
       bg2 =      0.1859  ;%(-0.02101, 0.3928)
       ag3 =     -0.2617  ;%(-0.3905, -0.1329)
       bg3 =      0.2617  ;%(0.1311, 0.3923)
       ag4 =      0.4435  ;%(0.2737, 0.6134)
       bg4 =      0.9948  ;%(0.85, 1.14)
       ag5 =     -0.2047  ;%(-0.355, -0.05449)
       bg5 =     -0.4625  ;%(-0.5912, -0.3338)
       wg =       6.841  ;%(6.819, 6.862)      

end
function sys=mdlDerivatives(t,x,u)
global g m1 m2 a1 a2 
global ag0 ag1 ag2 ag3 ag4 ag5 bg1 bg2 bg3 bg4 bg5 wg 
 global g m1 m2 a1 a2 

D = [(m1 + m2)*a1^2 + m2*a2^2 + 2*m2*a1*a2*cos(x(3)), m2*a2^2 + m2*a1*a2*cos(x(3));...
    m2*a2^2 + m2*a1*a2*cos(x(3)), m2*a2^2];
% Mbar = M - delM;
% C = [ -m2*a1*a2*(2*x(2)*x(4)+(x(4))^2)*sin(x(4));...
%     m2*a1*a2*(x(2))^2*sin(x(3))] ; 
C = [-x(4)*m2*a1*a2*sin(x(3)),-(x(2)+x(4))*m2*a1*a2*sin(x(3));(x(2))*m2*a1*a2*sin(x(3)),0]; 
    
G = [(m1+m2)*g*a1*cos(x(1))+ m2*g*a2*cos(x(1)+x(3));...
    m2*g*a2*cos(x(1)+x(3))];  
N = C+ G ; 


dq=[x(2);x(4)];
F=0.02*sign(dq);%
grf = ag0 + ag1*cos(t*wg) + bg1*sin(t*wg) + ...
               ag2*cos(2*t*wg) + bg2*sin(2*t*wg) + ag3*cos(3*t*wg) + bg3*sin(3*t*wg) + ...
               ag4*cos(4*t*wg) + bg4*sin(4*t*wg) + ag5*cos(5*t*wg) + bg5*sin(5*t*wg)...;

told=[0.1*sin(t);0.1*sin(t)];%

tor_grf= [0.1;1]*grf;
tol=[u(1);u(2)]; 

S=inv(D)*(tol+tor_grf-C*dq-G-F-told);

sys(1)=x(2);
sys(2)=S(1);
sys(3)=x(4);
sys(4)=S(2);
end
function sys=mdlOutputs(t,x,u)
global adk0 adk1 bdk1 adk2 bdk2 adk3 bdk3 adk4 bdk4 adk5 bdk5 wdk
global ada0 ada1 bda1 ada2 bda2 ada3 bda3 ada4 bda4 ada5 bda5 wda

global g m1 m2 a1 a2 Fai
D = [(m1 + m2)*a1^2 + m2*a2^2 + 2*m2*a1*a2*cos(x(3)), m2*a2^2 + m2*a1*a2*cos(x(3));...
    m2*a2^2 + m2*a1*a2*cos(x(3)), m2*a2^2];
% Mbar = M - delM;
C = [-x(4)*m2*a1*a2*sin(x(3)),-(x(2)+x(4))*m2*a1*a2*sin(x(3));(x(2))*m2*a1*a2*sin(x(3)),0]; 
G = [(m1+m2)*g*a1*cos(x(1))+ m2*g*a2*cos(x(1)+x(3));...
    m2*g*a2*cos(x(1)+x(3))];  
dq=[x(2);x(4)];
F=0.2*sign(dq);

qd1= adk0 + adk1*cos(wdk*t) + bdk1*sin(wdk*t) +  adk2*cos(2*wdk*t) +...
    bdk2*sin(2*wdk*t) + adk3*cos(3*wdk*t) + bdk3*sin(3*wdk*t)+adk4*cos(4*wdk*t) +...
    bdk4*sin(4*wdk*t) +adk5*cos(5*wdk*t) + bdk5*sin(5*wdk*t) ;

d_qd1= -adk1*wdk*sin(wdk*t) + bdk1*wdk*cos(wdk*t) -  adk2*2*wdk*sin(2*wdk*t) +...
    bdk2*2*wdk*cos(2*wdk*t) - adk3*3*wdk*sin(3*wdk*t) + bdk3*3*wdk*cos(3*wdk*t)- adk4*4*wdk*sin(4*wdk*t)+...
    4*wdk*bdk4*cos(4*wdk*t) -adk5*5*wdk*sin(5*wdk*t) + bdk5*5*wdk*cos(5*wdk*t) ;


dd_qd1= - adk1*wdk^2*cos(wdk*t) - bdk1*wdk^2*sin(wdk*t) -  adk2*wdk^2*4*cos(2*wdk*t) -...
    bdk2*wdk^2*4*sin(2*wdk*t) - adk3*wdk^2*9*cos(3*wdk*t) - bdk3*wdk^2*9*sin(3*wdk*t)-adk4*wdk^2*16*cos(4*wdk*t) -...
    bdk4*wdk^2*16*sin(4*wdk*t) -adk5*wdk^2*25*cos(5*wdk*t) - bdk5*wdk^2*25*sin(5*wdk*t) ;


%ankle--link2 
qd2= ada0 + ada1*cos(wda*t) + bda1*sin(wda*t) +  ada2*cos(2*wda*t) +...
    bda2*sin(2*wda*t) + ada3*cos(3*wda*t) + bda3*sin(3*wda*t)+ada4*cos(4*wda*t) +...
    bda4*sin(4*wda*t) +ada5*cos(5*wda*t) + bda5*sin(5*wda*t) ;

d_qd2=-ada1*wda*sin(wda*t) + bda1*wda*cos(wda*t) -  ada2*2*wda*sin(2*wda*t) +...
    bda2*2*wda*cos(2*wda*t) - ada3*3*wda*sin(3*wda*t) + bda3*3*wda*cos(3*wda*t)- ada4*4*wda*sin(4*wda*t)+...
    4*wda*bda4*cos(4*wda*t) -ada5*5*wda*sin(5*wda*t) + bda5*5*wda*cos(5*wda*t) ;

dd_qd2= -ada1*wda^2*cos(wda*t) - bda1*wda^2*sin(wda*t) -  ada2*wda^2*4*cos(2*wda*t) -...
    bda2*wda^2*4*sin(2*wda*t) - ada3*wda^2*9*cos(3*wda*t) - bda3*wda^2*9*sin(3*wda*t)-ada4*wda^2*16*cos(4*wda*t) -...
    bda4*wda^2*16*sin(4*wda*t) -ada5*wda^2*25*cos(5*wda*t) - bda5*wda^2*25*sin(5*wda*t) ;
q1=x(1);
d_q1=dq(1);
q2=x(3);
d_q2=dq(2);
q=[q1;q2];
e1=qd1-q1;
e2=qd2-q2;
de1=d_qd1-d_q1;
de2=d_qd2-d_q2;
e=[e1;e2];
de=[de1;de2];
Fai=20*eye(2); % lamda 
dqd=[d_qd1;d_qd2];
dqr=dqd+Fai*e;
ddqd=[dd_qd1;dd_qd2];
ddqr=ddqd+Fai*de;
% f=D*ddqr+C*dqr+G+F;
f=D*ddqr+C*dqr+G+F;
f_norm=norm(f);

sys(1)=x(1);
sys(2)=x(2);
sys(3)=x(3);
sys(4)=x(4);
sys(5)=f_norm;
sys(6)= f(1);
sys(7)= f(2);
sys(8)= e1
sys(9)= e2
sys(10)= de1
sys(11)= de2
end