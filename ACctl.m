function [sys,x0,str,ts] = spacemodel(t,x,u,flag)
switch flag,
case 0,
[sys,x0,str,ts]=mdlInitializeSizes;
case 1,
sys=mdlDerivatives(t,x,u);
case 3,
sys=mdlOutputs(t,x,u);
case {2,4,9}
sys=[];
otherwise
error(['Unhandled flag = ',num2str(flag)]);
end
end
function [sys,x0,str,ts]=mdlInitializeSizes
 global Kv Kp
 global  g m1t m2t a1 a2 m1 m2
 global node c1 c2 b Fai b1 b2 c
node=30;



c =[-0.4:0.027:0.4 ; -12:0.827:12 ;-1:0.068:1 ; -5:0.344:5 ; -70:4.8:70;...
    -0.3:0.02068:0.3 ; -15:1.0344:15 ;-0.3:0.02068:0.3 ; -3.5:0.24:3.5 ; -70:4.8:70] ;
b = [0.1:16.67:500] ;  %30 so far best working 
global ag0 ag1 ag2 ag3 ag4 ag5 bg1 bg2 bg3 bg4 bg5 wg 

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





Fai=20*eye(2);%lambda
g=9.8;
m1 = 2.6366;%5 ; 
m2 = 0.82;%2.6366+0.82 ;
a1 = 0.19;%0.38 a2 = 0.0608;%0.0608+0.38 ;

m1t = 2.85;%2.6366;%5 ; 
m2t = 1;%0.82;%2.6366+0.82 ;



Kv= 6*eye(2);

sizes = simsizes;
sizes.NumContStates = 2*node;
sizes.NumDiscStates = 0;
sizes.NumOutputs = 5;%3
sizes.NumInputs = 21;%
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1;
sys = simsizes(sizes);
x0 = 2*node;%
str = [];
ts = [0 0];
end
function sys=mdlDerivatives(t,x,u)
global node c b Fai c1 c2 b1 b2
global ag0 ag1 ag2 ag3 ag4 ag5 bg1 bg2 bg3 bg4 bg5 wg 
 
qd1=u(1);
d_qd1=u(2);
dd_qd1=u(3);
qd2=u(4);
d_qd2=u(5);
dd_qd2=u(6);
q1=u(7);
d_q1=u(8);
q2=u(9);
d_q2=u(10);
f1 = u(12); 
f2 = u(13); 
% f1_ctl = u(14); 
% f2_ctl = u(15) ; 
eC1=u(18);
eC2=u(19);
eA1=u(20);
eA2=u(21);
eC = [eC1;eC2] ;
eA = [eA1;eA2] ;
q=[q1;q2];
e1=qd1-q1;
e2=qd2-q2;
de1=d_qd1-d_q1;
de2=d_qd2-d_q2;
e=[e1;e2];
de=[de1;de2];
r=de+Fai*e;
qd=[qd1;qd2];
dqd=[d_qd1;d_qd2];
dqr=dqd+Fai*e;
ddqd=[dd_qd1;dd_qd2];
ddqr=ddqd+Fai*de;

z = [e(1);de(1);qd(1);dqd(1);ddqd(1);e(2);de(2);qd(2);dqd(2);ddqd(2)]
for j=1:1:node

 h(j)=exp(-norm(z-c(:,j))^2/(b(j)*b(j)));
end
F=15*eye(node);
for i=1:1:node
sys(i)=22*h(i)*r(1)-22*x(1)*norm(eA); % including eA eA
sys(i+node)=22*h(i)*r(2)-22*x(2)*norm(eA);


end
end

function sys=mdlOutputs(t,x,u)
% global node c b Fai
global node c b Fai Kv c1 c2 b1 b2 
global ag0 ag1 ag2 ag3 ag4 ag5 bg1 bg2 bg3 bg4 bg5 wg 
qd1=u(1);
d_qd1=u(2);
dd_qd1=u(3);
qd2=u(4);
d_qd2=u(5);
dd_qd2=u(6);
q1=u(7);
d_q1=u(8);
q2=u(9);
d_q2=u(10);
f1 = u(12); 
f2 = u(13); 
% f1_ctl = u(14); 
% f2_ctl = u(15) ; 
eC1=u(18);
eC2=u(19);
eA1=u(20);
eA2=u(21);
eC = [eC1;eC2] ;
eA = [eA1;eA2] ;


q=[q1;q2];
e1=qd1-q1;
e2=qd2-q2;
de1=d_qd1-d_q1;
de2=d_qd2-d_q2;
e=[e1;e2];
de=[de1;de2];
r=de+Fai*e;
qd=[qd1;qd2];
dqd=[d_qd1;d_qd2];
dqr=dqd+Fai*e;
ddqd=[dd_qd1;dd_qd2];
ddqr=ddqd+Fai*de;
W_f1=[x(1:node)]';
W_f2=[x(node+1:node*2)]';

z = [e(1);de(1);qd(1);dqd(1);ddqd(1);e(2);de(2);qd(2);dqd(2);ddqd(2)]
for j=1:1:node

 h(j)=exp(-norm(z-c(:,j))^2/(b(j)*b(j)));
end
% fn=[W_f1*h1';
% W_f2*h2'];
fn=[W_f1*h';
W_f2*h'];

if fn(1)>35 
    fn(1)=35 ;
end
if fn(2)>3.5 
    fn(2)=3.5 ;
end






epN=0.20;bd=0.1;
v=-(epN+bd)*sign(r);
grf = ag0 + ag1*cos(t*wg) + bg1*sin(t*wg) + ag2*cos(2*t*wg) + bg2*sin(2*t*wg) + ag3*cos(3*t*wg) + bg3*sin(3*t*wg) + ...
               ag4*cos(4*t*wg) + bg4*sin(4*t*wg) + ag5*cos(5*t*wg) + bg5*sin(5*t*wg);
tor_grf= grf.*[0.1;1];
tol=fn+Kv*r-v- tor_grf;

fn_norm=norm(fn);
sys(1)=tol(1);
sys(2)=tol(2);
sys(3)=fn_norm;
sys(4)=fn(1); 
sys(5) = fn(2); 
end