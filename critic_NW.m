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

% end sfuntmpl

%
%=============================================================================
% mdlInitializeSizes
% Return the sizes, initial conditions, and sample times for the S-function.
%=============================================================================
%
function [sys,x0,str,ts]=mdlInitializeSizes
global c_critic b_critic nodecritic
c_critic =[-0.4:0.027:0.4 ; -12:0.827:12 ;-1:0.068:1 ; -5:0.344:5 ; -70:4.8:70;...
    -0.3:0.02068:0.3 ; -15:1.0344:15 ;-0.3:0.02068:0.3 ; -3.5:0.24:3.5 ; -70:4.8:70] ;
b_critic = [0.1:16.67:500] ;  %30 so far best working 
nodecritic = 30 ;

%
sizes = simsizes;
sizes.NumContStates  = 2*nodecritic;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 6;
sizes.NumInputs      = 19;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1;   % at least one sample time is needed
sys = simsizes(sizes);
x0  = 2*nodecritic;
str = [];
ts  = [0 0];

end
function sys=mdlDerivatives(t,x,u)
global c_critic b_critic nodecritic alpha qm1 qm2 dqm1 dqm2 J
 qm1 = 1.13 ; qm2 = 0.17 ; dqm1= 6.65; dqm2  =3.68;

alpha = 0.97 ; 
q1=u(1);
d_q1=u(2);
q2=u(3);
d_q2=u(4);
f1 = u(6); 
f2 = u(7);
f = [f1;f2] ; 
e1 = u(8); 
e2 = u(9);
de1 = u(10);
de2 = u(11);
J1old =  u(12);
J2old = u(13);
qd1 = u(14);
d_qd1 = u(15);
dd_qd1 = u(16);
qd2 = u(17);
d_qd2 = u(18);
dd_qd2 = u(19);
e =[e1;e2] ; 
de = [de1;de2] ;
qd=[qd1;qd2];
dqd=[d_qd1;d_qd2];
f = [f1;f2] ; 
J_Prev = [J1old;J2old] ; 
r = 0.01*(e+(6*de));

z_critic = [e(1);de(1);qd(1);dqd(1);f(1);e(2);de(2);qd(2);dqd(2);f(2)]
for j=1:1:nodecritic
   h_critic(j)=exp(-norm(z_critic-c_critic(:,j))^2/(b_critic(j)*b_critic(j)));
end

S = [-0.5*((qd1-q1)/qm1)^2 -0.5*((d_qd1-d_q1)/dqm1)^2 ; ...
     -0.5*((qd2-q2)/qm2)^2 -0.5*((d_qd2-d_q2)/dqm2)^2 ]
eC = [J_Prev - S] - alpha.*J 
for i=1:1:nodecritic
sys(i)=alpha*25*h_critic(i)*r(1)-5*x(1)*norm(eC(1));
sys(i+nodecritic)=alpha*25*h_critic(i)*r(2)-5*x(2)*norm(eC(2)); 
end

end
% end mdlDerivatives

function sys=mdlOutputs(t,x,u)
global c_critic b_critic nodecritic qm1 qm2 alpha dqm1 dqm1 dqm2 alpha J
alpha = 0.97 ; 
 qm1 = 1.13 ; qm2 = 0.17 ; dqm1= 6.65; dqm2  =3.68;% maxmimum angles
q1=u(1);
d_q1=u(2);
q2=u(3);
d_q2=u(4);
f1 = u(6); 
f2 = u(7);

e1 = u(8); 
e2 = u(9);
de1 = u(10);
de2 = u(11);
 
J1old =  u(12);
J2old = u(13);
qd1 = u(14);
d_qd1 = u(15);
dd_qd1 = u(16);
qd2 = u(17);
d_qd2 = u(18);
dd_qd2 = u(19);
e =[e1;e2] ; 
de = [de1;de2] ;
qd=[qd1;qd2];
dqd=[d_qd1;d_qd2];
f = [f1;f2] ; 
J_Prev = [J1old;J2old] 
Wcritic_J1=[x(1:nodecritic)]';
Wcritic_J2=[x(nodecritic+1:nodecritic*2)]'
S = [(-0.5*(e1/qm1)^2)-(0.5*(de1/dqm1)^2) ; (-0.5*(e2/qm2)^2)-(0.5*(de2/dqm2)^2)]
    
z_critic = [e(1);de(1);qd(1);dqd(1);f(1);e(2);de(2);qd(2);dqd(2);f(2)]
for j=1:1:nodecritic
 h_critic(j)=exp(-norm(z_critic-c_critic(:,j))^2/(b_critic(j)*b_critic(j)))
end
% fn=[W_f1*h1';
% W_f2*h2'];
J=[Wcritic_J1*h_critic';Wcritic_J2*h_critic'];

eC =  (J_Prev - S) - (alpha.*J)
eA = [0;0] - J  


sys(1)=J(1);
sys(2)=J(2);
sys(3)=eC(1);
sys(4)=eC(2);
sys(5)=eA(1);
sys(6)=eA(2);

end

