close all; clc;

figure(1);
subplot(211);
plot(out.tout,out.qd.Data(:,1),'r',out.tout,out.q.Data(:,1),'k:','linewidth',2);
xlabel('time(s)');ylabel('Angle tracking for knee (link 1)');
legend('ideal angle for knee(link 1)','angle tracking for (link 1)');
subplot(212);
plot(out.tout,out.qd.Data(:,4),'r',out.tout,out.q.Data(:,3),'k:','linewidth',2);
xlabel('time(s)');ylabel('Angle tracking for link 2');
legend('ideal angle for ankle_(link 2)','angle tracking for ankle(link 2)');



