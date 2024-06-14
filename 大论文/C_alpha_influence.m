 clc;
 clear all;
 close all;
 % figure out the influence of different a
 f = @(phiF,C,a) phiF+C.*sin(phiF+atan(a));  %solve the phi0
 subplot(2,1,1);
 C = 2.7; 
 a1 = 2;  % a means the alpha
 a2 = 4;
 a3 = 7;
 phiF=0:0.1:30;
 phi01=f(phiF,C,a1);
 phi02=f(phiF,C,a2);
 phi03=f(phiF,C,a3); 
 
 plot(phi01,phiF);
 hold on;
 plot(phi02,phiF,':')
 hold on;
 plot(phi03,phiF,'-.');
 legend('c=0.7,a=2','c=0.7,a=4','c=0.7,a=7')
 
 % figure out the influence of different C
 subplot(2,1,2);
 C1 = 0.7;
 C2 = 2.7;
 C3 = 4.7;
 C4 = 6.7;
 a = 5;
 phiF = 0:0.1:30;
 
 phi01 = f(phiF,C1,a);
 plot(phi01,phiF);
 hold on;
 
 phi02 = f(phiF,C2,a);
 plot(phi02,phiF,':');
 hold on;
 
 phi03 = f(phiF,C3,a);
 plot(phi03,phiF,'-.');
 hold on;
 
 phi04 = f(phiF,C4,a);
 plot(phi04,phiF,'--');
 hold on;
 legend('c=9,a=3','c=6,a=3','c=4,a=3','c=0.8,a=3');
 
 
 
 