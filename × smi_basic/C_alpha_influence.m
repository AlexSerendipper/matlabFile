 % figure out the influence of different a
 f = @(phiF,C,a) phiF+C.*sin(phiF+atan(a));  %solve the phi0
 subplot(2,1,1);
 C = 6; 
 a = 2;  %a means the alpha
 phiF=0:0.1:30;
 phi0=f(phiF,C,a);
 plot(phi0,phiF);
 hold on;
 C = 6;
 a = 5;
 phiF = 0:0.1:30;
 phi0=f(phiF,C,a);
 plot(phi0,phiF,':')
 hold on;
 C = 6; 
 a = 8;
 phiF = 0:0.1:30;
 phi0 = f(phiF,C,a);
 plot(phi0,phiF,'-.');
 legend('c=6,a=2','c=6,a=5','c=6,a=8')
 
 % figure out the influence of different C
 subplot(2,1,2);
 C = 9;
 a = 3;
 phiF = 0:0.1:30;
 phi0 = f(phiF,C,a);
 plot(phi0,phiF);
 hold on;
 C = 6;
 a = 3;
 phiF = 0:0.1:30;
 phi0 = f(phiF,C,a);
 plot(phi0,phiF,':');
 hold on;
 C = 4;
 a = 3;
 phiF = 0:0.1:30;
 phi0 = f(phiF,C,a);
 plot(phi0,phiF,'-.');
 hold on;
 C = 0.8;
 a = 3;
 phiF = 0:0.1:30;
 phi0 = f(phiF,C,a);
 plot(phi0,phiF,'--');
 hold on;
 legend('c=9,a=3','c=6,a=3','c=4,a=3','c=0.8,a=3');
 

 
 