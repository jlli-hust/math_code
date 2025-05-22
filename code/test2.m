clc;clear all;

%原始股：
% epsilon=0.1; w0=10; k2=178.0715; kesi=6.1938; f=10; dp=0; di=0; p=0.6; dd=0; g=1.5; k=-63.31; tau=0;

%时滞分数阶PID
epsilon=0.1; w0=10; k2=60; kesi=2; f=10; dp=10; di=60; dd=10; q=1.5; k=-10; tau=0.1;

W=0.15:0.005:3;
L1=length(W);

p1=0.1;
for i=1:L1
    w(i)=w0*W(i);
    sigma(i)=w0^2*(W(i)^2-1)/epsilon;
    a11(i)=di*w(i)^(-p1)*sin(p1*pi/2+w(i)*tau)-dd*w(i)^q*sin(q*pi/2-w(i)*tau)+dp*sin(w(i)*tau)-kesi*w(i)*k2^2/(kesi^2*w(i)^2+(k+k2-w(i)^2)^2);
    a12(i)=di*w(i)^(-p1)*cos(p1*pi/2+w(i)*tau)+dd*w(i)^q*cos(q*pi/2-w(i)*tau)+dp*cos(w(i)*tau)+k2-sigma(i)-k2^2*(k+k2-w(i)^2)/(kesi^2*w(i)^2+(k+k2-w(i)^2)^2);
    a1(i)=100*f/(a11(i)^2+a12(i)^2)^0.5;
end

p2=0.3;
for i=1:L1
    w(i)=w0*W(i);
    sigma(i)=w0^2*(W(i)^2-1)/epsilon;
    a11(i)=di*w(i)^(-p2)*sin(p2*pi/2+w(i)*tau)-dd*w(i)^q*sin(q*pi/2-w(i)*tau)+dp*sin(w(i)*tau)-kesi*w(i)*k2^2/(kesi^2*w(i)^2+(k+k2-w(i)^2)^2);
    a12(i)=di*w(i)^(-p2)*cos(p2*pi/2+w(i)*tau)+dd*w(i)^q*cos(q*pi/2-w(i)*tau)+dp*cos(w(i)*tau)+k2-sigma(i)-k2^2*(k+k2-w(i)^2)/(kesi^2*w(i)^2+(k+k2-w(i)^2)^2);
    a2(i)=100*f/(a11(i)^2+a12(i)^2)^0.5;
end

p3=0.5;
for i=1:L1
    w(i)=w0*W(i);
    sigma(i)=w0^2*(W(i)^2-1)/epsilon;
    a11(i)=di*w(i)^(-p3)*sin(p3*pi/2+w(i)*tau)-dd*w(i)^q*sin(q*pi/2-w(i)*tau)+dp*sin(w(i)*tau)-kesi*w(i)*k2^2/(kesi^2*w(i)^2+(k+k2-w(i)^2)^2);
    a12(i)=di*w(i)^(-p3)*cos(p3*pi/2+w(i)*tau)+dd*w(i)^q*cos(q*pi/2-w(i)*tau)+dp*cos(w(i)*tau)+k2-sigma(i)-k2^2*(k+k2-w(i)^2)/(kesi^2*w(i)^2+(k+k2-w(i)^2)^2);
    a3(i)=100*f/(a11(i)^2+a12(i)^2)^0.5;
end

p4=0.7;
for i=1:L1
    w(i)=w0*W(i);
    sigma(i)=w0^2*(W(i)^2-1)/epsilon;
    a11(i)=di*w(i)^(-p4)*sin(p4*pi/2+w(i)*tau)-dd*w(i)^q*sin(q*pi/2-w(i)*tau)+dp*sin(w(i)*tau)-kesi*w(i)*k2^2/(kesi^2*w(i)^2+(k+k2-w(i)^2)^2);
    a12(i)=di*w(i)^(-p4)*cos(p4*pi/2+w(i)*tau)+dd*w(i)^q*cos(q*pi/2-w(i)*tau)+dp*cos(w(i)*tau)+k2-sigma(i)-k2^2*(k+k2-w(i)^2)/(kesi^2*w(i)^2+(k+k2-w(i)^2)^2);
    a4(i)=100*f/(a11(i)^2+a12(i)^2)^0.5;
end

p5=0.9;
for i=1:L1
    w(i)=w0*W(i);
    sigma(i)=w0^2*(W(i)^2-1)/epsilon;
    a11(i)=di*w(i)^(-p5)*sin(p5*pi/2+w(i)*tau)-dd*w(i)^q*sin(q*pi/2-w(i)*tau)+dp*sin(w(i)*tau)-kesi*w(i)*k2^2/(kesi^2*w(i)^2+(k+k2-w(i)^2)^2);
    a12(i)=di*w(i)^(-p5)*cos(p5*pi/2+w(i)*tau)+dd*w(i)^q*cos(q*pi/2-w(i)*tau)+dp*cos(w(i)*tau)+k2-sigma(i)-k2^2*(k+k2-w(i)^2)/(kesi^2*w(i)^2+(k+k2-w(i)^2)^2);
    a5(i)=100*f/(a11(i)^2+a12(i)^2)^0.5;
end

plot(W, a1,'color',[1.00,0.00,0.00])
hold on
plot(W, a2,'color',[0.00,0.00,1.00])
hold on
plot(W, a3,'color',[0.49,0.18,0.56])
hold on
plot(W, a4,'color',[0.93,0.69,0.13])
hold on
plot(W, a5,'color',[0.30,0.75,0.93])
xlabel('\it\Omega','FontSize',15,'FontName','Times New Noman')
ylabel('\itT_A','FontSize',15,'FontName','Times New Roman')
legend('\lambda=0.1', '\lambda=0.3', '\lambda=0.5', '\lambda=0.7', '\lambda=0.9')
set(gca,'xtick',0:0.5:3)
set(gca,'ytick',0:0.5:3.5)
axis([0 3 0 3.5])