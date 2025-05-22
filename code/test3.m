clc;clear all;

%原始股：
epsilon=0.1; w0=10; k2=178.0715; kesi=6.1938; f1=10; dp1=0; di1=0; p=0.6; dd=0; q=1.5; k=-63.31; tau=0;

%时滞分数阶PID
% epsilon=0.1; w0=10; k2=178.0715; kesi=6.1938; f=10; dp=60; di=60; p=0.8; dd=80; q=1.1; k=-63.31; tau=0.06;

W=0.01:0.01:3;
L1=length(W);

for i=1:L1
    w(i)=w0*W(i);
    sigma(i)=w0^2*(W(i)^2-1)/epsilon;
    a11(i)=di1*w(i)^(-p)*sin(p*pi/2+w(i)*tau)-dd*w(i)^q*sin(q*pi/2-w(i)*tau)+dp1*sin(w(i)*tau)-kesi*w(i)*k2^2/(kesi^2*w(i)^2+(k+k2-w(i)^2)^2);
    a12(i)=di1*w(i)^(-p)*cos(p*pi/2+w(i)*tau)+dd*w(i)^q*cos(q*pi/2-w(i)*tau)+dp1*cos(w(i)*tau)+k2-sigma(i)-k2^2*(k+k2-w(i)^2)/(kesi^2*w(i)^2+(k+k2-w(i)^2)^2);
    a1(i)=100*f1/(a11(i)^2+a12(i)^2)^0.5;
end

%分数阶PID
epsilon=0.1; w0=10; k2=178.0715; kesi=6.1938; f=10; dp2=80; di2=71.3425; p2=1; dd2=80; q2=1; k=-63.31; tau2=0;

W=0.01:0.01:3;
L1=length(W);

for i=1:L1
    w(i)=w0*W(i);
    sigma(i)=w0^2*(W(i)^2-1)/epsilon;
    a11(i)=di2*w(i)^(-p2)*sin(p2*pi/2+w(i)*tau2)-dd2*w(i)^q2*sin(q2*pi/2-w(i)*tau2)+dp2*sin(w(i)*tau2)-kesi*w(i)*k2^2/(kesi^2*w(i)^2+(k+k2-w(i)^2)^2);
    a12(i)=di2*w(i)^(-p2)*cos(p2*pi/2+w(i)*tau2)+dd2*w(i)^q2*cos(q2*pi/2-w(i)*tau2)+dp2*cos(w(i)*tau2)+k2-sigma(i)-k2^2*(k+k2-w(i)^2)/(kesi^2*w(i)^2+(k+k2-w(i)^2)^2);
    a2(i)=100*f/(a11(i)^2+a12(i)^2)^0.5;
end

%分数阶PID
epsilon=0.1; w0=10; k2=178.0715; kesi=6.1938; f=10; dp3=73.9937; di3=69.8127; p3=0.0012; dd3=75.4832; q3=0.9596; k=-63.31; tau3=0;

W=0.01:0.01:3;
L1=length(W);

for i=1:L1
    w(i)=w0*W(i);
    sigma(i)=w0^2*(W(i)^2-1)/epsilon;
    a11(i)=di3*w(i)^(-p3)*sin(p3*pi/2+w(i)*tau3)-dd3*w(i)^q3*sin(q3*pi/2-w(i)*tau3)+dp3*sin(w(i)*tau3)-kesi*w(i)*k2^2/(kesi^2*w(i)^2+(k+k2-w(i)^2)^2);
    a12(i)=di3*w(i)^(-p3)*cos(p3*pi/2+w(i)*tau3)+dd3*w(i)^q3*cos(q3*pi/2-w(i)*tau3)+dp3*cos(w(i)*tau3)+k2-sigma(i)-k2^2*(k+k2-w(i)^2)/(kesi^2*w(i)^2+(k+k2-w(i)^2)^2);
    a3(i)=100*f/(a11(i)^2+a12(i)^2)^0.5;
end

%时滞分数阶PID
epsilon=0.1; w0=10; k2=178.0715; kesi=6.1938; f=10; dp4=76.7567; di4=65.6025; p4=0.39409; dd4=78.0157; q4=1.1646; k=-63.31; tau4=0.067315;

W=0.01:0.01:3;
L1=length(W);

for i=1:L1
    w(i)=w0*W(i);
    sigma(i)=w0^2*(W(i)^2-1)/epsilon;
    a11(i)=di4*w(i)^(-p4)*sin(p4*pi/2+w(i)*tau4)-dd4*w(i)^q4*sin(q4*pi/2-w(i)*tau4)+dp4*sin(w(i)*tau4)-kesi*w(i)*k2^2/(kesi^2*w(i)^2+(k+k2-w(i)^2)^2);
    a12(i)=di4*w(i)^(-p4)*cos(p4*pi/2+w(i)*tau4)+dd4*w(i)^q4*cos(q4*pi/2-w(i)*tau4)+dp4*cos(w(i)*tau4)+k2-sigma(i)-k2^2*(k+k2-w(i)^2)/(kesi^2*w(i)^2+(k+k2-w(i)^2)^2);
    a4(i)=100*f/(a11(i)^2+a12(i)^2)^0.5;
end

plot(W, a1,'color',[1.00,0.00,0.00])
hold on
plot(W, a2,'color',[0.00,0.00,1.00])
hold on
plot(W, a3,'color',[0.49,0.18,0.56])
hold on
plot(W, a4,'color',[0.93,0.69,0.13])
xlabel('\it\Omega','FontSize',15,'FontName','Times New Noman')
ylabel('\itT_A','FontSize',15,'FontName','Times New Roman')
legend('Passive', 'PID', 'FOPID', 'DFOPID')
axis([0 3 0 3])