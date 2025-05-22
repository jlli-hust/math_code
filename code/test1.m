clc;clear all;
epsilon=0.1; w0=10; k2=60; kesi=2; f=10; dp=10; di=10; p=0.6; dd=10; q=1.5; k=-10; W=1; w=w0*W; sigma=w0^2*(W^2-1)/epsilon;
h=0.01;
t=0:h:300;
tau=0.1; 
%下为数值结果
if tau==0
    m=1;
else
    if mod(tau,h)==0
        m=round(tau/h);
    else
        m=1+round(tau/h);
    end
end

y1(m)=1;
y2(m)=0;
y3(m)=0;
y4(m)=0;
y5(m)=0;
y6(m)=0;
Cl(1)=1;
Cll(1)=1;

for i=m:length(t)-1
    i
    y1(i+1)=y3(i)*h+y1(i);
    y2(i+1)=y4(i)*h+y2(i);
    y3(i+1)=(epsilon*f*cos(w*i*h)-w0^2*y1(i+1)-epsilon*k2*(y1(i+1)-y2(i+1))-epsilon*dp*y1(i-m+1)-epsilon*di*y5(i-m+1)-epsilon*dd*y6(i-m+1))*h+y3(i);
    y4(i+1)=(-k*y2(i+1)-k2*(y2(i+1)-y1(i+1))-kesi*y4(i))*h+y4(i);
    vl=0;
    vll=0;
    for j=1:i
        Cl(j+1) = (1-(1+p)/j)*Cl(j);
        vl=vl+Cl(j+1)*y5(i+1-j);
        j=j+1;
    end
    y5(i+1)=y1(i)*h^p-vl;
    for j=1:i
        Cll(j+1) = (1-(1+0.5)/j)*Cll(j);
        vll=vll+Cll(j+1)*y6(i+1-j);
        j=j+1;
    end
%     y6(i+1)=y3(i)*h^0.5-vll;
    y6(i+1)=(epsilon*f*cos(w*i*h)-w0^2*y1(i+1)-epsilon*k2*(y1(i+1)-y2(i+1))-epsilon*dp*y1(i-m+1)-epsilon*di*y5(i-m+1)-epsilon*dd*y6(i-m+1))*h^0.5-vll;
end

%下为解析结果
a1=f/((di*w^(-p)*sin(p*pi/2+w*tau)-dd*w^q*sin(q*pi/2-w*tau)+dp*sin(w*tau)-kesi*w*k2^2/(kesi^2*w^2+(k+k2-w^2)^2))^2+(di*w^(-p)*cos(p*pi/2+w*tau)+dd*w^q*cos(q*pi/2-w*tau)+dp*cos(w*tau)+k2-sigma-k2^2*(k+k2-w^2)/(kesi^2*w^2+(k+k2-w^2)^2))^2)^0.5;
theta1=atan((di*w^(-p)*sin(p*pi/2+w*tau)-dd*w^q*sin(q*pi/2-w*tau)+dp*sin(w*tau)-kesi*w*k2^2/(kesi^2*w^2+(k+k2-w^2)^2))/(di*w^(-p)*cos(p*pi/2+w*tau)+dd*w^q*cos(q*pi/2-w*tau)+dp*cos(w*tau)+k2-sigma-k2^2*(k+k2-w^2)/(kesi^2*w^2+(k+k2-w^2)^2)));
% theta1=3.2+atan((di*w^(-p)*sin(p*pi/2+w*tau)-dd*w^q*sin(q*pi/2-w*tau)+dp*sin(w*tau)-kesi*w*k2^2/(kesi^2*w^2+(k+k2-w^2)^2))/(di*w^(-p)*cos(p*pi/2+w*tau)+dd*w^q*cos(q*pi/2-w*tau)+dp*cos(w*tau)+k2-sigma-k2^2*(k+k2-w^2)/(kesi^2*w^2+(k+k2-w^2)^2)));

t1=0:0.01:300;
L1=length(t1);

for i=1:L1
y11(i)=a1*cos(w*t1(i)+theta1);
end

% plot(t,y1);
step=2
t=downsample(t,step)
y111= downsample(y1,step)
plot(t,y111,'*','MarkerEdgeColor','b');
hold on
plot(t1, y11)
xlabel('\itt','FontSize',15,'FontName','Times New Noman')
ylabel('\itx_1','FontSize',15,'FontName','Times New Roman')
set(gca,'xtick',290:2:300)
set(gca,'ytick',-0.05:0.025:0.05)
legend('Num. sol.','App. Anal. sol.')
axis([290 300 -0.05 0.05])



