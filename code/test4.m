clc;clear all;

%FOPID
epsilon=0.05; w0=10; k2=178.0715; kesi=6.1938; f=10; dp=76.7567; di=65.6025; p=0.39; dd=78.0157; q=1.1646; k=-63.31; W=1; w=w0*W; sigma=w0^2*(W^2-1)/epsilon; tau=0.067315;

h=0.01;
t=0:h:100; 

N=length(t);
randn('seed',N);
Y=randn(N,1);
b=fix((t)*100)+1;


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

y1(m)=0;
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
    y3(i+1)=(Y(fix((i*h)*100)+1)-w0^2*y1(i+1)-epsilon*k2*(y1(i+1)-y2(i+1))-epsilon*dp*y1(i-m+1)-epsilon*di*y5(i-m+1)-epsilon*dd*y6(i-m+1))*h+y3(i);
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
    y6(i+1)=(Y(fix((i*h)*100)+1)-w0^2*y1(i+1)-epsilon*k2*(y1(i+1)-y2(i+1))-epsilon*dp*y1(i-m+1)-epsilon*di*y5(i-m+1)-epsilon*dd*y6(i-m+1))*h^0.5-vll;
end

%uncontrol
w01=10;

h1=0.01;
t1=0:h1:100; 

N1=length(t1);
randn('seed',N1);
Y1=randn(N1,1);
b1=fix((t1)*100)+1;





%下为数值结果
y11(1)=0;
y31(1)=0;

for i=1:length(t1)-1
    i
    y11(i+1)=y31(i)*h1+y11(i);
    y31(i+1)=(Y1(fix((i*h1)*100)+1)-w01^2*y11(i+1))*h1+y31(i);
end


plot(t1,y11);
hold on
plot(t,y1,'color',[1.00,0.00,0.00]);
xlabel('\itt','FontSize',15,'FontName','Times New Noman')
ylabel('\itx_1','FontSize',15,'FontName','Times New Roman')
set(gca,'xtick',0:10:100)
set(gca,'ytick',-0.12:0.04:0.12)
axis([0 100 -0.12 0.12])
legend('Uncontrol', 'DFOPID')

