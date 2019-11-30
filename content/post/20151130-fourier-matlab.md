---
date: "2015-08-03T21:13:14-05:00"
title: Fourier Matlab script
---

```matlab
clc, close all

% load('fourier_data.mat')
% load('fourier_qiang.mat')
%%  fourier main script
%   Writed By Dong 2014/07/22

% x=[1 1.4 1.0 1.4 1 1.4 1.0 1.4];
% x=[0 3 2 0 -1];
% x=[1 1 1 1 0 0 0 0];
x=[3.304 2.5535 -0.22396 0.56711 0.06798 0.69711 -0.14572 0.7031, ...
    -1.635 -1.1464 -0.55064 -0.70925 -0.18781 -0.55497 -0.62887 -2.1101];
% x=x_resi;
% x=y;

n=length(x);

% Fs=n/(2*pi);
Fs = 12;

t  = (0:n-1)/Fs;
tt = (0:0.1:n)/Fs;

M=floor(n/2);
d=fft(x);

a0 = 2*d(1)/n;
an = 2*real(d(2:M+1))/n;
bn = -2*imag(d(2:M+1))/n;
if n==2*M,an(end)=an(end)/2;end

[an;bn];

%% 主频率分析
f = (0:n-1)*(Fs/n);     % Frequency range
power = d.*conj(d)/n;   % Power of the DFT
figure('pos',[267   149   740   444])

subplot(211)

f=f(2:M+1);
power=power(2:M+1);
plot(f,power)

xlabel('Frequency (Hz)')
ylabel('Power')
title('{\bf Periodogram}')


[val,index]=sort(power,'descend');
subplot(212)

plot(cumsum(val)/sum(val),'bo','markersize',3,'markerfaceColor','b')
xlabel('harmonic wave number'),ylabel('Variance contribution')
set(gca,'ygrid','on')

%%  拟合结果
result=a0/2;
for i=1:M
    result=result+an(i)*cos(i*tt*2*pi/n*Fs)+bn(i)*sin(i*tt*2*pi/n*Fs);
end

%%  去除较为明显的三个周期性成分
result2=a0/2;
nWave=3;        %三个主要简谐波
f=f(index);T=1./f(1:nWave);

T               %周期
for i=1:3
    k=index(i);
    result2=result2+an(k)*cos(k*t*2*pi/n*Fs)+bn(k)*sin(k*t*2*pi/n*Fs);
end

%%  图件输出
figure('pos',[267   149   740   444],'numbertitle','off','name','Fourier')
% t=(0:n-1)*2*pi/n;
subplot(211)

handle=plot(t,x,'r.','markersize',6,'markerfaceColor','r');hold on
plot(t,result2,'b-','linewidth',2)

handle_legend=legend('Original Data','Fourier Curve','location','NW');
set(handle_legend,'fontsize',8)

uistack(handle_legend,'down',2)
uistack(handle,'top');
set(gca,'fontname','微软雅黑','fontsize',12)

y_lim=get(gca,'ylim');

grid on
title('Fourier Transform','fontname','微软雅黑','fontsize',12)

subplot(212)
fourier_resi=x'-result2;
plot(t,fourier_resi)
title('Residual','fontname','微软雅黑','fontsize',12)
```

``` c
#include <stdio.h>
int main(void)
{
  printf("Hello world!");
  return 0;
}
```

```r
wave.3 <- 0.5 * wave.1 + 0.25 * wave.2
plot(xs,wave.3,type="l") 
title("Eg complex wave")
abline(h=0,lty=3)

# repeat.xs     <- seq(-2*pi,0,pi/100)
wave.3.repeat <- 0.5*sin(3*repeat.xs) + 0.25*sin(10*repeat.xs)
# plot(xs,wave.3,type="l") 
title("Repeating pattern")
points(repeat.xs,wave.3.repeat,type="l",col="red") 
abline(h=0,v=c(-2*pi,0),lty=3)
```
