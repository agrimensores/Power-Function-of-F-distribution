function p=powerf(alfa,h,f,NP)
%%
%powerf :yields the power of the test for statistics following
%        F-distribution (complementary of CDF of non-central F-distribution)
%--------------------------------------------------------------------------
%OUTPUT :
%       p    = Power of the test
%INPUT  :
%       alfa = Significance level
%       h    = First degrees of freedom
%       f    = Second degrees of freedom
%       NP   = Non-centrality parameter
%
%--------------------------------------------------------------------------
%Written by Cuneyt Aydin (2020-2021, 2023)
%       Yildiz Technical University-Civil Engineering Faculty
%       Geodesy Division, Istanbul/Turkey
%       Contact: caydin78@gmail.com;caydin@yildiz.edu.tr
%       https://avesis.yildiz.edu.tr/caydin
%--------------------------------------------------------------------------
%References:
%[1] Aydin C and Gunes O (2024). Power function of F-distribution: 
% Revisiting its computation and solution for geodetic studies, 
% Journal of Geodesy, 98, https://doi.org/10.1007/s00190-024-01905-7
%--------------------------------------------------------------------------

%%

if f>10^12 %this is for saving time. To see whether it works even for f>10^10, remove this line
    f=10^12; 
end
done=0;i=0;epsL=10^-20;
F=f_inv(1-alfa,h,f);
delta=NP/2;a=h/2;b=f/2;x=h*F/(f+h*F);
T=1-alfa;logT=log(T);m=1-alfa;%the first term
sumP=exp(logT-delta);

while ~done
    i=i+1;
    c=betaa(x,a+i-1,b);
    m=m-c;
    if m>0
        logT=logT+log(delta)+log(m)-log(m+c)-log(i);
        P=exp(logT-delta);
        sumP=sumP+P;
        if i>delta
            if P<epsL
                done=1;
            end
        end
    end
    if m<0
        done=1;
    end
end
p=exp(log(sumP)); %CDF
p=1-p; %power
%%
function c=betaa(x,a,b)
%logarithm of gamma(a+b)/(gamma(a)*gamma(b))*x^a*(1-x)^b
d=(a)*log(x)+b*log(1-x)+gln(a+b)-gln(a+1)-gln(b);
c=exp(d);
%%
function g=gln(s)
%ln of gamma function (Stirling Approximation-Abromowitz and Stegun-1972, p.
%257)
%Didonato and Jarnagin (1967)
if s>100
    g=(s-0.5)*log(s-1)-(s-1)+0.5*log(2*pi)+1/(12*(s-1))-1/(360*(s-1)^3)+1/(1260*(s-1)^5);
else
    h=mod(s,1);
    if h==0.5
        s=s-0.5;g=log(sqrt(pi))-s*log(2);
        for i=1:s
            g=g+log(2*s-2*i+1);
        end
    end
    if h==0
        g=0;
        for i=1:(s-1)
            g=g+log(i);
        end
    end
end
%%
function F=f_inv(CL,h,f)
%Upper alpha percentage value of F-distribution
%CL=1-alpha
%h=first degrees of freedom
%f=second degrees of freedom
%Note=The program is stable for CL<0.99 and f>1. For CL getting close to 1
%the percentage value go to infitinity, and so big errors may occur. However,
%these errors do not affect on the power of the test which is the main aim 
%of the program "powerf"

Fa=0;Fb=10^10;eps2=10^-12;
if CL>0.99
    eps2=10^-5;
end
Fo=(Fa+Fb)/2;done=0;
while done==0
    c=cdf0(Fo,h,f);
    if c<CL
        Fa=Fo;
    else
        Fb=Fo;
    end
    Fo=(Fa+Fb)/2;
    if abs(Fa-Fb)<eps2
        done=1;F=Fo;
    end
end
%%
function p=cdf0(Fo,h,f)
%CDF of central F-distribution
%Fo=value
%h=first degrees of freedom
%f=second degrees of freedom

a=h/2;b=f/2;x=h*Fo/(f+h*Fo);

if x*(a+b+2)<(a+1)
    p=icbeta(x,a,b); %incomplete beta function ratio(Zelen and Severo Algorithm)
else
    p=1-icbeta(1-x,b,a);
end
%%
function pp=icbeta(x,p,q)

%icbeta:  Incomplete beta function ratio while x*(p+q-2)<(p-1).
%         Otherwise, pp=1-icbeta(1-x,q,p). 
%
%         This function is created from Phien (1990) for computing the upper 
%         percentage value of F-distribution. It gives the most efficient 
%         solution among many algorithms for incomplete beta function.
%
%Created: by C. Aydin (on April-2021)
%         Yildiz Technical University, Geodesy Division
%         Istanbul, Turkey
%         caydin78@gmail.com;caydin@yildiz.edu.tr
%         https://avesis.yildiz.edu.tr/caydin
%
%References:
%     [1] Phien, HN, (1990). A note on the computation of the Incomplete
%         beta function, Adv. Eng. Software, 12/1, pp.39-44.
%     [2] Abramowitz, M and Stegun, IA, (1972). Handbook of mathematical
%         functions, Dover Publications, New York.

eps1=10^-12;
An=1;Bn=1;Az=1;

qab=p+q;qap=p+1;qan=p-1;
Bz=1-qab*x/qap;
done=0;n=0;
while done==0
    n=n+1;
    D=n*(q-n)*x/((qan+2*n)*(p+2*n));
    Ap=Az+D*An;
    Bp=Bz+D*Bn;
    D=-(p+n)*(qab+n)*x;
    D=D/((p+2*n)*(qap+2*n));
    App=Ap+D*Az;
    Bpp=Bp+D*Bz;
    Aold=Az;
    An=Ap/Bpp;Bn=Bp/Bpp;
    Az=App/Bpp;Bz=1;
    if abs(Az-Aold)<eps1*abs(Az)
       done=1;pp=Az;
    end
end
Bab=gln(p)+gln(q)-gln(p+q);
c=exp(p*log(x)+q*log(1-x)-Bab)/p; %constant

pp=pp*c;
