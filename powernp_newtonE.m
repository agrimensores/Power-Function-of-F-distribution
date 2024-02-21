function np=powernp_newtonE(alfa,power,h,f)
%%%
%powernp :  Gives the non-centrality parameter for a power function 
%           with (alfa, power, h and f)--Newton-Raphson Modified algorithm
%--------------------------------------------------------------------------
%OUTPUT  :
%           np      = Lower Bound of the Non-Centrality Parameter
%INPUT   :
%           alfa    = Significance level
%           power   = Power of the test
%           h       = First degrees of freedom (numerator)
%           f       = Second degrees of freedom (denominator)
%
%NOTE      : Take f=10^8 for computing the noncentrality parameter of
%            chi-square distribution
%--------------------------------------------------------------------------
%Written by Cuneyt Aydin (2021-2023)
%           Yildiz Technical University-Civil Engineering Faculty
%           Geodesy Division, Istanbul/Turkey
%           Contact: caydin78@gmail.com;caydin@yildiz.edu.tr
%           https://avesis.yildiz.edu.tr/caydin
%--------------------------------------------------------------------------
%References:
%[1] Aydin C and Gunes O (2023), Power function of F-distribution:
%Revisiting its computation and solution for geodetic studies.
%[2] Guirgus GH (1990), A note on computing the noncentrality 
%parameter of the noncentral F-distribution, Communications in Statistics
%-Simulation and Computation,19:4, 1497-1511.
%[3] Ding CG (1997),On using newton's method for computing the
%noncentrality parameter of the noncentral F distribution, Communications 
%in Statistics - Simulation and Computation, 26:1, 259-26.
%--------------------------------------------------------------------------


if power<alfa
    error('power should be bigger than or equal to alfa')
end
if (alfa<1*10^-5)||(power>0.9999)
    error('alfa should be bigger than 10^-5, power should be smaller than 0.9999')
end
if (h<=0)||(f<=0)
    error('h and f should be bigger than zero')
end
if (floor(h)-h<0)||(floor(f)-f<0)
    error('h and f should be integers')
end

powerd= @(L)(1-powerf(alfa,h,f,L))-(1-power);

F=f_inv(1-alfa,h,f);
L=LNPapp_c(F,alfa,power,h,f);

if L>10^6
    error('Solution is out of capacity')
end

eps2=10^-6;

while true

    power_dif=powerd(L);
    power_dif_deriv=(powerd(L+eps2)-power_dif)/eps2; %numerical derivative
    L_next=L+((power_dif+1-power)/power_dif_deriv)*log((1-power)/(power_dif+1-power)); %modified Newton-Raphson
    if L_next>10^6
        error('Solution is out of capacity')
    end

    if abs(L_next-L)<eps2
        L=L_next;
        break;
    end
    L=L_next;
end
np=L;
%%
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
%[1] Aydin C and Gunes O (2023), Power function of F-distribution:
%Revisiting its computation and solution for geodetic studies.
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
        T=(T*delta/i)*m/(m+c);
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
function delta=LNPapp_c(c_value,alfa,power,h,f)

% Approximation for non-centrality parameter computation using
% normal-distribution approximation
eps2=10^-5;
la=0;lb=1;

if alfa>power
    error('alfa must be smaller than power')
end

while lb-la>eps2

    lo=(la+lb)/2;
    powerd=ffcdf(c_value,h,f,lo)-(1-power);
    if abs(powerd)<eps2
        break;
    end
    if powerd<0
        lb=lo;
    else
        la=lo;
        if lo>lb*0.9
            lb=lb*2;
        end
    end
end
delta=(la+lb)/2;


%%
function c=ffcdf(x,h,f,L)
%Approximate solution for CDF of non-central F-distribution 
%--------------------------------------------------------------------------
%References:
%[1] Severo NC and Zelen M (1960), Normal approximation to the chi-square
%and non-central F probability functions, Biometrika, 47(3-4),411-416.
%[2] Tiku ML, (1965),Laguerre series forms of non-central chi-square and 
%F distributions, Biometrika,52 (3-4), 415-427.
%[3] Vazquez-Leal H, Castaneda-Sheissa R, Filobello-Nino U, 
%Sarmiento-Reyes A, Sanchez Orea H, (2012), High accurate simple 
%approximation of normal distribution integral, Mathematical Problems in 
%Engineering, vol. 2012, Article ID 124029, 22 pages,
%https://doi.org/10.1155/2012/124029.
%--------------------------------------------------------------------------
%approximation of Severo and Zelen (1960), Biometrica
x1=(h*x/(h+L))^(1/3);
x2=1-2/(9*f);
x4=2*(h+2*L)/(9*(h+L)^2);
x3=1-x4;
x5=(1-x2)*((h*x/(h+L))^(2/3));
xx=(x1*x2-x3)/sqrt(x4+x5);

%approximation of Tiku 1965, Biometrika
%     x1=sqrt(2*f-1);
%     x2=h*x/f;
%     x3=2*(h+L);
%     x4=(h+2*L)/(h+L);
%     xx=(x1*sqrt(x2)-sqrt(x3-x4))/sqrt(x2+x4);

%approximation of cdf of ....mathematical problems in engineering (high
%accurate simple approximation of normal distribution
b1=exp(-358*xx/23+111*atan(37*xx/294));
bbo=1/(b1+1);
c=bbo;
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
%the percentage value go to infinity, and so big errors may occur. However,
%these errors do not affect on the power of the test which is the main aim
%of the program "powerf"

Fa=0;Fb=10^10;Fo=(Fa+Fb)/2;eps2=10^-12;
if CL>0.99
    eps2=10^-5;
end
done=0;
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
    if Fo>10^10
        error('solution is out of capacity')
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
    p=icbeta(x,a,b); %incomplete beta function (Zelen and Severo Algorithm)
else
    p=1-icbeta(1-x,b,a);
end
%%
function pp=icbeta(x,p,q)

%icbeta:  Incomplete beta function while x*(p+q-2)<(p-1).
%         Otherwise, pp=1-icbeta(1-x,q,p).
%
%         This function is created from Phien (1990) for computing the upper
%         percentage value of F-distribution. It gives the most eff_invcient
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

function p=powerf2(alfa,h,f,NP)
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
%[1] Aydin C and Gunes O (2023), Power function of F-distribution:
%Revisiting its computation and solution for geodetic studies.
%--------------------------------------------------------------------------

%%

if f>10^12 %this is for saving time. To see whether it works even for f>10^10, remove this line
    f=10^12; 
end
done=0;i=0;epsL=10^-20;
F=f_inv(1-alfa,h-2,f);
F=F*(h-2)/(h-2+2);
delta=NP/2;a=h/2;b=f/2;x=h*F/(f+h*F);
T=1-alfa;logT=log(T);m=1-alfa;%the first term
sumP=exp(logT-delta);

while ~done
    i=i+1;
    c=betaa(x,a+i-1,b);
    m=m-c;
    if m>0
        logT=logT+log(delta)+log(m)-log(m+c)-log(i);
        T=(T*delta/i)*m/(m+c);
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