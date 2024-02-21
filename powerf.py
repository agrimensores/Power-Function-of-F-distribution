import math

def powerf(alfa, h, f, np):
#
#powerf :yields the power of the test for statistics following
#        F-distribution (complementary of CDF of non-central F-distribution)
#------------------------------------------------------------------------------
#OUTPUT :
#       p    = Power of the test
#INPUT  :
#       alfa = Significance level
#       h    = First degrees of freedom
#       f    = Second degrees of freedom
#       NP   = Non-centrality parameter
#
#------------------------------------------------------------------------------
#Written by Cuneyt Aydin (2020-2021, 2023)
#       Yildiz Technical University-Civil Engineering Faculty
#       Geodesy Division, Istanbul/Turkey
#       Contact: caydin78@gmail.com;caydin@yildiz.edu.tr
#       https://avesis.yildiz.edu.tr/caydin
#------------------------------------------------------------------------------
#References:
#[1] Aydin C and Gunes O (2023), Power function of F-distribution:
#Revisiting its computation and solution for geodetic studies.
#------------------------------------------------------------------------------
    if f>1e12: #this is for saving time. To see whether it works even for f>10^10, remove this line
        f=1e12

    done = 0
    i = 0
    epsL = 1e-20
    ff=f_inv(1-alfa,h,f) #upper alfa percentage value of F-distribution
    if np<1e-20:
        np=1e-24
        
    delta = np / 2
    a = h / 2
    b = f / 2
    x = h * ff / (f + h * ff)
    t = 1 - alfa
    logt = math.log(t)
    m = 1 - alfa  # the first term
    sump = math.exp(logt-delta)
    while done == 0:
        i = i + 1
        c = betaa(x, a + i - 1, b)
        m = m - c

        if m > 0:
            logt = logt + math.log(delta) + math.log(m) - math.log(m + c) - math.log(i)
            pp = math.exp(logt - delta)
            sump = sump + pp
            if i > delta:
                if pp < epsL:
                    done = 1
        if m <= 0:
            done = 1
        if sump < 1e-20:
            sump=1e-20
    p = math.exp(math.log(sump))  # CDF
    p = 1 - p  # power

    return p
#------------------------------------------------------------------------------
def betaa(x, a, b):
#logarithm of gamma(a+b)/(gamma(a)*gamma(b))*x^a*(1-x)^b
    d = (a) * math.log(x) + b * math.log(1 - x) + gln(a + b) - gln(a + 1) - gln(b)
    c = math.exp(d)

    return c
#------------------------------------------------------------------------------

def gln(s):
#ln of gamma function (Stirling Approximation-Abromowitz and Stegun-1972, p.257)
#Didonato and Jarnagin (1967)
    if s > 100:
        g = (s - 0.5) * math.log(s - 1) - (s - 1) + 0.5 * math.log(2 * math.pi) + 1 / (12 * (s - 1)) - 1 / (
                360 * (s - 1) ** 3) + 1 / (1260 * (s - 1) ** 5)
    else:
        h = float(s) - int(s)
        if h == 0.5:
            s = s - 0.5
            g = math.log(math.sqrt(math.pi)) - s * math.log(2)
            for i in range(1, int(s) + 1):
                g = g + math.log(2 * s - 2 * i + 1)
        if h == 0:
            g = 0
            for i in range(1, int(s)):
                g = g + math.log(i)

    return g
#------------------------------------------------------------------------------


def icbeta(x, p, q):
#icbeta:  Incomplete beta function ratio while x*(p+q-2)<(p-1).
#         Otherwise, pp=1-icbeta(1-x,q,p). 
#
#         This function is created from Phien (1990) for computing the upper 
#         percentage value of F-distribution. It gives the most efficient 
#         solution among many algorithms for incomplete beta function.
#
#Created: by C. Aydin (on April-2021)
#         Yildiz Technical University, Geodesy Division
#         Istanbul, Turkey
#         caydin78@gmail.com;caydin@yildiz.edu.tr
#         https://avesis.yildiz.edu.tr/caydin
#
#References:
#     [1] Phien, HN, (1990). A note on the computation of the Incomplete
#         beta function, Adv. Eng. Software, 12/1, pp.39-44.
#     [2] Abramowitz, M and Stegun, IA, (1972). Handbook of mathematical
#         functions, Dover Publications, New York.

    eps1 = 1e-12
    an = 1
    bn = 1
    az = 1
    qab = p + q
    qap = p + 1
    qan = p - 1
    bz = 1 - qab * x / qap
    done = 0
    n = 0
    bab = gln(p) + gln(q) - gln(p + q)
    c = math.exp(p * math.log(x) + q * math.log(1 - x) - bab) / p  # constant
    while done == 0:
        n = n + 1
        d = n * (q - n) * x / ((qan + 2 * n) * (p + 2 * n))
        ap = az + d * an
        bp = bz + d * bn
        d = -(p + n) * (qab + n) * x
        d = d / ((p + 2 * n) * (qap + 2 * n))
        app = ap + d * az
        bpp = bp + d * bz
        aold = az
        an = ap / bpp
        bn = bp / bpp
        az = app / bpp
        bz = 1
        if abs(az - aold) < (eps1 * abs(az)):
            done = 1
            pp = az
    pp = pp * c

    return pp

#------------------------------------------------------------------------------
def cdf0(fo, h, f):
# cdf0: CDF of central F-distribution
#      Fo=value
#      h=first degrees of freedom
#      f=second degrees of freedom
    a = h / 2
    b = f / 2
    x = h * fo / (f + h * fo)
    if x * (a + b + 2) < (a + 1):
        p = icbeta(x, a, b)
    else:
        p = 1 - icbeta(1 - x, b, a)

    return p

#------------------------------------------------------------------------------
def f_inv(cl, h, f):
#Upper alpha percentage value of F-distribution
#CL=1-alpha
#h=first degrees of freedom
#f=second degrees of freedom
#Note=The program is stable for CL<0.99 and f>1. For CL getting close to 1
#the percentage value go to infitinity, and so big errors may occur. However,
#these errors do not affect on the power of the test which is the main aim 
#of the program "powerf"
    fa = 0
    fb = 1e10
    eps2 = 1e-12
    if cl > 0.99:
        eps2 = 1e-8
    fo = (fa + fb) / 2
    done = 0
    while done == 0:
        c = cdf0(fo, h, f)
        if c < cl:
            fa = fo
        else:
            fb = fo
        fo = (fa + fb) / 2
        if abs(fa - fb) < eps2:
            done = 1
            ff = fo

    return ff
#------------------------------------------------------------------------------




