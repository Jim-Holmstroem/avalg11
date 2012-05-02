from fractions import gcd

def f(x,c,n):
    return (x**2+c)%n

def pollardrho(n,c):
    x=2
    y=2
    d=1
    while(d==1):
        x=f(x,c,n)
        y=f(f(y,c,n),c,n)
        d=gcd(abs(x-y),n)
    return d


def test_pollardrho(n):
    c=1
    factor=n
    while(factor==n):
        while(c%n==0 or c%n==-2):
            c=c+1
        print "n=",n,",c=",c,
        factor = pollardrho(n,c)
        c=c+1
        
    return factor
