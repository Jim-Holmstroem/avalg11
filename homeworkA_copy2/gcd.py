def gcd(a,b):
  if b > a:
    print "sorting"
    return gcd(b, a)
  elif b == 0:
    print "b=0"
    return a
  elif a%2==0 and b%2==0: #a and b even
    print "a and b even"
    return 2*gcd(a/2, b/2)
  elif a%2==0: #a even
    print "a even"
    return gcd(a/2, b)
  elif b%2==0: #b even
    print "b even"
    return gcd(a, b/2)
  else:
    print "gcd(b, a-b) /w (a,b)=(",a,",",b,")"
    return gcd(b, a-b)
