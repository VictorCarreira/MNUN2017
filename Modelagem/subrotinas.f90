real function fon_d2g(t,fc,fa)
!*** function fonte sismica ***
!     Fonte: derivada segunda da Gaussiana

implicit none

real*4    t,fc,pi,fa
parameter (pi=3.141593)

fon_d2g = -pi*(pi*fc*t)*(pi*fc*t)
fon_d2g = exp(fon_d2g)
fon_d2g = fa*fon_d2g*(1.-2.*pi*(pi*fc*t)*(pi*fc*t))

end
