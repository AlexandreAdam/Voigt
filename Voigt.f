      function voigt(y,x)
      implicit real*8 (a-h,o-z)

c     generate hjerting functions
c     copied straight out of j.humlicek, 1979, j.quant. spec.
c          and rad. transfer, 21,309
c
c     voigt is voigt(y,x) -- the real part of w(z)=exp(-z**2)*erfc(-iz) in
c     the upper 1/2 plane y>=0
c
c	voigt(a,v) = hjerting fn where a=gamma/4pi dnu(dop)
c        dnu(dop) = sqrt(2kt/m)/lambda0
c         v = (nu-nu0)/dnu(dop)   no. of doppler widths from line center
c
c     checked by comparison with finn and mugglestone, 1965, mnras,
c       129,221
c
c     error < o(10**-5)  larger than advertised, but ok
c

c     copied form jill bechtold (steward observatory)

      dimension t(6),c(6),s(6)

      data(t(i),i=1,6)/.314240376,.947788391,1.59768264,2.27950708,
     *3.02063703,3.8897249/
      data(c(i),i=1,6)/1.01172805,-.75197147,1.2557727e-2,
     *1.00220082e-2,-2.42068135e-4,5.00848061e-7/
      data(s(i),i=1,6)/1.393237,
     *.231152406,-.155351466,6.21836624e-3,9.19082986e-5,-6.27525958e-7/

      voigt=0.
      y1=y+1.5
      y2=y1*y1

      if(y.gt.0.85.or.abs(x).lt.18.1*y+1.65) go to 2

c     region ii

      if(abs(x).lt.12.) voigt=exp(-x*x)
      y3=y+3.

      do 1 i=1,6
 	 r=x-t(i)
 	 r2=r*r
 	 d=1./(r2+y2)
 	 d1 = y1*d
	 d2 = r * d
	 voigt=voigt+y*(c(i)*(r*d2-1.5*d1)+s(i)*y3*d2)/ (r2+2.25)
	 r=x+t(i)
	 r2=r*r
	 d=1./(r2+y2)
	 d3=y1*d
	 d4=r*d
	 voigt=voigt+y*(c(i)*(r*d4-1.5*d3)-s(i)*y3*d4)/(r2+2.25)
    1 continue

      return

c     region i

    2 do 3 i=1,6
	 r=x-t(i)
	 d=1./(r*r+y2)
	 d1=y1*d
	 d2=r*d
	 r=x+t(i)
	 d=1./(r*r+y2)
         d3=y1*d
	 d4=r*d
	 voigt=voigt+c(i)*(d1+d3)-s(i)*(d2-d4)
    3 continue

      return
      end

