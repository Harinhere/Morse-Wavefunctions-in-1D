c Main program  
      implicit double precision(a-h,o-z)
      open(13,file='mors.out')
c sqm is sqrt(2.*am),am=reduced mass	  
      sqm=sqrt(85462.9378)	  
      am=sqm*sqm/2.
c V=an*(exp(-alf*r)-1)**2
      alf=.8507
      an=.1382
c nv1 is the vibrational level
      nv1=1	 
      rstep=0.005
      rmin=-2.
      rmax=2.
      r1e=0.
      irs=(rmax-rmin)/rstep
      do 1 ir=0,irs
	    r=rmin+ir*rstep
	    rr=r-r1e
            call wavemors(an,alf,am,nv1,rr,ev1,fm1)		
            write(13,*) r,fm1
 1    continue
      end
c subroutine wavemors returns the eigenfunctions and eigenvalues	  
      subroutine wavemors(an,alf,ams,nv,x,ev,fm0)	  
      implicit double precision(a-h,o-z)    
      sqm=sqrt(2.*ams)
c eigen energy
      r1=.5+nv
c Lambda
     ak=2.*sqm*sqrt(an)/alf	  
c omega
     omg=2.*alf*sqrt(an)/sqm
     ev=r1*omg*(1.-r1/ak)
c initial eigenfunction ps0 for recurrence relation
     pi=2.*asin(1.d0)
     s2pi=sqrt(2.*pi)
     s=ak-2.*nv-1.
     y=ak*exp(-alf*x)
     z=ak-nv
     r2=s*log(y)/2.-(z-.5)*log(z)/2.
     r2=r2-y/2.+z/2.
     pf=exp(r2)*sqrt(alf*s/s2pi)
c gamma(2*lambda)=gamma(s), use the asymptotic series expansion
c  8 terms in the expansion
          t1=1.
          t2=(1./12.)/z	  
	  t3=(1./288.)/z**2
	  t4=(139./51840.)/z**3
	  t5=(571./2488320.)/z**4
	  t6=(163879./209018880.)/z**5
	  t7=(5246819./75246796800.)/z**6
	  t8=(534703531./902961561600.)/z**7
c psi0
      pf=pf/sqrt(t1+t2+t3-t4-t5+t6+t7-t8)       	  
c now do the recursion for factorial
      do 4 n=1,nv
	  gn=n
 4    pf=pf*sqrt(gn)
c recursion relation for Laguerre polynomials
      fm_1=0.
      fm0=1.	  
      do 2 n=1,nv
	  cf1=2.*n+s-1-y
	  cf2=s+n-1
	  fm1=(cf1*fm0-cf2*fm_1)/n
	  fm_1=fm0
	  fm0=fm1
 2    continue
      fm0=pf*fm0
 20   return  
      end
