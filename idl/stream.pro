function stream,q=q,P=P,M1=M1,tstream=tstream,n=n,DIST=dist,Tb=Tb,OMEGAB=omegab,K0=k0,$
                X0=x0,Y0=y0,Z0=z0,VX0=vx0,VY0=vy0,VZ0=vz0,rdisc=rdisc,discbrfac=discbrfac,$
                kindex=kindex,SUB=SUB
;,DRAG=drag
;+s=stream(
;                Q=mass ratio (secondary/primary)            
;                P=orbital period (days)                     
;                M1=primary mass (units of Msun)             
;                t=number of orbital periods to cover (0.25) 
;                n=number of points to calculate (t*12000)   
;                X0 starting position of stream relative to
;                Y0 centre of mass in units of a.
;                Z0 
;                dist=distance between points along to marked with
;                     circle in units of a
;
;Calculate path of accretion stream in a close binary. Calculation
;done in corotating frame.
  
print,'Calculating stream trajectory ......'

  if not keyword_set(n) then n=1000L
  if not keyword_set(discbrfac) then discbrfac=1d0
  if not keyword_set(rdisc) then rdisc=0d0
  if n_elements(kindex) eq 0 then begin
      if keyword_set(k0) then begin
          print,'Must set kindex.'
          return,0
      endif else kindex=0d0
  endif
  if not keyword_set(k0) then k0=0d0
  if not keyword_set(Tb) then begin
      if keyword_set(omegab) then Tb=2d0*!dpi/double(omegab) else Tb=0d0
  endif
  if not keyword_set(tstream) then tstream=0d0
  if not keyword_set(SUB) then SUB=500L

  x=dblarr(n)
  y=dblarr(n)
  z=dblarr(n)
  vix=dblarr(n)
  viy=dblarr(n)
  viz=dblarr(n)
  vkx=dblarr(n)
  vky=dblarr(n)
  vkz=dblarr(n)
  d=dblarr(n)
  lmag=dblarr(n)
  
  if keyword_set(X0) then x[0]=X0 else x[0]=(binary_xl1(q)-0.001d0)*binary_a(q=q,M1=M1,P=P)
  if keyword_set(Y0) then y[0]=Y0 else y[0]=0d0
  if keyword_set(Z0) then z[0]=Z0 else z[0]=0d0
  if keyword_set(VX0) then vix[0]=VX0 else vix[0]=0d0
  if keyword_set(VY0) then viy[0]=VY0 else viy[0]=0d0
  if keyword_set(VZ0) then viz[0]=VZ0 else viz[0]=0d0
  
  dummy=call_external(!LIBBINARY,"idl_stream",$
                      double(q),double(M1),double(P)*24d0*3600d0,$
                      double(Tb),double(k0), double(kindex), double(rdisc),double(discbrfac),$
                      double(tstream),long(n),long(SUB),$
                      x,y,z,vix,viy,viz,vkx,vky,vkz,d,lmag,/i_value)

  a=binary_a(q=q,M1=M1,P=P)
  if not keyword_set(dist) then dist=0.1d0 else dist=double(dist)
  nk=floor(max(d)/dist/a)+1     ;Number of points where vk is calculated 
  ik=intarr(nk)                 ;Indices of vi corresponding to vk points
  for i=0L,nk-1 do begin
    k=(where(d/dist/a ge double(i)))[0]
    ik[i]=k
  endfor
  
return,{x:x,y:y,z:z,vix:vix,viy:viy,viz:viz,vkx:vkx,vky:vky,vkz:vkz,ik:ik,d:d,lmag:lmag}

end

