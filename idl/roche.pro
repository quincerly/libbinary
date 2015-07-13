function unit_vector,theta,phi
; Return unit vector in direction (theta,phi)
;	theta = angle between direction and x-axis
;	phi   = azimuthal angle about x-axis measured from z to y ..

x=cos(theta)
y=sin(theta)*sin(phi)
z=sin(theta)*cos(phi)

return,[x,y,z]

end


function roche,q=q,n_angles=n_angles,r_steps=r_steps,DISPLAY=display,XZ=xz

; All distance in units of a (binary separation)
; Distance of L1 point from centre of secondary

if not keyword_set(n_angle) then n_angles=60
if not keyword_set(r_steps) then r_steps=100

;print,'Calculating Roche geometry ......'

xl1=binary_xl1(q)		

pl1=[xl1,0d0,0d0]		;Position of L1 point

pM1=[-q/(1d0+q),0d0,0d0]	;Position of M1
pM2=[1/(1d0+q),0d0,0d0]		;Position of M2

potl_l1=roche_potl(pl1,q)	;Roche potential at L1 point

phi=!dpi/2d0
if keyword_set(xz) then phi=0d0

;Calculate primary Roche lobe
x1=dblarr(n_angles)
y1=dblarr(n_angles)
z1=dblarr(n_angles)
s1=(dindgen(r_steps)+1)/r_steps*(q/(1d0+q)+xl1)
p1=dblarr(r_steps)
for i=0,n_angles-1 do begin
	dirn1=unit_vector(2.0*!pi*double(i)/double(n_angles),phi)
	for j=0,r_steps-1 do p1[j]=roche_potl(pM1+s1[j]*dirn1,q) 
	r=interpol(s1,p1,potl_l1)
	x1[i]=pM1[0]+r*dirn1[0]
	y1[i]=pM1[1]+r*dirn1[1]
	z1[i]=pM1[2]+r*dirn1[2]
endfor

;Calculate secondary Roche lobe
x2=dblarr(n_angles)
y2=dblarr(n_angles)
z2=dblarr(n_angles)
s2=(dindgen(r_steps)+1)/r_steps*(1d0/(1d0+q)-xl1)
p2=dblarr(r_steps)
for i=0,n_angles-1 do begin
	dirn2=unit_vector(2.0*!pi*double(i)/double(n_angles),phi)
	for j=0,r_steps-1 do p2[j]=roche_potl(pM2+s2[j]*dirn2,q) 
	r=interpol(s2,p2,potl_l1)
	x2[i]=pM2[0]+r*dirn2[0]
	y2[i]=pM2[1]+r*dirn2[1]
	z2[i]=pM2[2]+r*dirn2[2]
endfor

;if keyword_set(display) then begin
;	np=50
;	rp=dblarr(np,np)
;	mask=dblarr(np,np)
;	xx=dindgen(np)/double(np-1)*1.2d0-(q/(1d0+q)+0.1)
;	yy=(dindgen(np)/double(np-1)-0.5d0)*1.2d0
;	for i=0,np-1 do begin
;		for j=0,np-1 do begin
;			rp[i,j]=roche_potl([xx[i],yy[j],0d0],q)
;		endfor
;	endfor
;	mask[*,*]=1d0
;	mask[where(rp lt potl_l1)]=0d0

;	bitmap,rp,xrange=[min(xx),max(xx)],yrange=[-0.6,0.6],pow=0.5d0,mask=mask
;	oplot,[x1,x1[0]],[y1,y1[0]]
;	oplot,[x2,x2[0]],[y2,y2[0]]
;endif

return,{x1:[x1,x1[0]],y1:[y1,y1[0]],z1:[z1,z1[0]],x2:[x2,x2[0]],y2:[y2,y2[0]],z2:[z2,z2[0]]}

end
