; Example use of stream() in IDL.

P=0.1 ; Orbital period (days)
M1=0.7 ; Compact object mass (solar masses)
q=0.3 ; mass ration M2/M1
tstream=0.4 ; Time to run stream integration (orbits)

; Calculates orbital separation in metres
a=binary_a(q=q,M1=M1,P=P)

; Calculates stream data
str=stream(q=q,M1=M1,P=P,tstream=0.4,n=1000)

; The integration is a simple Eulerian integration, with 500*n
; integration steps. The returned arrays contain only every 500th
; point (ie. n points). n defaults to 1000 if not specified. 
; By default, the stream starts the integration 99.9% of the way from
; the centre of mass to the L1 point (to ensure the stream doesn't
; fall into the donor Roche lobe).
; Stream returns and IDL structure containing the results.
; The useful elements are:
; x,y,z - arrays containing the trajectory of the stream
;         relative to the centre of mass.
;         x is in the direction of M1 -> M2 in orbital plane
;         y is perpendicular to x also in the plane
;         z is in direction of orbital angular velocity vector
;         (all in metres)
; vix,viy,viz - arrays containing velocities along stream in
;               inertial frame, resolved along the x,y,z directions
;               (all in metres/second)
; d - array containing integrated distance along the stream path (in metres)
; Elements can be address as [stream_structure].[name],
; e.g. str.x, str.vix in this example.

; Plot stream trajectory
plot,str.x/a,str.y/a,xrange=[-1.,1.],yrange=[-1.,1.]

; Pause for 5 seconds
wait,5

; Plot stream velocities
plot,str.vix,str.viy,xrange=[-2d6,2d6],yrange=[-2d6,2d6]

end 
