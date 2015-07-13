function roche_potl,p,q
; Calculate Roche potl at position p[x,y,z] in plane of disc in J/kg / (G*M1*a)
; Origin=centre of mass
; x = line of centres -> secondary
; y = perp to x, in orbital plane
; z = normal to orb plane in direction of 'omega vector'
; q = mass ratio M2/M1
; All distance in units of a (binary separation)

potl=call_external(!LIBBINARY,"idl_roche_potl",$
                double(p[0]),double(p[1]),double(p[2]),double(q),/d_value) 

return,potl

end
