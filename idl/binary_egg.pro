function binary_egg,q
;+ RL1=binary_RL1(q)
;
;	RL1 in units of binary separation a
;       relative to c.o.m. -> secondary
;-

if n_elements(q) gt 1 then begin
    egg=q*0d0
    for i=0,n_elements(q)-1 do egg[i]=call_external(!LIBBINARY,"idl_eggleton_radius",$
                                                 double(q[i]),/d_value)
endif else egg=call_external(!LIBBINARY,"idl_eggleton_radius",$
                             double(q),/d_value)

return,egg

end
