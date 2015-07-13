function binary_xl1,q
;+ RL1=binary_RL1(q)
;
;	position of L1 in units of binary separation a
;       along x axis relative to c.o.m. -> secondary
;-

n=n_elements(q)

if n ne 1 then begin
    xl1=dblarr(n)
    for i=0L,n-1 do begin
        xl1t=call_external(!LIBBINARY,"idl_roche_xl1",$
                          double(q[i]),/d_value)
        xl1[i]=xl1t
    endfor
endif else xl1=call_external(!LIBBINARY,"idl_roche_xl1",$
                double(q),/d_value)

return,xl1

end

