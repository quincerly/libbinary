function binary_a,M1=m1,q=q,P=P
;+ a=binary_a(M1=m1,q=q,P=P)
;
;	M1 in solar masses
;	q
;	P in days
;	a in metres
;-

a=call_external(!LIBBINARY,"idl_binary_sep",$
                double(q),double(M1),double(P)*24d0*3600d0,/d_value)

return,a
                     
end
