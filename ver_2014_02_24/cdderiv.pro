;+
;program to compute the defocus derivative of a set of
;Nf through focus images by estimating a linear comb.
;of M - f(z)'s for intensity. The linear exp. coeff. 
;are computed using SVD.
;@author, Charudatta Phatak, cmu, 05.03.07
;-
function cdderiv,imstack,defocus,mask=mask,verbose=verbose
compile_opt idl2

;parameters setup..
ndims=size(imstack,/dim)
dim=ndims[0]
Nf=ndims[2] ;No. of Defocus images
M=3 ;No. of terms in the expansion to be used (a0+a1*z+a2*z^2+....)
def=defocus;
if (keyword_Set(verbose)) then print, def
;def=(findgen(Nf)-float(Nf/2))*defocus
a0=fltarr(dim,dim)&a1=fltarr(dim,dim)&a2=a1
q=array_indices(mask,where(mask eq 1.0,count))
print,'Total no. of pixels for dIdz cal.',count
sing_val_tot=0
sing_val=0
for c=0,count-1 do begin
	;print,'.',format='(A1,$)'
	if keyword_set(verbose) then if (c mod 1000 eq 0) then if (c mod 10000 eq 0) then print,c else print,'.',format='(A1,$)'
	;for y=0,dim-1 do begin
	x=q[0,c]&y=q[1,c]
	;xmin=x-16 > 0& xmax=x+16 < dim-1
	;ymin=y-16 > 0& ymax=y+16 < dim-1
	;A=fltarr(M,Nf)&b=fltarr(Nf)
	;for nn=0,Nf-1 do begin
	;	sig=stddev(imstack[xmin:xmax,ymin:ymax,nn])
	;	if (sig eq 0.0) then sig=1.0
	;	A[0,nn]=1.0/sig&A[1,nn]=def[nn]/sig&A[2,nn]=def[nn]^2/sig
	;	b[nn]=imstack[x,y,nn]/sig
	;endfor
	;svdc,A,w,u,v
	;zinds=where(w eq 0.0)
	;if (zinds[0]  ne -1) then print,'W was zero for',x,y
	;res=reform(((u#b)/w)#v)

	;use of SVDFIT function to get the values
	xsvdfit=findgen(Nf)-float(Nf/2)
	xsvdfit=def
	ysvdfit=fltarr(Nf)
	for nn=0,Nf-1 do ysvdfit[nn]=imstack[x,y,nn]
	res=svdfit(xsvdfit,ysvdfit,M,/double,status=sing_val)
	a0[x,y]=res[0]&a1[x,y]=res[1]&a2[x,y]=res[2]
sing_val_tot += sing_val
endfor
;endfor
print,'Singular values found:',sing_val_tot
return,a1

end
