function mkmask,image,sig,blur=blur
device,cursor_standard=43

dim=size(image,/dim)
dim=dim[0]
d2=dim/2
wsi=dim/4
print,dim
if (~keyword_Set(blur)) then begin
window,8,xsi=wsi,ysi=wsi,retain=2
tvscl,rebin(image,wsi,wsi)
print,'Select LL and TR corners of mask'
cursor,x1,y1,4,/dev
plots,[x1,x1],[0,dim-1],/dev,color=200
plots,[0,dim-1],[y1,y1],/dev,color=200
cursor,x2,y2,4,/dev
plots,[x2,x2],[0,dim-1],/dev,color=200
plots,[0,dim-1],[y2,y2],/dev,color=200
mask=fltarr(dim,dim)
scl=dim/wsi
mask[x1*scl:x2*scl,y1*scl:y2*scl]=1.0
wdelete,8
endif else mask=image
if (n_elements(sig) eq 0) then wexp=200.0 else wexp=sig
x=findgen(dim)-float(d2)
y=x
!EXCEPT=0
s=exp(-x^2/wexp)#exp(-y^2/wexp)
norm=total(s)
fm=fft(mask,-1)
fs=fft(s,-1)
mask=float(shift(float(fft(fm*fs,1)),d2,d2)*dim*dim/norm)
fm=0 & fs=0 &x=0 & y=0 & s=0
!EXCEPT=1
return,mask
end
