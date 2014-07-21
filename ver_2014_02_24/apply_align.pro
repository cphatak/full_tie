pro apply_align,fd,cleanup=cleanup
;
; this routine applies the alignment stored in ###.align to the 
; data stored in ###.decon, and produces a file ###.aligned; 
; optionally, it will also remove the ###.decon file, since this
; can easily be recreated using MTF_decon.pro.
;
num = fd.num

param = fltarr(num,4)
param[0,0:3] = [0.0,0.0,0.0,1.0]

a=0.0 & b=0.0 & c=0.0 & d=0.0
openr,1,fd.path+fd.prefix+'.align'
for i=1,num-1 do begin
	readf,1,a,b,c,d
	param[i,0:3] = [a,b,c,d]
endfor
close,1

; next, apply these scaling and translation parameters to each of the 
; images in prefix.decon, and create a new (huge) file prefix.aligned,
; which will, in turn be read by the TIE program...

; first, get the array size and create the interpolation arrays
print,format="('determining image size ',$)"
openr,10,fd.path+fd.prefix+'.decon'
q = assoc(10,intarr(3))
d = q[0]
dim = d[0]
close,10
print,format="(I5,'x',I5,' pixels')",dim,dim

window,0,xsi=dim/4,ysi=dim/4,title='aligned images',retain=2
window,1,xsi=dim/4,ysi=dim/4,title='common mask',retain=2

print,'creating interpolation arrays'
line = findgen(dim)
bintx = replicate(1.0,dim)##line
binty = replicate(1.0,dim)#line 

print,'starting alignments'
openr,10,fd.path+fd.prefix+'.decon'
q=assoc(10,fltarr(dim,dim),12)

openw,15,fd.path+fd.prefix+'.aligned'
qq = assoc(15,intarr(3))
qq[0] = [dim,dim,num]
qq=assoc(15,fltarr(dim,dim),12)

print,'storing in-focus image (no alignment needed)'
im = q[0]
qq[0] = im
infocus = im
print,'allocating common mask'
mask = replicate(1.0,dim,dim)
for i=1,num-1 do begin
	print,format="('starting image ',I5,' ',$)",i+1
  print, param[i,*]
	im = q[i]
	b = rot(im,param[i,2],param[i,3],missing=-1.0,cubic=-0.5,/pivot)
	;b = rot(im,0.0,1.0,missing=-1.0,cubic=-0.5,/pivot)
	im = interpolate(b,bintx+param[i,0],binty + param[i,1],cubic=-0.5,missing=-1.0)
	s = where(im eq -1.0,count)
	if (count gt 0) then mask[s] = 0.0
	wset,0
	tvscl,rebin(im,dim/4,dim/4)
	wset,1
	tvscl,rebin(mask,dim/4,dim/4)
	qq[i] = im
	print,' -> done'
endfor
close,10

print,'saving mask to .aligned file as last image'
qq[num] = mask
close,15

print,'Total number of pixels in image = ',long(dim)^2
t = total(mask)
print,'Total number of pixels in mask = ',long(t)
print,'Percentage of pixels available for reconstruction = ',t/float(dim)^2*100.0

wset,1
tvscl,rebin(infocus * mask,dim/4,dim/4)

if keyword_set(cleanup) then begin
	ans = ''
	print,'Keyword >cleanup< was detected.'
	read,'Are you sure that you want to delete the file '+fd.prefix+'.decon  ?  (y/n) ',ans
	if (ans eq 'yes') then begin
		spawn,'/bin/rm '+fd.path+fd.prefix+'.decon'
		print,fd.path+fd.prefix+'.decon deleted'
		print,' (This file can always be reconstructed using MTF_decon.pro) '
	endif
endif

wdelete,0
wdelete,1

end
