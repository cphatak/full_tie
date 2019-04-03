pro colorwheel,dimx,dimy,bx,by,cimage,fname
;
; implement the color wheel representation for 
; magnetic induction maps
;
mmax = max(sqrt(bx^2+by^2))
bx=bx/mmax
by=by/mmax
cmag = sqrt((bx^2+by^2))*255.0
cimage = intarr(3,dimx+100,dimy)
red = bytarr(dimx,dimy)
green = bytarr(dimx,dimy)
blue = bytarr(dimx,dimy)
cang = 255.0*bx/cmag
sang = sqrt(1.0-cang^2)
; first green component
q=where((bx lt 0.0) and (by ge 0.0))
green(q)=byte(cmag(q)*abs(cang(q)))
q=where((bx ge 0.0) and (by lt 0.0))
green(q)=byte(cmag(q)*abs(sang(q)))
q=where((bx lt 0.0) and (by lt 0.0))
green(q)=byte(cmag(q))
; then red  
q=where((bx ge 0.0) and (by lt 0.0))
red(q)=byte(cmag(q))
q=where((bx ge 0.0) and (by ge 0.0))
red(q)=byte(cmag(q)*abs(cang(q)))
q=where((bx lt 0.0) and (by lt 0.0))
red(q)=byte(cmag(q)*abs(sang(q)))
; finally blue
q=where(by ge 0.0)
blue(q)=byte(cmag(q)*abs(sang(q)))
cimage(0,0:dimx-1,0:*)=red(0:*,0:*)
cimage(1,0:dimx-1,0:*)=green(0:*,0:*)
cimage(2,0:dimx-1,0:*)=blue(0:*,0:*)
;
; add color legend
;
rd= 50
red = bytarr(2*rd,2*rd)
green = bytarr(2*rd,2*rd)
blue = bytarr(2*rd,2*rd)
line=findgen(2*rd)-rd
bx=fltarr(2*rd,2*rd)
for i=0,2*rd-1 do bx(0:*,i)=line(0:*)
by=fltarr(2*rd,2*rd)
for i=0,2*rd-1 do by(i,0:*)=line(0:*)
mmax = max([abs(bx),abs(by)])
bx=bx/mmax
by=by/mmax
tr = shift(dist(2*rd),rd,rd)
qtr = where(tr gt rd)
cmag = sqrt((bx^2+by^2))*215.0
cang = 215.0*bx/cmag
sang = sqrt(1.0-cang^2)
; first green component
q=where((bx lt 0.0) and (by ge 0.0))
green(q)=byte(cmag(q)*abs(cang(q)))
q=where((bx ge 0.0) and (by lt 0.0))
green(q)=byte(cmag(q)*abs(sang(q)))
q=where((bx lt 0.0) and (by lt 0.0))
green(q)=byte(cmag(q))
; then red 
q=where((bx ge 0.0) and (by lt 0.0))
red(q)=byte(cmag(q))
q=where((bx ge 0.0) and (by ge 0.0))
red(q)=byte(cmag(q)*abs(cang(q)))
q=where((bx lt 0.0) and (by lt 0.0))
red(q)=byte(cmag(q)*abs(sang(q)))
; finally blue
q=where(by ge 0.0)
blue(q)=byte(cmag(q)*abs(sang(q)))
red(qtr)=0
green(qtr)=0
blue(qtr)=0
;
ypos = (dimy-2*rd)/2
cimage(0,dimx:dimx+2*rd-1,ypos:ypos+2*rd-1)=red(0:*,0:*)
cimage(1,dimx:dimx+2*rd-1,ypos:ypos+2*rd-1)=green(0:*,0:*)
cimage(2,dimx:dimx+2*rd-1,ypos:ypos+2*rd-1)=blue(0:*,0:*)
;
write_tiff,fname,reverse(cimage,3)
;
end

