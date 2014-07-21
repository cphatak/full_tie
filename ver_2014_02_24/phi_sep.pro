;+
;program to separate the magnetic phase and electrostatic phase shift
;and store the data. Generates tiff files for phase, bxt, byt, and colorwheel.
;@author, Charudatta Phatak, ANL, 11/16/09.
;-

pro phi_sep,unflip,flip
compile_opt idl2

;load the phase files
;unflipped
openr,1,unflip.path+unflip.prefix+'.phase'
q=assoc(1,intarr(3))
dim=q[0]
num=dim[2]
dim=dim[0]
d2=dim/2
q=assoc(1,fltarr(dim,dim),12)
phiarru = fltarr(dim,dim,num)
for ii=0,num-1 do phiarru[*,*,ii]=q[ii]
close,1

;flipped phase
openr,1,flip.path+flip.prefix+'.phase'
q=assoc(1,fltarr(dim,dim),12)
phiarrf = fltarr(dim,dim,num)
for ii=0,num-1 do phiarrf[*,*,ii]=q[ii]
close,1

;for rotating the array
ang=54
openr,1,flip.path+flip.prefix+'_final.aligned'
;ndim=1024
;nd2=ndim/2
q=assoc(1,fltarr(dim,dim),12)
mm = q[flip.num]
close,1
;erode mask to eliminate errors on edge
mm=erode(mm,replicate(1,5,5))
mm2 = rot(mm,ang,cubic=-0.5,missing=0.0)
;mm = mm[nd2-d2:nd2+d2-1,nd2-d2:nd2+d2-1]

;gradient kernels
kx=([[0.0,0.0,0.0],[-1.0,0.0,1.0],[0.0,0.0,0.0]])
ky=([[0.0,1.0,0.0],[0.0,0.0,0.0],[0.0,-1.0,0.0]])

mphi = 0.5*(phiarru-phiarrf);*(-1)
ephi = 0.5*(phiarru+phiarrf);*(-1)

;mask=mkmask(total(mphi,3),2.0)

dimx=dim & dimy=dim
for ii=0,num-1 do begin
	;mphi[*,*,ii] *= mask
	fn = flip.path+'/images/'+flip.prefix+'_phisep'+string(ii,format='(I01)')
	write_tiff,fn+'_mphase.tiff',(reverse(reform(mphi[*,*,ii]),2)),/float
	write_tiff,fn+'_ephase.tiff',(reverse(reform(ephi[*,*,ii]),2)),/float
	bxt = convol(reform(mphi[*,*,ii]),ky,/center,/edge_truncate)*mm
	byt = convol(reform(mphi[*,*,ii]),kx,/center,/edge_truncate)*mm
	write_tiff,fn+'_bxt.tiff',(reverse(bxt,2)),/float
	write_tiff,fn+'_byt.tiff',(reverse(byt,2)),/float
	colorwheel,dimx,dimy,bxt,byt,cim,fn+'_color.tiff'
  ;mmphi=rot(reform(mphi[*,*,ii]),ang,cubic=-0.5,missing=0.0)
  ;bxt = convol(mmphi,ky,/center,/edge_truncate)*mm2
  ;byt = convol(mmphi,kx,/center,/edge_truncate)*mm2
  ;colorwheel,dimx,dimy,bxt,byt,cim,fn+'_rot_color22.tiff'
endfor

;save the magnetic phase shift data..
openw,1,flip.path+flip.prefix+'_magphi.data'
q=assoc(1,intarr(3))
q[0]=[dim,dim,num]
q=assoc(1,fltarr(dim,dim),12)
for ii=0,num-1 do q[ii] = mphi[*,*,ii]
close,1

openw,1,flip.path+flip.prefix+'_ephi.data'
q=assoc(1,intarr(3))
q[0]=[dim,dim,num]
q=assoc(1,fltarr(dim,dim),12)
for ii=0,num-1 do q[ii] = ephi[*,*,ii]
q[num]=mm
close,1

stop
;computing the net magnetization along the x-y axes of array
print,'Computing net magnetization'
openr,1,unflip.path+unflip.prefix+'.aligned'
q=assoc(1,fltarr(dim,dim),12)
mask = q[0]
close,1
openr,1,flip.path+flip.prefix+'_final.aligned'
q=assoc(1,fltarr(dim,dim),12)
mm = q[flip.num]
close,1
;erode mask to eliminate errors on edge
mm=erode(mm,replicate(1,5,5))
mask *= mm

;get bxt & byt
bxt = convol(reform(mphi[*,*,num-1]),ky,/center,/edge_truncate)
byt = convol(reform(mphi[*,*,num-1]),kx,/center,/edge_truncate)

window,0,xsi=dim,ysi=dim
mask = smooth(bytscl(mask) lt 170,3)*mm
tvscl,mask
wait,0.4
lmask = label_region(mask,/ulong)
lmax = max(lmask)
print,lmax
;next get the bxt and byt images and for each island
;determine (bxt,byt) and then assign it +1 or -1 (mx,my)
;and compute the net magnetization.

;mx & my arrays measured along the length of the lattice
mx = 0.0
my = 0.0
;colormask similar to simulations
colormask = bytarr(3,dim,dim)

;loop over all the islands
for ii=1,lmax do begin
  if (ii mod 10 eq 0) then print,ii else print,format='(A1,$)','.'
  ;get the island
  isle = lmask eq ii
  ;get bxt and byt for the island
  bxt_is = total(bxt*isle)
  byt_is = total(byt*isle)
  ;4 cases
  if ((bxt_is gt 0) AND (byt_is gt 0)) then begin
    my += 1.0
    colormask[1,*,*] += isle*255
    colormask[2,*,*] += isle*255
  endif
  if ((bxt_is gt 0) AND (byt_is lt 0)) then begin
    mx += 1.0
    colormask[0,*,*] += isle*255
  endif
  if ((bxt_is lt 0) AND (byt_is gt 0)) then begin
    mx -= 1.0
    colormask[1,*,*] += isle*128
  endif
  if ((bxt_is lt 0) AND (byt_is lt 0)) then begin
    my -= 1.0
    colormask[0,*,*] += isle*255
    colormask[1,*,*] += isle*255
  endif
endfor
print,'Mx, My, |M|',mx,my,sqrt(mx^2+my^2)
tvscl,colormask,true=1
xyouts,0.05,0.95,'Mx='+string(mx,format='(F7.2)')+',My='+string(my,format='(F7.2)')+$
  ',|M|='+string(sqrt(mx^2+my^2),format='(F6.2)'),charsize=3,/normal
cc=tvrd(/order,true=1)
write_tiff,fn+'_colormask_netmu.tiff',cc

stop
end

