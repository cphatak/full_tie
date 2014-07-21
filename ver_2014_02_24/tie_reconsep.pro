@mdgsolvetie
@simulationparameters
@cdderiv
;+
;program to perform TIE reconstructions on series of 
;through focus images at different foci. This uses the
;mdgsolvetie to compute the phase and store the phase for
;each series as a data file. The image and the surfaceplot
;are also saved alongwith.
;@author, Charudatta Phatak, CMU, 12/04/08.
;-
;modified for direct reconstruction of mphase & ephase

pro TIE_reconsep,fd,flipfile=flipfile,dcdfile=dcdfile

dim=fd.dim ;2048
d2=dim/2
scl=1.0
ndim=dim/scl
nd2=ndim/2
close,/all

;defocus values
;defstep = 10000.0 ;series 1
;defvals=fd.defstep*(findgen(fd.num/2)+1.0)
;defstep=fd.defstep
;defvals = [1105920,2211840]*0.18
defvals = [3133440]*0.18
;defvals = [23040,115200,322560,506880]*0.18
;defvals = [552960,1105920,2027520,3133440,5160960]*0.18

;No. of images
num=fd.num

;get the scope paramters and TIE pointers
;dd=4.13 ;pixel size in nm
dd=fd.delta*1000.0
delta=dd
EE=200.0D3 ;V
cs=0.2D6 ;nm
cc=4.0D6 ;nm
thetac=0.0001D0 ;rad
sigV=1.0D-6 ;sig(V)/V
sigI=1.0D-6 ;sig(I)/I
sigE=2.0D-6 ;sig(E)/E
spread= cc * sqrt(sigV^2 + sigI^2 + sigE^2)

wsi=1024
window,0,xsi=wsi,ysi=wsi,retain=2;,xpos=2000
;array to store all the phase

;rotational angle
rot_ang = 0.0

;load all the tfseries images
;unflipped
openr,1,fd.path+fd.prefix+'.aligned'
q=assoc(1,fltarr(dim,dim),12)
unstack=fltarr(dim,dim,num)
for ii=0,num-1 do unstack[*,*,ii]=median(rot(q[ii],rot_ang,cubic=-0.5,missing=0.0),3)
unflip_mask = rot(q[num],rot_ang,missing=0.0)
close,1

;flipped
openr,1,flipfile.path+flipfile.prefix+'_final.aligned'
q=assoc(1,fltarr(dim,dim),12)
flstack = fltarr(dim,dim,num)
for ii=0,num-1 do flstack[*,*,ii]=median(rot(q[ii],rot_ang,cubic=-0.5,missing=0.0),3)
flip_mask = rot(q[num],rot_ang,missing=0.0)
close,1

;common mask
mask = unflip_mask*flip_mask
mask=erode(mask,replicate(1.0,9,9))
;select the central contiguous area
ll=label_Region(mask)
mask=ll eq ll[d2,d2]

;stop

;option to select a smaller 1K by 1K area from the images - 
ans=''
cen=0
;change to smaller dimensions..
ndim=512
nd2=ndim/2
;xs=0 & xe=dim-1
;ys=0 & ye=dim-1
read,prompt='Do you want to select a smaller '+strtrim(string(ndim),2)+' by '+strtrim(string(ndim),2)+' region? (y/n): ',ans
if (ans eq 'y') then begin
 read,prompt='Do you want to select a central region (0) or custom region (1) : ',cen
  if (cen eq 1) then begin
    print,'Select the center of the area:'
    sat=0
    while (sat eq 0) do begin
		;erase
    tvscl,rebin(unstack[*,*,0]+flstack[*,*,0],wsi,wsi)
    cursor,xx,yy,4,/dev
    plots,xx,yy,/dev,psym=7,symsize=2,color=200
    ;xx=936.0/2 & yy=1218.00/2.0
    ;xx = 1018/2. & yy=1088/2.0
    ;xx=1380./2 & yy=1516./2
    xyouts,xx+10,yy+10,string(xx*2)+','+string(yy*2),/dev,color=200,charsize=2.0
    xx*=2.0 & yy*=2.0
    if (xx lt dim/2) then begin
			xs = 0 > (xx-ndim/2)
			xe = (xx+ndim/2-1) > (ndim-1)
		endif else begin
			xs = (dim-1-ndim+1) < (xx-ndim/2)
			xe = (xx+ndim/2-1) < (dim-1)
		endelse
		if (yy lt dim/2) then begin
			ys = 0 > (yy-ndim/2)
			ye = (yy+ndim/2-1) > (ndim-1)
		endif else begin
			ys = (dim-1-ndim+1) < (yy-ndim/2)
			ye = (yy+ndim/2-1) < (dim-1)
		endelse
     plots,[xs/2,xe/2,xe/2,xs/2,xs/2],[ys/2,ys/2,ye/2,ye/2,ys/2],/dev,color=200
     ;plots,[xs,xe,xe,xs,xs],[ys,ys,ye,ye,ys],/dev,color=200
     read,prompt='Are you satisfied with the selection (0/1):',sat
    cc=tvrd(/order,true=1)
    write_tiff,flipfile.path+'/images/'+flipfile.prefix+'_selarea.tiff',cc
   endwhile
 ;xs*=2 & xe*=2 & ys*=2 &ye*=2
 endif else begin
	xs = d2-nd2 & xe = d2+nd2-1
	ys = d2-nd2 & ye = d2+nd2-1
endelse
unstack=unstack[xs:xe,ys:ye,0:*]
flstack=flstack[xs:xe,ys:ye,0:*]
mask=mask[xs:xe,ys:ye,0:*]
dim=ndim
d2=dim/2
endif

;check if user wants to use CDDERIV or not.
chk_deriv=''
read,prompt='Do you want to use polynomial fit for calculating intensity derivative? (y/n):',chk_deriv

if (chk_deriv eq 'n') then begin
;ndim=512
phiarr_m=fltarr(dim,dim,num/2)
phiarr_e=fltarr(dim,dim,num/2)
;Loop to start reconstruction for each defocus value
for k=1,num/2 do begin
        ;Linear Focal series
	;defval = k*defstep
	;Quadratic Focal series
	;defval = k^2*defstep
	;Preset Defocus series
	defval = defvals[k-1]
	
	print,'Reconstructing phase for defocus =',defval

	;load the images
	imstack=fltarr(dim,dim,5)
	;unflip, underfocus +-
  	imstack[*,*,0] = unstack[*,*,k]*mask
	;infocus
  	imstack[*,*,2] = (unstack[*,*,0]+flstack[*,*,0])/2.0*mask
  	;imstack[*,*,2] = (unstack[*,*,0]);+flstack[*,*,0])/2.0*mask
	; unflip, overfocus ++
  	imstack[*,*,3]=unstack[*,*,num/2+k]*mask
	;flip, underfocus --
  	imstack[*,*,1] = flstack[*,*,k]*mask
	;flip, overfocus -+
  	imstack[*,*,4]=flstack[*,*,num/2+k]*mask
	

;initialize the pointers for getting various parameters
	params=simulationparameters(dim=dim,delta=dd,EE=ee,cs=cs,$
		thetac=thetac,spread=spread,defstep=defval)
	pimage=params[0]
	pscope_ltem=params[1]
	pphyscon=params[2]
	ptie=params[3]

  ;### make sure the integrated intensity in each image is the same ###
	t = total(total(imstack,1),1)
	print,'image totals : ',t
	tm = max(t)
	t = tm/t
	for kk=0,4 do imstack[0:*,0:*,kk] *= t[kk]
	t = total(total(imstack,1),1)
	print,'image totals : ',t
	;t = total(total(fimstack,1),1)
	;print,'image totals : ',t
	;tm = max(t)
	;t = tm/t
	;for kk=0,2 do fimstack[0:*,0:*,kk] *= t[kk]
	;t = total(total(fimstack,1),1)
	;print,'image totals : ',t


	;####scale image intensity to interval [0,1] and make
	; sure that in-focus image is shifted away from zero intensity####
  ;stop
	imax = max(imstack)
	imstack /= imax
	;for kk=0,2 do begin
  ;  temp = reform(imstack[*,*,kk])
  ;  temp /= max(temp)
  ;  temp = (temp-mean(temp))/stddev(temp)
  ;  imstack[*,*,kk] = temp
  ;endfor
  imstack += 0.001D0
	imstack[0:*,0:*,2] += 1.D0-mask[*,*]
  ;imstack[*,*,1] *= hotpix
  ;fimax = max(fimstack)
  ;fimstack /= fimax
  ;fimstack += 0.001D
  ;fimstack[*,*,1] += 1.D0-mask[*,*]

	;## compute intensity derivatives
	dIdz_m = (reform(imstack[*,*,3]) - reform(imstack[*,*,0]) $
    - reform(imstack[*,*,4]) + reform(imstack[*,*,1]))/2.0
	dIdz_e = (reform(imstack[*,*,3])-reform(imstack[*,*,0]) $
    + reform(imstack[*,*,4]) - reform(imstack[*,*,1]))/2.0
	;## next, make sure that the derivative has zero total "energy" 
	tot = total(dIdz_m)/total(mask)
	dIdz_m -= tot
	tot = total(dIdz_e)/total(mask)
	dIdz_e -= tot
	
	; in-focus image 
	(*ptie).ima_if=reform(imstack[0:*,0:*,2],dim,dim)

	tvscl,rebin(dIdz_m,wsi/2,wsi/2)
	tvscl,rebin(dIdz_e,wsi/2,wsi/2),wsi/2,0
	tvscl,rebin(imstack[*,*,0],wsi/4,wsi/4),0
	tvscl,rebin(imstack[*,*,1],wsi/4,wsi/4),1
	tvscl,rebin(imstack[*,*,2],wsi/4,wsi/4),2

	;solve for each phase
	;e-phase shift
	(*ptie).ima_dIdz = dIdz_e
	;now perform phase reconstruction
  	;mdgsolvetie,'laplacian',ptie,/symmetrize;,tikhonov=2.0
  	mdgsolvetie,'laplacian',ptie,tikhonov=2.0
	;mdgsolvetie,'inverse gradient',ptie;,/symmetrize 
	
	;display the phase
	phi_e=(*ptie).phase
	et_x=(*ptie).bxt
	et_y=(*ptie).byt

	tvscl,rebin(phi_e,wsi/2,wsi/2),0
	plot,phi_e[*,nd2],/noerase
	;phiarr[*,*,k-1]=phi

	;m-phase shift
	(*ptie).ima_dIdz = dIdz_m
	;now perform phase reconstruction
  	;mdgsolvetie,'laplacian',ptie,/symmetrize;,tikhonov=2.0
  	mdgsolvetie,'laplacian',ptie,tikhonov=2.0
	;mdgsolvetie,'inverse gradient',ptie;,/symmetrize 
	
	;display the phase
	phi_m=(*ptie).phase
	bt_x=(*ptie).bxt
	bt_y=(*ptie).byt
	bbt=sqrt(bt_x^2+bt_y^2)

	tvscl,rebin(phi_m,wsi/2,wsi/2),0
	plot,phi_m[*,nd2],/noerase
	phiarr_m[*,*,k-1] = phi_m
	phiarr_e[*,*,k-1] = phi_e
	
	;save the images
	;fn=fd.path+'/images/'+fd.prefix+'_sep'+string(k,format='(I01)')
	fn=flipfile.path+'/images/'+flipfile.prefix+'_sep'+string(k,format='(I01)')
	write_tiff,fn+'_mphase.tiff',(Reverse(phi_m,2)),/float,units=3,xresol=1/delta,yresol=1/delta
	write_tiff,fn+'_ephase.tiff',(Reverse(phi_e,2)),/float,units=3,xresol=1/delta,yresol=1/delta
	write_tiff,fn+'_btx.tiff',(Reverse(bt_x,2)),/float,units=3,xresol=1/delta,yresol=1/delta
	write_tiff,fn+'_bty.tiff',(Reverse(bt_y,2)),/float,units=3,xresol=1/delta,yresol=1/delta
	write_tiff,fn+'_bbt.tiff',(Reverse(bbt,2)),/float,units=3,xresol=1/delta,yresol=1/delta
	write_tiff,fn+'_dIdzE.tiff',(reverse(didz_e,2)),/float,units=3,xresol=1/delta,yresol=1/delta
	write_tiff,fn+'_dIdzM.tiff',(reverse(didz_m,2)),/float,units=3,xresol=1/delta,yresol=1/delta
	colorwheel,dim,dim,bt_x,bt_y,cim,fn+'_colorB.tiff'
	
;	;free pointers
	ptr_free,ptr_valid()
endfor

;store the mag phase data
openw,1,flipfile.path+'/images/'+flipfile.prefix+'_tiesep.mphase'
q=assoc(1,intarr(3))
q[0]=[dim,dim,num/2]
q=assoc(1,fltarr(dim,dim),12)
for kk=1,num/2 do q[kk-1]=phiarr_m[*,*,kk-1]
close,1
openw,1,flipfile.path+'/images/'+flipfile.prefix+'_tiesep.ephase'
q=assoc(1,intarr(3))
q[0]=[dim,dim,num/2]
q=assoc(1,fltarr(dim,dim),12)
for kk=1,num/2 do q[kk-1]=phiarr_e[*,*,kk-1]
close,1

;end of chk_deriv
endif else begin
print,'Computing the phase using polynomial fit for derivative..'


;concatenate the unflip and flip arrays
unstack = [[[reverse(unstack[*,*,0:num/2],3)]],[[unstack[*,*,num/2+1:*]]]]
flstack = [[[reverse(flstack[*,*,0:num/2],3)]],[[flstack[*,*,num/2+1:*]]]]
imstack = [[[unstack]],[[flstack]]]

;going to use poly fit for derivative.
;Linear Focal series
;defval = (findgen(num)-float(num/2))*defstep;*defvals[0]
;Quadratic Focal series
;tt = (findgen(num)-float(num/2))
;defval = tt*abs(tt)*defstep
;Preset Focal series
defval = [[-reverse(defvals)],0,[defvals]]
defstep = defvals[0]

for ii=0,2*num-1 do imstack[*,*,ii] *= mask

;initialize the pointers for getting various parameters
	params=simulationparameters(dim=dim,delta=dd,EE=ee,cs=cs,$
		thetac=thetac,spread=spread,defstep=defstep)
	pimage=params[0]
	pscope_ltem=params[1]
	pphyscon=params[2]
	ptie=params[3]

  ;### make sure the integrated intensity in each image is the same ###
	t = total(total(imstack,1),1)
	print,'image totals : ',t
	tm = max(t)
	t = tm/t
	for kk=0,num*2-1 do imstack[0:*,0:*,kk] *= t[kk]
	t = total(total(imstack,1),1)

	;####scale image intensity to interval [0,1] and make
	; sure that in-focus image is shifted away from zero intensity####
	imax = max(imstack)
	imstack /= imax
	;for kk=0,2 do begin
  ;  temp = reform(imstack[*,*,kk])
  ;  temp /= max(temp)
  ;  temp = (temp-mean(temp))/stddev(temp)
  ;  imstack[*,*,kk] = temp
  ;endfor
  imstack += 0.001D0
  imstack[*,*,num/2] += 1.D0-mask[*,*]
  imstack[*,*,num/2*3] += 1.D0-mask[*,*]
  inf_img = (imstack[*,*,num/2]+imstack[*,*,num/2*3])/2.0

;## compute intensity derivatives
	;dIdz_un = cdderiv(reform(imstack[*,*,0:num-1]),defval,mask=mask,/verbose)*defval*num	
	;dIdz_fl = cdderiv(reform(imstack[*,*,num:*]),defval,mask=mask,/verbose)*defval*num
	;dIdz_m = (dIdz_un-dIdz_fl)/2.0
	;dIdz_e = (dIdz_un+dIdz_fl)/2.0
	dIdz_m = cdderiv(imstack,[defval,reverse(defval)],mask=mask,/verbose)*defstep*(num*2+1)/2.0
	dIdz_e = cdderiv(imstack,[defval,defval],mask=mask,/verbose)*defstep*(num*2+1)/2.0
	;Using pre-calculated dIdz_m and dIdz_e
;	dIdz_m = reverse(read_image(flipfile.path+'images/ftfs2_deriv_dIdzM.tiff'),2);/(num*2+1)*4.0
;	dIdz_e = reverse(read_image(flipfile.path+'images/ftfs2_deriv_dIdzE.tiff'),2);/(num*2+1)*4.0


;## next, make sure that the derivative has zero total "energy" 
	tot = total(dIdz_m)/total(mask)
	dIdz_m -= tot
	tot = total(dIdz_e)/total(mask)
	dIdz_e -= tot
	
	; in-focus image 
	(*ptie).ima_if=inf_img

	;solve for each phase
	;e-phase shift
	(*ptie).ima_dIdz = dIdz_e
	;now perform phase reconstruction
  	mdgsolvetie,'laplacian',ptie;,/symmetrize;,tikhonov=2.0
  	;mdgsolvetie,'laplacian',ptie,tikhonov=2.0
	;mdgsolvetie,'inverse gradient',ptie;,/symmetrize 
	
	;display the phase
	phi_e=(*ptie).phase
	et_x=(*ptie).bxt
	et_y=(*ptie).byt

	tvscl,congrid(phi_e,wsi/2,wsi/2),0
	plot,phi_e[*,nd2],/noerase
	;phiarr[*,*,k-1]=phi

	;m-phase shift
	(*ptie).ima_dIdz = dIdz_m
	;now perform phase reconstruction
  	mdgsolvetie,'laplacian',ptie;,/symmetrize;,tikhonov=2.0
  	;mdgsolvetie,'laplacian',ptie,tikhonov=2.0
	;mdgsolvetie,'inverse gradient',ptie;,/symmetrize 
	
	;display the phase
	phi_m=(*ptie).phase
	bt_x=(*ptie).bxt
	bt_y=(*ptie).byt
	bbt = sqrt(bt_x^2+bt_y^2)

	tvscl,congrid(phi_m,wsi/2,wsi/2),0
	plot,phi_m[*,nd2],/noerase
	
	;stop	
	;save the images
	fn=flipfile.path+'/images/'+flipfile.prefix+'_deriv'
	write_tiff,fn+'_mphase.tiff',(Reverse(phi_m,2)),/float,units=3,xresol=1/delta,yresol=1/delta
	write_tiff,fn+'_ephase.tiff',(Reverse(phi_e,2)),/float,units=3,xresol=1/delta,yresol=1/delta
	write_tiff,fn+'_btx.tiff',(Reverse(bt_x,2)),/float,units=3,xresol=1/delta,yresol=1/delta
	write_tiff,fn+'_bty.tiff',(Reverse(bt_y,2)),/float,units=3,xresol=1/delta,yresol=1/delta
	write_tiff,fn+'_bbt.tiff',(Reverse(bbt,2)),/float,units=3,xresol=1/delta,yresol=1/delta
	write_tiff,fn+'_dIdzE.tiff',(reverse(didz_e,2)),/float,units=3,xresol=1/delta,yresol=1/delta
	write_tiff,fn+'_dIdzM.tiff',(reverse(didz_m,2)),/float,units=3,xresol=1/delta,yresol=1/delta
	resz=64
	set_plot,'ps'
	device,filename=fn+'_vecplotB.eps',/encap,xsi=6,ysi=6,/inch
	velovect,rebin(bt_x,resz,resz),rebin(bt_y,resz,resz),len=2
	device,/close
	set_plot,'x'
	colorwheel,dim,dim,bt_x,bt_y,cim,fn+'_colorB.tiff'

	;store the mag phase data
	openw,1,flipfile.path+'/images/'+flipfile.prefix+'_derivsep.phase'
	q=assoc(1,intarr(2))
	q[0]=[dim,dim,2]
	q=assoc(1,fltarr(dim,dim),12)
	q[0]=phi_e
	q[1]=phi_m
	close,1
	
;	;free pointers
	ptr_free,ptr_valid()

endelse ;end of chk_Deriv

stop

end
