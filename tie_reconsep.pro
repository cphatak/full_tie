@mdgsolvetie
@simulationparameters
@cdderiv
;-------------------------------------------------------------------------------------------------
;+
;NAME: TIE_RECONSEP.PRO
;
;CALLIN SEQUENCE: result = tie_reconsep(unflip=unflip,flip=flip,data=data)
;
;PURPOSE: This is the modified tie_reconsep as a function to work with full_tie.pro
;	  This routine will perform TIE based phase retrieval using the modified 
;	  reconstruction method (Ultra paper 2014). There are options to use 3 point
;	  derivative or SVD fitting to polynomial. The routine will compute the phase
;	  shift and automatically save the images   (32-bit TIFF files) of phase, Bxt, 
;	  Byt, BB, inf_image, colorplot. The actual reconstructed phase will be saved
;	  in an assoc IDL file for later usage and reference.  
;
;REQUIRES: mdgsolvetie.pro - core reconstruction routine
;	   simulationparameters.pro - to get the values of physical parameters
;	   cdderiv.pro - requried to compute the derivative using SVD fitting method.
;
;VERSION: 3.0 - to work with full_tie.pro
;
;AUTHOR: CD Phatak, ANL, Sep 8th, 2014.
;
;-
;-------------------------------------------------------------------------------------------------

function TIE_reconsep,unflip=unflip,flip=flip,data=data

dim=unflip.dim ;2048
d2=dim/2
scl=1.0
ndim=dim/scl
nd2=ndim/2
close,/all

;defocus values
defvals = unflip.defvals;[3133440]*0.18

;No. of images
num=unflip.num

;get the scope paramters and TIE pointers
;dd=4.13 ;pixel size in nm
dd=unflip.delta;*1000.0
delta=dd;*1000.
EE=data.ee ;V
cs=data.cs ;nm
cc=4.0D6 ;nm
thetac=0.0001D0 ;rad
sigV=1.0D-6 ;sig(V)/V
sigI=1.0D-6 ;sig(I)/I
sigE=2.0D-6 ;sig(E)/E
spread= cc * sqrt(sigV^2 + sigI^2 + sigE^2)

;TIE recon. parameters
set=data.tiesetting
if (data.tiemethod eq 0) then met='laplacian' else met = 'inverse gradient'

wsi=1024
window,0,xsi=wsi,ysi=wsi,retain=2;,xpos=2000
;array to store all the phase

;rotational angle
rot_ang = 0.0

;load all the tfseries images
;unflipped
openr,1,unflip.path+unflip.prefix+'.decon'
q=assoc(1,fltarr(dim,dim),12)
unstack=fltarr(dim,dim,num)
for ii=0,num-1 do unstack[*,*,ii]=(rot(q[ii],rot_ang,cubic=-0.5,missing=0.0))
unflip_mask = rot(q[num],rot_ang,missing=0.0)
close,1

;flipped
openr,1,flip.path+flip.prefix+'_final.aligned'
q=assoc(1,fltarr(dim,dim),12)
flstack = fltarr(dim,dim,num)
for ii=0,num-1 do flstack[*,*,ii]=rot(q[ii],rot_ang,cubic=-0.5,missing=0.0)
flip_mask = rot(q[num],rot_ang,missing=0.0)
close,1

;common mask
mask = unflip_mask*flip_mask
mask=erode(mask,replicate(1.0,15,15))
;select the central contiguous area
ll=label_Region(mask)
mask=ll eq ll[d2,d2]

;stop

;option to select a smaller 2K by 2K area from the images - 
ans=''
cen=0
;change to smaller dimensions..
ndim=1024*2
nd2=ndim/2
;xs=0 & xe=dim-1
;ys=0 & ye=dim-1
read,prompt='Do you want to select a smaller '+strtrim(string(ndim),2)+' by '+strtrim(string(ndim),2)+' region? (y/n): ',ans
if (ans eq 'y') then begin
 read,prompt='Do you want to select a central region (0) or custom region (1) : ',cen
  if (cen eq 1) then begin
    cen_select=1
  	read,prompt='Do you want to enter coordinates for center(0) or select by clicking(1): ',cen_select
  	if (cen_select eq 0) then begin
  	xx=0.0 & yy=0.0
  	read,prompt='Enter coordinate for x: ',xx
  	read,prompt='Enter coordinate for y: ',yy
  	tvscl,rebin(unstack[*,*,0]+flstack[*,*,0],wsi,wsi)
  	xx /= 2.0 & yy /= 2.0
  	plots,xx,yy,/dev,psym=7,symsize=2,color=200
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
     cc=tvrd(/order,true=1)
     endif else begin
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
   endwhile
   endelse
   fld_inf = file_info(flip.path+'/images/')
    if (fld_inf.exists ne 1) then spawn,'mkdir '+flip.path+'/images/'
    write_tiff,flip.path+'/images/'+flip.prefix+'_selarea.tiff',cc
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
if (data.deriv eq 0) then chk_deriv='n' else chk_deriv='y'
;read,prompt='Do you want to use polynomial fit for calculating intensity derivative? (y/n):',chk_deriv

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
  	;mdgsolvetie,'laplacian',ptie,tikhonov=2.0
	;mdgsolvetie,'inverse gradient',ptie;,/symmetrize
	if (set[0]*set[1] eq 1) then mdgsolvetie,met,ptie,tikhonov=data.tikval,/symmetrize
	if (set[0]+set[1] eq 0) then mdgsolvetie,met,ptie else begin
		if (set[0] eq 0) then mdgsolvetie,met,ptie,/symmetrize
		if (set[1] eq 0) then mdgsolvetie,met,ptie,tikhonov=data.tikval
	endelse
	
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
  	;mdgsolvetie,'laplacian',ptie,tikhonov=2.0
	;mdgsolvetie,'inverse gradient',ptie;,/symmetrize 
	if (set[0]*set[1] eq 1) then mdgsolvetie,met,ptie,tikhonov=data.tikval,/symmetrize
	if (set[0]+set[1] eq 0) then mdgsolvetie,met,ptie else begin
		if (set[0] eq 0) then mdgsolvetie,met,ptie,/symmetrize
		if (set[1] eq 0) then mdgsolvetie,met,ptie,tikhonov=data.tikval
	endelse
	
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
	fld_inf = file_info(flip.path+'/images/')
	if (fld_inf.exists ne 1) then spawn,'mkdir '+flip.path+'/images/'
	;fn=fd.path+'/images/'+fd.prefix+'_sep'+string(k,format='(I01)')
	fn=flip.path+'/images/'+flip.prefix+'_sep'+string(k,format='(I01)')
	write_tiff,fn+'_mphase.tiff',(Reverse(phi_m,2)),/float,units=3,xresol=1/delta,yresol=1/delta
	write_tiff,fn+'_ephase.tiff',(Reverse(phi_e,2)),/float,units=3,xresol=1/delta,yresol=1/delta
	write_tiff,fn+'_btx.tiff',(Reverse(bt_x,2)),/float,units=3,xresol=1/delta,yresol=1/delta
	write_tiff,fn+'_bty.tiff',(Reverse(bt_y,2)),/float,units=3,xresol=1/delta,yresol=1/delta
	write_tiff,fn+'_bbt.tiff',(Reverse(bbt,2)),/float,units=3,xresol=1/delta,yresol=1/delta
	write_tiff,fn+'_dIdzE.tiff',(reverse(didz_e,2)),/float,units=3,xresol=1/delta,yresol=1/delta
	write_tiff,fn+'_dIdzM.tiff',(reverse(didz_m,2)),/float,units=3,xresol=1/delta,yresol=1/delta
	if (k eq 1) then write_tiff,fn+'_infimg.tiff',Reverse(reform(imstack[*,*,2])*mask,2),/float,units=3,xresol=1/delta,yresol=1/delta
	colorwheel,dim,dim,bt_x,bt_y,cim,fn+'_colorB.tiff'
    ;stop
	
;	;free pointers
	ptr_free,ptr_valid()
endfor

;store the mag phase data
openw,1,flip.path+'/images/'+flip.prefix+'_tiesep.mphase'
q=assoc(1,intarr(3))
q[0]=[dim,dim,num/2]
q=assoc(1,fltarr(dim,dim),12)
for kk=1,num/2 do q[kk-1]=phiarr_m[*,*,kk-1]
close,1
openw,1,flip.path+'/images/'+flip.prefix+'_tiesep.ephase'
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
;stop
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
	deriv_file_info = file_info(flip.path+'/images/'+flip.prefix+'_deriv_dIdzM.tiff')
	deriv_file_ans='n'
	if (deriv_file_info.exists eq 1) then read,prompt='PolyFit Intensity derivative exists. Use existing data (y/n):',deriv_file_ans
	if (deriv_file_ans eq 'n') then begin
	;dIdz_un = cdderiv(reform(imstack[*,*,0:num-1]),defval,mask=mask,/verbose)*defval*num	
	;dIdz_fl = cdderiv(reform(imstack[*,*,num:*]),defval,mask=mask,/verbose)*defval*num
	;dIdz_m = (dIdz_un-dIdz_fl)/2.0
	;dIdz_e = (dIdz_un+dIdz_fl)/2.0
	dIdz_m = cdderiv(imstack,[defval,reverse(defval)],mask=mask,/verbose)*defstep*(num*2+1)/2.0
	dIdz_e = cdderiv(imstack,[defval,defval],mask=mask,/verbose)*defstep*(num*2+1)/2.0
	endif else begin
	;Using pre-calculated dIdz_m and dIdz_e
	dIdz_m = reverse(read_image(flip.path+'/images/'+flip.prefix+'_deriv_dIdzM.tiff'),2);/(num*2+1)*4.0
	dIdz_e = reverse(read_image(flip.path+'/images/'+flip.prefix+'_deriv_dIdzE.tiff'),2);/(num*2+1)*4.0
	endelse


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
  	;mdgsolvetie,'laplacian',ptie;,/symmetrize;,tikhonov=2.0
  	;mdgsolvetie,'laplacian',ptie,tikhonov=2.0
	;mdgsolvetie,'inverse gradient',ptie;,/symmetrize
	if (set[0]*set[1] eq 1) then mdgsolvetie,met,ptie,tikhonov=data.tikval,/symmetrize
	if (set[0]+set[1] eq 0) then mdgsolvetie,met,ptie else begin
		if (set[0] eq 0) then mdgsolvetie,met,ptie,/symmetrize
		if (set[1] eq 0) then mdgsolvetie,met,ptie,tikhonov=data.tikval
	endelse
	
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
  	;mdgsolvetie,'laplacian',ptie;,/symmetrize;,tikhonov=2.0
  	;mdgsolvetie,'laplacian',ptie,tikhonov=2.0
	;mdgsolvetie,'inverse gradient',ptie;,/symmetrize 
	if (set[0]*set[1] eq 1) then mdgsolvetie,met,ptie,tikhonov=data.tikval,/symmetrize
	if (set[0]+set[1] eq 0) then mdgsolvetie,met,ptie else begin
		if (set[0] eq 0) then mdgsolvetie,met,ptie,/symmetrize
		if (set[1] eq 0) then mdgsolvetie,met,ptie,tikhonov=data.tikval
	endelse
	
	;display the phase
	phi_m=(*ptie).phase
	bt_x=(*ptie).bxt
	bt_y=(*ptie).byt
	bbt = sqrt(bt_x^2+bt_y^2)

	tvscl,congrid(phi_m,wsi/2,wsi/2),0
	plot,phi_m[*,nd2],/noerase
	
	;stop	
	;save the images
	fld_inf = file_info(flip.path+'/images/')
	if (fld_inf.exists ne 1) then spawn,'mkdir '+flip.path+'/images/'
	fn=flip.path+'/images/'+flip.prefix+'_deriv'
	write_tiff,fn+'_mphase.tiff',(Reverse(phi_m,2)),/float,units=3,xresol=1/delta,yresol=1/delta
	write_tiff,fn+'_ephase.tiff',(Reverse(phi_e,2)),/float,units=3,xresol=1/delta,yresol=1/delta
	write_tiff,fn+'_btx.tiff',(Reverse(bt_x,2)),/float,units=3,xresol=1/delta,yresol=1/delta
	write_tiff,fn+'_bty.tiff',(Reverse(bt_y,2)),/float,units=3,xresol=1/delta,yresol=1/delta
	write_tiff,fn+'_bbt.tiff',(Reverse(bbt,2)),/float,units=3,xresol=1/delta,yresol=1/delta
	write_tiff,fn+'_dIdzE.tiff',(reverse(didz_e,2)),/float,units=3,xresol=1/delta,yresol=1/delta
	write_tiff,fn+'_dIdzM.tiff',(reverse(didz_m,2)),/float,units=3,xresol=1/delta,yresol=1/delta
	write_tiff,fn+'_infimg.tiff',Reverse(reform(inf_img)*mask,2),/float,units=3,xresol=1/delta,yresol=1/delta
	resz=64
	set_plot,'ps'
	device,filename=fn+'_vecplotB.eps',/encap,xsi=6,ysi=6,/inch
	velovect,rebin(bt_x,resz,resz),rebin(bt_y,resz,resz),len=2
	device,/close
	set_plot,'x'
	colorwheel,dim,dim,bt_x,bt_y,cim,fn+'_colorB.tiff'

	;store the mag phase data
	openw,1,flip.path+'/images/'+flip.prefix+'_derivsep.phase'
	q=assoc(1,intarr(2))
	q[0]=[dim,dim,2]
	q=assoc(1,fltarr(dim,dim),12)
	q[0]=phi_e
	q[1]=phi_m
	close,1
	
;	;free pointers
	ptr_free,ptr_valid()

endelse ;end of chk_Deriv
wdelete,0
return,1
end
