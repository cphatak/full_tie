@mdgsolvetie
;+
;program to perform TIE reconstructions on series of 
;through focus images at different foci. This uses the
;mdgsolvetie to compute the phase and store the phase for
;each series as a data file. The image and the surfaceplot
;are also saved alongwith.
;@author, Charudatta Phatak, CMU, 12/04/08.
;-

pro TIE_recon,fd,flipfile=flipfile,dcdfile=dcdfile,single=single

dim=fd.dim ;2048
d2=dim/2
scl=1.0
ndim=dim/scl
nd2=ndim/2
close,/all

;defocus values
;defstepsize = 90.0*0.18 ;nm
;defstep value
;defstep = 10.0 ;series 1
;defvals=[1013760]*0.18
;defvals = [103680,207360,506880]*0.18
;defvals = [28800,51840,92160]*0.18
;defvals = [11520,23040,40320]*0.18
;defvals = [1105920,3133440]*0.18
;defvals = [11520,25920]*0.18
defvals=[1000*0.18]
defstep=fd.defstep

;No. of images
num=fd.num

;get the scope paramters and TIE pointers
dd=fd.delta*1000 ;pixel size in nm
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

;load all the tfseries images
;unflipped
openr,1,fd.path+fd.prefix+'.aligned'
q=assoc(1,fltarr(dim,dim),12)
unstack=fltarr(dim,dim,num)
for ii=0,num-1 do unstack[*,*,ii]=q[ii]
unflip_mask = q[num]
close,1

;flipped
flip_mask = replicate(1.0,dim,dim)
if (~keyword_Set(single)) then begin
openr,1,flipfile.path+flipfile.prefix+'_final.aligned'
q=assoc(1,fltarr(dim,dim),12)
flstack = fltarr(dim,dim,num)
for ii=0,num-1 do flstack[*,*,ii]=q[ii]
flip_mask = q[num]
close,1
endif



;if using with DCD
if (keyword_Set(dcdfile)) then begin
  openr,10,dcdfile.path+dcdfile.prefix+'_final.aligned'
  q=assoc(10,fltarr(dim,dim),12)
  dcd_mask=q[num]
  close,1
endif else dcd_mask = replicate(1.0,dim,dim)

;common mask
mask = unflip_mask*flip_mask*dcd_mask
mask=erode(mask,replicate(1.0,9,9))
;select the central most contigous region.
lmask=label_Region(mask)
mask=lmask eq lmask[d2,d2]

;selecting smaller area
if (~keyword_Set(single)) then image = (unstack[*,*,0]+flstack[*,*,0])/2.0*mask else $
	image = unstack[*,*,0]*mask
ans=''
cen=0
ndim=1024
nd2=ndim/2
;xs=0 & xe=dim-1
;ys=0 & ye=dim-1
read,prompt='Do you want to select a smaller 1K by 1K region? (y/n): ',ans
if (ans eq 'y') then begin
 read,prompt='Do you want to select a central region (0) or custom region (1) : ',cen
  if (cen eq 1) then begin
    print,'Select the center of the 512 by 512 area'
    sat=0
    while (sat eq 0) do begin
		;erase
    tvscl,rebin(image,wsi,wsi)
    cursor,xx,yy,4,/dev
    ;xx=1200.0/2 & yy=1386.00/2.0
    xx = 727 & yy= 683
    plots,xx,yy,/dev,psym=7,symsize=2,color=200
    xyouts,xx+10,yy+10,string(xx)+','+string(yy),/dev,color=200,charsize=2.0
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
     read,prompt='Are you satisfied with the selection (0/1):',sat
    cc=tvrd(/order,true=1)
   endwhile
 ;xs*=2 & xe*=2 & ys*=2 &ye*=2
 endif else begin
	xs = d2-nd2 & xe = d2+nd2-1
	ys = d2-nd2 & ye = d2+nd2-1
endelse
unstack=unstack[xs:xe,ys:ye,0:*]
if (~keyword_set(single)) then flstack=flstack[xs:xe,ys:ye,0:*]
mask=mask[xs:xe,ys:ye,0:*]
dim=ndim
d2=dim/2
endif
print,dim

;check if user wants to use CDDERIV or not.
chk_deriv=''
read,prompt='Do you want to use polynomial fit for calculating intensity derivative? (y/n):',chk_deriv

;arrays to hold the phase shifts
unphiarr=fltarr(dim,dim,num/2)
flphiarr=fltarr(dim,dim,num/2)

;check if dcd_file option is set. If set, then we are directly 
;subtracting previously previously calculated e-phase from 
;the total phase shift. No need to reconstruct the flipped
;phase shift.

if (chk_deriv eq 'n') then begin

if (~keyword_set(dcdfile)) then begin
;Loop to start reconstruction for each defocus value
;reconstruction first for unflipped images
for k=1,num/2 do begin
  	defval = k^2*defstep;*defstepsize
	;defval = defvals[k-1]
	print,'Reconstructing phase for defocus =',defval

	;load the images
	imstack=fltarr(dim,dim,3)
	;underfocus
	imstack[*,*,0] = unstack[*,*,k]*mask
	;infocus
	imstack[*,*,1] = unstack[*,*,0]*mask
	;overfocus
	imstack[*,*,2] = unstack[*,*,num/2+k]*mask
	

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
	for kk=0,2 do imstack[0:*,0:*,kk] *= t[kk]
	t = total(total(imstack,1),1)
	print,'image totals : ',t

	;####scale image intensity to interval [0,1] and make
	; sure that in-focus image is shifted away from zero intensity####
  ;stop
	imax = max(imstack)
	imstack /= imax
  	imstack += 0.001D0
	imstack[0:*,0:*,1] += 1.D0-mask[*,*]

	; for now, simply use two images to get to the derivative
	dIdz = reform(imstack[*,*,2]) - reform(imstack[*,*,0])
	;using cdderiv to estimate derivative
	(*ptie).ima_dIdz = dIdz

	tvscl,rebin((*ptie).ima_dIdz,wsi,wsi)
	tvscl,rebin(imstack[*,*,0],wsi/4,wsi/4),0
	tvscl,rebin(imstack[*,*,1],wsi/4,wsi/4),1
	tvscl,rebin(imstack[*,*,2],wsi/4,wsi/4),2

	; next, make sure that the derivative has zero total "energy" 
	tot = total((*ptie).ima_dIdz)/total(mask)
	(*ptie).ima_dIdz -= tot
	;print,'total dIdz = ',total((*ptie).ima_dIdz)

	; in-focus image 

	(*ptie).ima_if=reform(imstack[0:*,0:*,1],dim,dim)
	
	;now perform phase reconstruction
	;mdgsolvetie,'laplacian',ptie;,/symmetrize;,tikhonov=2.0
 	mdgsolvetie,'laplacian',ptie,tikhonov=0.5
;	mdgsolvetie,'inverse gradient',ptie,/symmetrize 
	
	;display the phase
	phi=(*ptie).phase
	bxt=(*ptie).bxt
	byt=(*ptie).byt

	tvscl,rebin(phi,wsi,wsi)
	plot,phi[*,nd2],/noerase
	unphiarr[*,*,k-1]=phi

	;write the image and surface plot
	;if (defval lt 1000) then fn_def=string(defval,format='(I3)')
	;if (defval gt 1000 AND defval lt 10000) then fn_def=string(defval,format='(I4)')
	;if (defval gt 10000 AND defval lt 100000) then fn_def=string(defval,format='(I5)')
	;if (defval gt 100000 AND defval lt 1000000) then fn_def=string(defval,format='(I6)')
	;;phase_image
	;write_tiff,fd.path+'/images/'+fd.prefix+fn_def+'_phi.tiff',bytscl(reverse(phi,2))
	;write_tiff,fd.path+'/images/'+fd.prefix+fn_def+'_bxt.tiff',bytscl(reverse(bxt,2))
	;write_tiff,fd.path+'/images/'+fd.prefix+fn_def+'_byt.tiff',bytscl(reverse(byt,2))
;
;	;surface_plot
;	set_plot,'ps'
;	device,filename=fd.path+'/images/'+fd.prefix+fn_def+'.eps',/encapsulated,$
;		xsi=6,ysi=6,/inches
;	shade_surf,phi,charsize=2,pixels=2048
;	device,/close
;	set_plot,'x'
;	;colorimage of B
;	fnameB = fd.path+'/images/'+fd.prefix+fn_def+'_color.tiff'
;	colorwheel,dim,dim,bxt,byt,cim,fnameB
;	cdim = size(cim,/dim)
;	tvscl,congrid(cim,cdim[0],cdim[1]/4,cdim[2]/4),2,true=1
	
;	;free pointers
	ptr_free,ptr_valid()
	
endfor
;store the phase data
openw,1,fd.path+fd.prefix+'.phase'
q=assoc(1,intarr(3))
q[0]=[dim,dim,num/2]
q=assoc(1,fltarr(dim,dim),12)
for kk=1,num/2 do q[kk-1]=unphiarr[*,*,kk-1]
close,1

endif ;dcd_file check
endif else begin;chk_deriv check

print,'Computing the phase using polynomial fit for derivative..'
;going to use poly fit for derivative.
nn = (findgen(num)-float(num/2))
;Quadratic series
defvals = nn^2*defstep*(nn/abs(nn))
;Linear Series
;defvals = nn*defstep
defvals[num/2]=0.0
;concatenate the unflip array
imstack = [[[reverse(unstack[*,*,0:num/2],3)]],[[unstack[*,*,num/2+1:*]]]]
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
	for kk=0,num-1 do imstack[0:*,0:*,kk] *= t[kk]
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
  inf_img = (imstack[*,*,num/2])
;Compute the longitudinal derivative using polynomial fitting
dIdz = cdderiv(imstack,defvals,mask=mask,/verbose)*defstep*num/2.0
;dIdz = reverse(read_tiff(fd.path+'images/'+fd.prefix+'_deriv_dIdz.tiff'),2)
(*ptie).ima_dIdz = dIdz

	; next, make sure that the derivative has zero total "energy" 
	tot = total((*ptie).ima_dIdz)/total(mask)
	(*ptie).ima_dIdz -= tot
	;print,'total dIdz = ',total((*ptie).ima_dIdz)

	; in-focus image 

	(*ptie).ima_if=reform(imstack[0:*,0:*,1],dim,dim)
	
	;now perform phase reconstruction
	mdgsolvetie,'laplacian',ptie;,tikhonov=2.0
 	;mdgsolvetie,'laplacian',ptie,tikhonov=2.0
	;mdgsolvetie,'inverse gradient',ptie,/symmetrize 
	
	;display the phase
	phi=(*ptie).phase
	bxt=(*ptie).bxt
	byt=(*ptie).byt

	tvscl,rebin(phi,wsi,wsi)
	plot,phi[*,nd2],/noerase
;store the phase data
openw,1,fd.path+fd.prefix+'_derivfit.phase'
q=assoc(1,intarr(2))
q=[dim,dim]
q=assoc(1,fltarr(dim,dim))
q[0]=phi
close,1

fn=fd.path+'images/'+fd.prefix
write_tiff,fn+'_deriv_phase.tiff',reverse(phi,2),/float,units=3,xresol=1/delta,yresol=1/delta
write_tiff,fn+'_deriv_dIdz.tiff',reverse(dIdz,2),/float,units=3,xresol=1/delta,yresol=1/delta
write_tiff,fn+'_deriv_bxt.tiff',reverse(bxt,2),/float,units=3,xresol=1/delta,yresol=1/delta
write_tiff,fn+'_deriv_byt.tiff',reverse(byt,2),/float,units=3,xresol=1/delta,yresol=1/delta

endelse ;chk_Deriv end


;now reconstruct the flipped phase
if (~keyword_set(single)) then begin
for k=1,num/2 do begin
  ;defval = k*defstep*defstepsize
	defval = defvals[k-1]
	print,'Reconstructing phase for defocus =',defval

	;load the images
	imstack=fltarr(dim,dim,3)
	;underfocus
	imstack[*,*,0] = flstack[*,*,k]*mask
	;infocus
	imstack[*,*,1] = flstack[*,*,0]*mask
	;overfocus
	imstack[*,*,2] = flstack[*,*,num/2+k]*mask
	

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
	for kk=0,2 do imstack[0:*,0:*,kk] *= t[kk]
	t = total(total(imstack,1),1)
	print,'image totals : ',t

	;####scale image intensity to interval [0,1] and make
	; sure that in-focus image is shifted away from zero intensity####
  ;stop
	imax = max(imstack)
	imstack /= imax
  	imstack += 0.001D0
	imstack[0:*,0:*,1] += 1.D0-mask[*,*]

	; for now, simply use two images to get to the derivative
	dIdz = reform(imstack[*,*,2]) - reform(imstack[*,*,0])
	;using cdderiv to estimate derivative
	(*ptie).ima_dIdz = dIdz

	tvscl,rebin((*ptie).ima_dIdz,wsi,wsi)
	tvscl,rebin(imstack[*,*,0],wsi/4,wsi/4),0
	tvscl,rebin(imstack[*,*,1],wsi/4,wsi/4),1
	tvscl,rebin(imstack[*,*,2],wsi/4,wsi/4),2

	; next, make sure that the derivative has zero total "energy" 
	tot = total((*ptie).ima_dIdz)/total(mask)
	(*ptie).ima_dIdz -= tot
	;print,'total dIdz = ',total((*ptie).ima_dIdz)

	; in-focus image 

	(*ptie).ima_if=reform(imstack[0:*,0:*,1],dim,dim)
	
	;now perform phase reconstruction
	;mdgsolvetie,'laplacian',ptie,/symmetrize;,tikhonov=2.0
 	mdgsolvetie,'laplacian',ptie,tikhonov=2.0
	;mdgsolvetie,'inverse gradient',ptie,/symmetrize 
	
	;display the phase
	phi=(*ptie).phase
	bxt=(*ptie).bxt
	byt=(*ptie).byt

	tvscl,rebin(phi,wsi,wsi)
	plot,phi[*,nd2],/noerase
	flphiarr[*,*,k-1]=phi

	;write the image and surface plot
	;if (defval lt 1000) then fn_def=string(defval,format='(I3)')
	;if (defval gt 1000 AND defval lt 10000) then fn_def=string(defval,format='(I4)')
	;if (defval gt 10000 AND defval lt 100000) then fn_def=string(defval,format='(I5)')
	;if (defval gt 100000 AND defval lt 1000000) then fn_def=string(defval,format='(I6)')
	;;phase_image
	;write_tiff,fd.path+'/images/'+fd.prefix+fn_def+'_phi.tiff',bytscl(reverse(phi,2))
	;write_tiff,fd.path+'/images/'+fd.prefix+fn_def+'_bxt.tiff',bytscl(reverse(bxt,2))
	;write_tiff,fd.path+'/images/'+fd.prefix+fn_def+'_byt.tiff',bytscl(reverse(byt,2))
;
;	;surface_plot
;	set_plot,'ps'
;	device,filename=fd.path+'/images/'+fd.prefix+fn_def+'.eps',/encapsulated,$
;		xsi=6,ysi=6,/inches
;	shade_surf,phi,charsize=2,pixels=2048
;	device,/close
;	set_plot,'x'
;	;colorimage of B
;	fnameB = fd.path+'/images/'+fd.prefix+fn_def+'_color.tiff'
;	colorwheel,dim,dim,bxt,byt,cim,fnameB
;	cdim = size(cim,/dim)
;	tvscl,congrid(cim,cdim[0],cdim[1]/4,cdim[2]/4),2,true=1
	
;	;free pointers
	ptr_free,ptr_valid()
	
endfor
;store the phase data
openw,1,flipfile.path+flipfile.prefix+'.phase'
q=assoc(1,intarr(3))
q[0]=[dim,dim,num/2]
q=assoc(1,fltarr(dim,dim),12)
for kk=1,num/2 do q[kk-1]=flphiarr[*,*,kk-1]
close,1
endif ; single tf-series 

;reconstructing the separate phase shifts using DCD_file input
if (keyword_set(dcdfile))then begin

;load the e-phase calculated previously
;openr,1,dcdfile.path+dcdfile.prefix+'_tiesep.ephase'
openr,1,dcdfile.path+dcdfile.prefix+'_ephi.data'
q=assoc(1,fltarr(dim,dim),12)
phi_e = q[0]
close,1

;compute the m-phase shift
phi_m = reform(flphiarr[*,*,0])-phi_e;*25.0

;gradient kernels
kx=([[0.0,0.0,0.0],[-1.0,0.0,1.0],[0.0,0.0,0.0]])
ky=([[0.0,1.0,0.0],[0.0,0.0,0.0],[0.0,-1.0,0.0]])

;compute mag. maps
bxt = convol(reform(phi_m[*,*]),ky,/center,/edge_truncate)
byt = convol(reform(phi_m[*,*]),kx,/center,/edge_truncate)

;save the images
fn=flipfile.path+'/images/'+flipfile.prefix+'_dcd'+string(k,format='(I01)')
	write_tiff,fn+'_mphase.tiff',(Reverse(phi_m,2)),/float,units=3,xresol=1/delta,yresol=1/delta
	write_tiff,fn+'_ephase.tiff',(Reverse(phi_e,2)),/float,units=3,xresol=1/delta,yresol=1/delta
	write_tiff,fn+'_btx.tiff',(Reverse(bxt,2)),/float,units=3,xresol=1/delta,yresol=1/delta
	write_tiff,fn+'_bty.tiff',(Reverse(byt,2)),/float,units=3,xresol=1/delta,yresol=1/delta
	colorwheel,dim,dim,bxt,byt,cim,fn+'_colorB.tiff'

endif

stop
end
