;pro align_series
;+
; NAME: align_series
;     
; PURPOSE: to align a Lorentz through-focus series
;     
;
; CALLING SEQUENCE: align_series
;     
;
; OPTIONAL INPUT: none
;
; OPTIONAL KEYWORD INPUT: none
; 
; METHOD: amoeba image alignment with DoG filter for translational and stretch parameters
;
; NOTES: 
;
; PROCEDURES/FUNCTIONS USED: widget-based, with simple user interface
;   
; PROJECT: DOE
;
; REVISION HISTORY
;      Written, MDG, 5/7/07
;-

;############################################
function getcm1,param
;############################################

common cmtarget,jj,st,bintx,binty,a,sz,rmask,blarge
common data_common, data 

if (data.nparams eq 4) then begin
  bl = rot(blarge,param[2],param[3],missing=0.0,cubic=-0.5,/pivot)
  c = interp2d(bl,bintx,binty,bintx+param[0],binty + param[1],cubic=-0.5,missing=0.0)
endif else begin
  bl=rot(blarge,param[2],missing=0.0,cubic=-0.5)
 c = interp2d(bl,bintx,binty,bintx*param[3]+binty*param[5]+param[0],$
 bintx*param[6]+binty*param[4]+param[1],missing=0.0,cubic=-0.5)
endelse
 
cc = smooth(c,3)
aa = smooth(a,3)

 ;no mag. & rot. correction
 ;bl = rot(blarge,0.0,1.0,missing=0.0,cubic=-0.5,/pivot)
 ;c = interpolate(bl,bintx+param[0],binty + param[1],cubic=-0.5,missing=0.0)

if (data.method eq 0) then begin
 ; compute the DOG filtered image
 cc = mdg_filter_2d(cc,'LoG',sigma=data.sigma)
 if (data.gtlt eq 1) then cc = float(cc GT data.thrs) * data.filter else cc = float(cc LT data.thrs) * data.filter
  aa = mdg_filter_2d(aa,'LoG',sigma=data.sigma)
 if (data.gtlt eq 1) then aa = float(aa GT data.thrs) * data.filter else aa = float(aa LT data.thrs) * data.filter

 ;cm = total(cc*aa)
 cm = -total((cc-aa)^2)
 endif else begin
 ;Using Mutual Information to align images
 ;use of mask
 ;mask_a = erode(a lt 1000, replicate(1,10,10))
 ;mask_c = erode(c lt 1000, replicate(1,10,10))
;stop

 aa=(aa-mean(aa))/stddev(aa)
 cc=(cc-mean(cc))/stddev(cc)
 ;aa=a
 ;cc=c
 imgs=bytarr(2,data.dim,data.dim)
 imgs[0,*,*] = bytscl(aa*data.filter)
 imgs[1,*,*] = bytscl(cc*data.filter)
 res = joint_entropy(imgs,mutual=cm,/nozero)
 endelse

 tvscl,rebin(0.5*(cc+aa)*data.filter,data.dim,data.dim)

return,-cm
end

;############################################
pro amdog,mode,seriesid=seriesid,imageid=imageid,shifts=shifts
;############################################
;
; use amoeba to align images...
;
common fd_common,fd
common cmtarget,jj,st,bintx,binty,a,sz,rmask,blarge
common data_common, data 
common widget_common,widget
common cm3,yratios
;
; arrays needed for interpolate function
intx = findgen(data.dim)
inty = intx

;line = findgen(data.dim)
;bintx = replicate(1.0,data.dim)##line
;line = findgen(data.dim)
;binty = replicate(1.0,data.dim)#line 
line = findgen(data.dim)-float(data.dim/2)
bintx = replicate(1.0,data.dim)##line
;line = findgen(data.dim)
binty = replicate(1.0,data.dim)#line 


WIDGET_CONTROL, SET_VALUE='Starting image alignment for image '+string(imageid), widget.status

;erode mask
strmsk=replicate(1.0,3,3)

; align each tilted image with respect to image a

;data.fits=fltarr(data.max_nota,4)
;data.fits[0,0:*] = [0.0,0.0,0.0,1.0]

openr,15,fd.path+fd.prefix+'.decon'
q=assoc(15,fltarr(data.dim,data.dim),12)
if (imageid eq (fd.num+1)/2) then a = q[0] else a = q[imageid-1]
;a = q[0]
close,15
wset,data.drawID
  ;a = float(erode(bytscl(a),strmsk,/gray))
  c = smooth(a,3)
  aa = c
  if (data.gtlt eq 1) then adog = float(mdg_filter_2d(c,'LoG',sigma=data.sigma) GT data.thrs)*data.filter else adog = float(mdg_filter_2d(c,'LoG',sigma=data.sigma) LT data.thrs)*data.filter
  ga = adog
  mean = total(ga)/float(data.dim)^2
  sd   = sqrt(total((ga-mean)^2))/float(data.dim)
  ga   = (ga-mean)/sd
  fa   = fft(ga*data.filter,-1)
; next, load the current image 
  openr,15,fd.path+fd.prefix+'.decon'
  q=assoc(15,fltarr(data.dim,data.dim),12)
  blarge = q[imageid]
  close,15
  ;blarge = float(erode(bytscl(blarge),strmsk,/gray))
  c = smooth(blarge,3)
  if (data.gtlt eq 1) then gb = float(mdg_filter_2d(c,'LoG',sigma=data.sigma) GT data.thrs)*data.filter else gb = float(mdg_filter_2d(c,'LoG',sigma=data.sigma) LT data.thrs)*data.filter
  mean = total(gb)/float(data.dim)^2
  sd   = sqrt(total((gb-mean)^2))/float(data.dim)
  gb    = (gb-mean)/sd
; then do a cross correlation 
  fb = conj(fft(gb * data.filter,-1))
  ab = smooth(abs(fft(fa*fb,1)),5)
  tvscl,rebin(shift(alog(ab+0.01),data.dim/2,data.dim/2),data.dim,data.dim)

; determine where the maximum is
  qmax = where (ab eq max(ab))
  roughx = qmax(0) mod data.dim
  roughy = qmax(0)/data.dim
  if (roughx gt data.dim/2) then roughx = roughx - data.dim
  if (roughy gt data.dim/2) then roughy = roughy - data.dim
  WIDGET_CONTROL,SET_VALUE='Image '+string(imageid,format="(I6)")+' rough alignment parameters ['+string(roughx)+','+string(roughy)+']', widget.status
  plots,(data.dim/2+[roughx-3,roughx+3])/2,(data.dim/2+[roughy,roughy])/2,/dev,color=230
  plots,(data.dim/2+[roughx,roughx])/2,(data.dim/2+[roughy-3,roughy+3])/2,/dev,color=230
  
  ftol = 1.0e-4
  ; transx, transy, rot, scale
  roughx = 0 & roughy = 0
  print,imageid,'Rough shifts=',-roughx,-roughy
  if (data.nparams eq 4) then begin
  pa = [-roughx,-roughy,0.0,1.0] 
  scl = [5.0,5.0,5.0,0.5]
  endif else begin
  pa=[-roughx,-roughy,0.0,1.0,1.0,0.0,0.0]
  scl=[5.0,5.0,5.0,0.5,0.5,0.05,0.05]
  endelse


  jj = imageid
  iter = 0
  itmax = 3000
  fmin = amoeba(ftol,function_name='getcm1',function_value=fv, ncalls = iter, nmax = itmax, p0=pa,scale=scl)
  print,imageid,' iterations = ',iter
  print,imageid,' parameter fit =',fmin
  ;data.fits[0,0:3] = fmin[0:3]
  shifts = fmin[0:*]
  WIDGET_CONTROL,SET_VALUE='Image '+string(jj,format="(I2)")+' alignment parameters ['+string(fmin[0])+','+string(fmin[1])+']', widget.status
end


;############################################
PRO get_mask_corners_event, ev
;############################################

common data_common,data
common cornercount,ccnt,num

  WIDGET_CONTROL, ev.TOP, GET_UVALUE=eventval

  if (ev.type eq 0) then begin 
    if (data.curcorner eq 1) then begin
      data.click[ccnt,0,0] = ev.X
      data.click[ccnt,0,1] = ev.Y
      plots,[0,data.dim-1],[data.click[ccnt,0,1],data.click[ccnt,0,1]],color=200,/dev
      plots,[data.click[ccnt,0,0],data.click[ccnt,0,0]],[0,data.dim-1],color=200,/dev
      data.curcorner=2
    end else begin
      data.click[ccnt,1,0] = ev.X
      data.click[ccnt,1,1] = ev.Y
      plots,[0,data.dim-1],[data.click[ccnt,1,1],data.click[ccnt,1,1]],color=200,/dev
      plots,[data.click[ccnt,1,0],data.click[ccnt,1,0]],[0,data.dim-1],color=200,/dev
      data.curcorner=1
      ccnt += 1
      if (ccnt eq num) then WIDGET_CONTROL, ev.TOP, /DESTROY
    endelse
  endif

END


;############################################
pro get_mask_corners,dummy
;############################################

common data_common,data
common cornercount,ccnt,num

  IF(XRegistered("get_mask_corners") NE 0) THEN RETURN
;
  base= WIDGET_BASE(TITLE='click on lower left and upper right mask corners for '+string(num,format="(I4)")+' masks', /COLUMN)

  draw = WIDGET_DRAW(base, XSIZE=data.dim/2, YSIZE=data.dim/2, /BUTTON_EVENTS)

  ccnt = 0
  ; Realize the widget hierarchy.
  WIDGET_CONTROL, base, /REALIZE

  ; Retrieve the widget ID of the draw widget. Note that the widget
  ; hierarchy must be realized before you can retrieve this value.
  WIDGET_CONTROL, draw, GET_VALUE=drawID

  ; Make the draw widget the current IDL drawable area.
  WSET, drawID

  im = data.images
  tvscl,rebin(im,data.dim/2,data.dim/2)
  ; Call XMANAGER to manage the widgets.
  XMANAGER, 'get_mask_corners', base
end


;############################################
PRO align_series_ev, event
;############################################
;
common fd_common,fd
common widget_common,widget
common data_common,data
common results_common, edm, ws, wsb
common cornercount,ccnt,num
;
formats = [' ',' ','(I2.2)','(I3.3)','(I4.4)','(I5.5)','(I6.6)','(I7.7)','(I8.8)']
;
WIDGET_CONTROL, event.id, GET_UVALUE = eventval         ;find the user value
IF N_ELEMENTS(eventval) EQ 0 THEN RETURN
CASE eventval OF
 'NUMMASK'   : begin 
             WIDGET_CONTROL, GET_VALUE=val, widget.nummask
             data.nummask = fix(val(0))
             WIDGET_CONTROL, SET_VALUE=string(val,format='(I3)'), widget.nummask
           endcase
 'IMNM'  : begin 
             WIDGET_CONTROL, GET_VALUE=val, widget.imnm
             if (fix(val[0]) ge fd.num) then print,'enter from 0-',fd.num-1 else $
		     data.imnm = fix(val(0))
             WIDGET_CONTROL, SET_VALUE=string(val,format='(I2)'), widget.imnm
	     ; compute the filter range and display in text widget
	     openr,10,fd.path+fd.prefix+'.decon'
	     q=assoc(10,fltarr(data.dim,data.dim),12)
	     im=q[data.imnm]
	     close,10
	     data.images[0:*,0:*]=im
	     tvscl,rebin(im,data.dim,data.dim)
	     imf = mdg_filter_2d(im,'LoG',sigma=data.sigma)
	     range = '['+string(min(imf))+' - '+string(max(imf))+' ]'
	     WIDGET_CONTROL, SET_VALUE=range, widget.range
           endcase
 'SIGMA' : begin 
             WIDGET_CONTROL, GET_VALUE=val, widget.sigma
             data.sigma = float(val(0))
             WIDGET_CONTROL, SET_VALUE=string(val,format='(F6.2)'), widget.sigma
	     ; compute the filter range and display in text widget
	     im = data.images[0:*,0:*]
	     imf = mdg_filter_2d(im,'LoG',sigma=data.sigma)
	     range = '['+string(min(imf))+' - '+string(max(imf))+' ]'
	     WIDGET_CONTROL, SET_VALUE=range, widget.range
           endcase
 'THRS'  : begin 
             WIDGET_CONTROL, GET_VALUE=val, widget.thrs
             data.thrs = float(val(0))
             WIDGET_CONTROL, SET_VALUE=string(val,format='(F6.2)'), widget.thrs
           endcase
 'SHX'  : begin 
             WIDGET_CONTROL, GET_VALUE=val, widget.shx
             data.shx = fix(val(0))
             WIDGET_CONTROL, SET_VALUE=string(val,format='(I5)'), widget.shx
           endcase
 'SHY'  : begin 
             WIDGET_CONTROL, GET_VALUE=val, widget.shy
             data.shy = fix(val(0))
             WIDGET_CONTROL, SET_VALUE=string(val,format='(I5)'), widget.shy
           endcase
 'NPARAMS' : begin
            WIDGET_CONTROL, GET_VALUE=val,widget.nparams
            if (val eq 0) then data.nparams=4 else data.nparams=7
            endcase
 'METHOD': begin
	   WIDGET_CONTROL,GET_VALUE=val,widget.method
	   if (val eq 0) then data.method=0 else data.method=1
	   endcase
 'LOADFILES' : begin 
	 ; load all images of the series
	     openr,10,fd.path+fd.prefix+'.decon'
	     q=assoc(10,fltarr(data.dim,data.dim),12)
	     images=q[0]
	     close,10
	     sa = size(images,/dimensions)
 	     data.images[0:sa[0]-1,0:sa[1]-1] = float(images[0:*,0:*])
	     tvscl,rebin(images,data.dim,data.dim)
	     images = 0
	 ; define the intensity range for the DoG filtered image 1
	     im = data.images[0:*,0:*]
	     data.imnm=fix(0)
	     imf = mdg_filter_2d(im,'LoG',sigma=data.sigma)
	     range = '['+string(min(imf))+' - '+string(max(imf))+' ]'
	     WIDGET_CONTROL, SET_VALUE=range, widget.range
	 ; at the end of this option, activate the DoG button
	     WIDGET_CONTROL, SENSITIVE=1, widget.goDoGGT
	     WIDGET_CONTROL, SENSITIVE=1, widget.goDoGLT
	     WIDGET_CONTROL, SENSITIVE=1, widget.mkmask
	     WIDGET_CONTROL, SENSITIVE=1, widget.doseries
             endcase
 'DOGDISPLAYGT' : begin 
	     im = data.images[0:*,0:*]
	     imf = mdg_filter_2d(im,'LoG',sigma=data.sigma) * data.filter
	     data.gtlt = 1
	     tvscl, rebin(imf GT data.thrs,data.dim,data.dim)
           endcase
 'DOGDISPLAYLT' : begin 
	     im = data.images[0:*,0:*]
	     imf = mdg_filter_2d(im,'LoG',sigma=data.sigma)* data.filter
	     data.gtlt = -1
	     tvscl, rebin(imf LT data.thrs,data.dim,data.dim)
           endcase
 'MKMASK' : begin
	     num = data.nummask
	     get_mask_corners

	     ; create the mask 
	     wexp = 200.D0
	     q = fltarr(data.dim,data.dim)
	     for ii=0,data.nummask-1 do q[2*data.click[ii,0,0]:2*data.click[ii,1,0],2*data.click[ii,0,1]:2*data.click[ii,1,1]] = 1.0

	     ; define the coordinate arrays used to compute the Gaussian
	     x = dindgen(data.dim)-float(data.dim/2)
	     y = dindgen(data.dim)-float(data.dim/2)
	     
	     ; compute the 2D Gaussian and its integral (to normalize afterwards)
	     !EXCEPT=0
	     s = exp(-x^2/wexp) # exp(-y^2/wexp)
	     norm = total(s)
	     
	     ; FFT both arrays
	     fq = fft(q,-1,/double)
	     fs = fft(s,-1,/double)
	     
	     ; multiply both  (to compute their convolution) and correct amplitude; then shift to center
	     data.filter = float(shift(float(fft(fq*fs,1,/double)),data.dim/2,data.dim/2) * data.dim * data.dim / norm)
	     !EXCEPT=1
	     wset,data.drawID
	     im = data.images[0:*,0:*]
	     tvscl,rebin(data.filter * im,data.dim,data.dim)

	     ; once the filter is defined, we can proceed with the alignment
	     WIDGET_CONTROL, SENSITIVE=1, widget.dotierecon
	   endcase
'DOSERIES' : begin 
	       if (data.nparams eq 4) then begin
         shfts = [0.0,0.0,0.0,1.0]
	       shifts = [0.0,0.0,0.0,1.0]
         endif else begin
         shfts = [0.0,0.0,0.0,1.0,1.0,0.0,0.0]
         shifts = [0.0,0.0,0.0,1.0,1.0,0.0,0.0]
         endelse
	       ; reset the output file
	       openw,1,fd.path+fd.prefix+'.align'
	       close,1
	       for jj=1,data.nota-1 do begin
		     ; reset for the second series
		     if (jj eq (fd.num+1)/2) then begin
           if (data.nparams eq 4) then shfts = [0.0,0.0,0.0,1.0] $
             else shfts = [0.0,0.0,0.0,1.0,1.0,0.0,0.0]
         endif
		     sid = 0
	       amdog,'series',seriesid=sid,imageid=jj,shifts=shifts
		     if (data.nparams eq 4) then begin
         shfts[0:2] += shifts[0:2]
		     shfts[3] *= shifts[3]
       endif else begin
         shfts[0:2] += shifts[0:2]
         shfts[3:4] *= shifts[3:4]
         shfts[5:6] += shifts[5:6]
       endelse
	       print,jj,'final shift=',shfts
		     openw,1,fd.path+fd.prefix+'.align',/append
		     printf,1,shfts
	       close,1
         endfor
	 ; store the stack in the workfolder
         WIDGET_CONTROL,SET_VALUE='Alignment parameters stored in '+fd.path+fd.prefix+'.align', widget.status
	       WIDGET_CONTROL, SENSITIVE=1, widget.goalign
         widget_control,sensitive=1,widget.manalign
           endcase

'DOTIERECON' : begin
	; here we take the images in threes and try to reconstruct the phase for each triplet, using TIE

	; first, set the simulation parameters
	EE = 200000.0
	defstep=4.0
	delta = 30.18/1024.0
	params = simulationparameters(EE=EE,defstep=defstep,delta=delta)
	pimage=params[0]
	pscope=params[1]
	pphyscon=params[2]
	ptie=params[3]
	print,'ptie array dimension = ',(*pimage).dim
	; loop over all triplets
	phase = dblarr(data.nota,data.dim/2,data.dim/2)
	for i=1,data.nota-2 do begin
                WIDGET_CONTROL,SET_VALUE='Reconstructing phase for image triplet '+string(i,format="(I3)"), widget.status
		q=reform(data.images[i,256:256+data.dim/2-1,256:256+data.dim/2-1],data.dim/2,data.dim/2) 
		q = q/max(q)
		r=reform(data.images[i-1,256:256+data.dim/2-1,256:256+data.dim/2-1],data.dim/2,data.dim/2) 
		r = r/max(r)
		s=reform(data.images[i+1,256:256+data.dim/2-1,256:256+data.dim/2-1],data.dim/2,data.dim/2) 
		s = s/max(s)
		t = [total(q),total(r),total(s)]
		(*ptie).ima_if= q
		q = (s*t[0]/t[2] - r*t[0]/t[1])
		(*ptie).ima_dIdz= q - total(q)/float(512)^2
		;tvscl,(*ptie).ima_dIdz,256,256
		print,'total dIdz = ',total((*ptie).ima_dIdz)
		; next, do the reconstruction using the inverse gradient approach
		if (i eq 1) then mdgsolvetie,'gradient',ptie,/symmetrize else mdgsolvetie,'gradient',ptie,/skip,/symmetrize
		phase[i,0:*,0:*] = (*ptie).phase[0:*,0:*]
		tvscl,phase[i,*,*],256,256
                WIDGET_CONTROL,SET_VALUE='Reconstructed phase for image triplet '+string(i,format="(I3)"), widget.status
	endfor
	;clear up all the pointers
	ptr_free, pimage
	ptr_free, pscope
	ptr_free, pphyscon
	ptr_free, ptie
	stop
	endcase
'MANALIGN' : begin
  	widget_control,set_value='Starting Manual Alignment..',widget.status
		;read the aligned data set
		openr,1,fd.path+fd.prefix+'.aligned'
		q=assoc(1,fltarr(data.dim,data.dim),12)
		imstack=fltarr(data.dim,data.dim,3)
		msk = q[data.nota]
		for ii=0,data.nota-1 do imstack[*,*,ii] =q[ii]*msk
		close,1
		;Calling Manual alignment routine
		man_param=fltarr(4,3)
		res = man_align(imstack=imstack,params=man_param)
		if (res eq -1) then widget_control,set_value='Manual alignment Aborted..',widget.status else begin
			openw,1,fd.path+fd.prefix+'_manual.align'
			printf,1,man_param
			close,1
			widget_control,set_value='Manual alignment paramters stored in '+fd.path+fd.prefix+'_manual.align',widget.status
			ans=''
			read,prompt='Do you want to apply the manual alignment (y/n)?',ans
			if (ans eq 'y') then begin
        widget_control,set_value='Appending '+fd.path+fd.prefix+'.align file..',widget.status
        fin_param = fltarr(fd.num-1,4)
        man_param2 = man_param
        man_param = fltarr(fd.num-1,4)
        man_param[0,*] = man_param2[*,0]
        man_param[1,*] = man_param2[*,2]
        a=0.0 & b=0.0 & c=0.0 & d=0.0
        openr,1,fd.path+fd.prefix+'.align'
        for i=0,fd.num-2 do begin
	        readf,1,a,b,c,d
	        fin_param[i,0:2] = [a+man_param[i,0],b+man_param[i,1],c+man_param[i,2]]
          fin_param[i,3] = d*man_param[i,3]
        endfor
        close,1
        fin_param=transpose(fin_param)
        print,fin_param
        openw,1,fd.path+fd.prefix+'.align'
        printf,1,fin_param
        close,1
        apply_align,fd
      endif
    endelse
      endcase
'GOALIGN' : begin
	       ;apply_align,fd
	       res=file_info(fd.path+fd.prefix+'.align')
	    if (res.exists eq 0) then begin
		    widget_control,set_value='Cannot find *.align file.Please run DoSeries',widget.status
		    break
	    endif
      ;if (Data.nparams eq 4) then begin
	;    apply_align,fd
    ;endif else begin
	    ;read the params
	    openr,1,fd.path+fd.prefix+'.align'
	    if (data.nparams eq 4) then fparams=fltarr(4,fd.num-1) else fparams=fltarr(7,fd.num-1)
	    readf,1,fparams
	    close,1
	    ;open the decon files
	    openr,10,fd.path+fd.prefix+'.decon'
	    q10=assoc(10,fltarr(fd.dim,fd.dim),12)
	    openw,15,fd.path+fd.prefix+'.aligned'
	    q15=assoc(15,intarr(3))
	    q15=[fd.dim,fd.dim,fd.num]
	    q15=assoc(15,fltarr(fd.dim,fd.dim),12)
	    ;store infocus image as is
	    print,'Storing infocus image as is..'
	    q15[0]=q10[0]
	    print,'Allocating mask'
	    mask = replicate(1.0,data.dim,data.dim)
	   line = findgen(data.dim)-float(data.dim/2)
	    bintx=replicate(1.0,data.dim)##line
      binty=replicate(1.0,data.dim)#line
	
	    for i=1,fd.num-1 do begin
		    print,format="('starting image ',I5,' ',$)",i
		    im=q10[i]
		    if (data.nparams eq 4) then begin
		    b=rot(im,fparams[2,i-1],fparams[3,i-1],missing=-1.0,cubic=-0.5,/pivot)
		    im = interp2d(b,bintx,binty,bintx+fparams[0,i-1],binty+fparams[1,i-1],cubic=-0.5,missing=-1.0)
		    endif else begin
		     b = rot(im,fparams[2,i-1],missing=-1.0,cubic=-0.5)
		     im = interp2d(b,bintx,binty,bintx*fparams[3,i-1]+binty*fparams[5,i-1]+fparams[0,i-1],$
 bintx*fparams[6,i-1]+binty*fparams[4,i-1]+fparams[1,i-1],missing=-1.0,cubic=-0.5)
 			endelse
		     s=where(im eq -1.0,count)
		    if (count gt 0) then mask[s]=0.0
        if (count gt 0) then im[s]=0.0
		    tvscl,im
		    q15[i]=im
	    endfor
	    close,10
	    print,'saving mask'
	    q15[fd.num]=mask
	    close,15
	    print,'Total number of pixels in image = ',long(data.dim)^2
	    t = total(mask)
	    print,'Total number of pixels in mask = ',long(t)
	    print,'Percentage of pixels available for reconstruction = ',t/float(data.dim)^2*100.0
    ;endelse

     WIDGET_CONTROL, SENSITIVE=1, widget.showseries
     widget_control, sensitive=1, widget.manalign

	endcase
'LOADALIGNED' : begin 
               WIDGET_CONTROL,SET_VALUE='disabled in this version', widget.status
	   endcase
'SHOWALIGNEDSERIES' : begin 
               ;WIDGET_CONTROL,SET_VALUE='to be rewritten ', widget.status
	       wset,data.drawID
		openr,1,fd.path+fd.prefix+'.aligned'
		q=assoc(1,fltarr(data.dim,data.dim),12)
		msk=q[data.nota]
		for xx=0,5 do begin
    ;  ;under to over
    ;  for inum=0,fd.num/2-1 do begin
    ;    tvscl,q[fd.num/2-inum]*msk
    ;    wait,0.2
    ;  endfor
    ;  tvscl,q[0]*msk
    ;  wait,0.2
    ;  for inum=0,fd.num/2-1 do begin
    ;    tvscl,q[(fd.num+1)/2+inum]*msk
    ;    wait,0.2
    ;  endfor
    ;  ;over to under
    ;  for inum=0,fd.num/2-1 do begin
    ;    tvscl,q[fd.num-1-inum]*msk
    ;    wait,0.2
    ;  endfor
    ;  tvscl,q[0]*msk
    ;  wait,0.2
    ;  for inum=0,fd.num/2-1 do begin
    ;    tvscl,q[1+inum]*msk
    ;    wait,0.2
    ;  endfor
		       tvscl,q[0]*msk
		       wait,0.2
		       tvscl,q[1]*msk
		       wait,0.2
		       tvscl,q[0]*msk
		       wait,0.2
		       tvscl,q[(fd.num+1)/2]*msk
		       wait,0.2
	       endfor
	       close,1
               WIDGET_CONTROL,SET_VALUE='Waiting for input.....   ', widget.status
	   endcase

 'QUIT'  : begin 
             WIDGET_CONTROL, widget.base, /DESTROY
           endcase
 else: MESSAGE, "Event User Value Not Found"
endcase
end


;############################################
pro align_series,fd
;############################################
;
common fd_common,filedata
common widget_common,widget
common data_common,data
common results_common, edm, ws, wsb
;
num=fd.num
dim = fd.dim
filedata = fd

widget = {widgetstruct, base:long(0), draw:long(0), status:long(0), nota:long(0), tinc:long(0), stan:long(0), imnm:long(0), sigma:long(0), range:long(0), $
           thrs:long(0), goalign:long(0), manalign:long(0), goloop:long(0), mkmask:long(0), goDoGGT:long(0), goDoGLT:long(0), gettpa:long(0), tpa:long(0), doseries:long(0),  $
	   edm:long(0), wshed:long(0), wshedb:long(0), carbides:long(0), nummask:long(0), showseries:long(0), shx:long(0), shy:long(0), dotierecon:long(0), nparams:long(0), method:long(0)}

data = {datastruct, max_nota:fix(25), dim:fix(dim), status:'waiting for input', nota:fix(4), tinc:float(6), stan:float(0), imnm:fix(0), sigma:float(1.0), $
	range:'This will contain the range', thrs:float(128), tpa:float(0), filenames:strarr(25), images:fltarr(dim,dim), aligned:fltarr(dim,dim), $
	click:intarr(25,2,2), filter:dblarr(dim,dim), curcorner:fix(1), drawID:long(0), fits:fltarr(25,4), gtlt:fix(1), grains:fltarr(dim,dim),  $
	dirname:strarr(4), nummask:fix(1), shx:fix(0), shy:fix(0), prevorone:fix(0), nparams:fix(4), method:fix(2) $
	}
;
fontstr='-adobe-new century schoolbook-bold-r-normal--14-100-100-100-p-87-iso8859-1'
fontstrlarge='-adobe-new century schoolbook-medium-r-normal--20-140-100-100-p-103-iso8859-1'
fontstrsmall='-adobe-new century schoolbook-medium-r-normal--14-100-100-100-p-82-iso8859-1'
;


data.nota = fd.num

data.filter = replicate(1.D0,data.dim,data.dim)

IF(XRegistered("align_series") NE 0) THEN RETURN
;
DEVICE, GET_SCREEN_SIZE=scr
ydim = fix(scr[1] * 0.9); set the viewing box to be 80% of the screen height
xdim = fix(ydim * 1.45); make it a rectangular widget with aspect ratio 1.45

widget.base= WIDGET_BASE(TITLE='Alignment of Lorentz through-focus series', /COLUMN, SCR_XSIZE=xdim, SCR_YSIZE=ydim)

; current filetype and image size
base2= WIDGET_BASE(widget.base,/ROW)
base3= WIDGET_BASE(base2,/COLUMN)
base4= WIDGET_BASE(base2,/COLUMN)

; create the input windows
infoblock= WIDGET_BASE(base3,/ROW,/FRAME)

descrip = WIDGET_BASE(infoblock,/COLUMN,/BASE_ALIGN_LEFT)
fields = WIDGET_BASE(infoblock,/COLUMN,/BASE_ALIGN_RIGHT)
; create the drawing area
widget.draw= WIDGET_DRAW(base4,COLOR_MODEL=2,RETAIN=2,X_SCROLL_SIZE=ydim*0.9,Y_SCROLL_SIZE=ydim*0.9,/FRAME,/SCROLL,XSIZE=dim, YSIZE=dim)
widget.status= WIDGET_TEXT(base4,VALUE=string(data.status,FORMAT='(A)'),XSIZE=20,UVALUE='STAT')
;
; button to start load files widget
buttons1 = WIDGET_BASE(fields,/ALIGN_CENTER,/ROW)
loadfiles = WIDGET_BUTTON(buttons1,/ALIGN_CENTER,VALUE='Load Files',UVALUE='LOADFILES',/FRAME,SENSITIVE=1)
loadfiles = WIDGET_BUTTON(buttons1,/ALIGN_CENTER,VALUE='Load Aligned Dataset',UVALUE='LOADALIGNED',/FRAME,SENSITIVE=1)
;
infoblock2= WIDGET_BASE(base3,/ROW,/FRAME)
descrip2 = WIDGET_BASE(infoblock2,/COLUMN,/BASE_ALIGN_LEFT)
fields2 = WIDGET_BASE(infoblock2,/COLUMN,/BASE_ALIGN_RIGHT)
; image number to consider for DoG display
item= WIDGET_LABEL(descrip2,VALUE='Image # for DoG parameters',font=fontstrsmall,YSIZE=30,XSIZE=300,/ALIGN_LEFT)
widget.imnm= WIDGET_TEXT(fields2,VALUE=string(data.imnm,FORMAT='(I2)'),XSIZE=8,UVALUE='IMNM',/EDITABLE)
; sigma
item= WIDGET_LABEL(descrip2,VALUE='Sigma',font=fontstrsmall,YSIZE=30)
widget.sigma= WIDGET_TEXT(fields2,VALUE=string(data.sigma,FORMAT='(F8.2)'),XSIZE=10,UVALUE='SIGMA',/EDITABLE)
; threshold parameter
widget.range= WIDGET_TEXT(descrip2,VALUE=string(data.range,FORMAT='(A)'),XSIZE=30,UVALUE='RANGE')
widget.thrs= WIDGET_TEXT(fields2,VALUE=string(data.thrs,FORMAT='(F8.4)'),XSIZE=8,UVALUE='THRS',/EDITABLE)
; button to display the DoG filtered image
buttons2 = WIDGET_BASE(fields2,/ALIGN_CENTER,/COLUMN)
widget.goDoGGT = WIDGET_BUTTON(buttons2,/ALIGN_CENTER,VALUE='DoG > THR',UVALUE='DOGDISPLAYGT',/FRAME,SENSITIVE=0)
widget.goDoGLT = WIDGET_BUTTON(buttons2,/ALIGN_CENTER,VALUE='DoG < THR',UVALUE='DOGDISPLAYLT',/FRAME,SENSITIVE=0)
nparams = widget_base(descrip2,/align_center,/column)
widget.nparams = CW_BGROUP(nparams,['4','7'],/EXCLUSIVE,COLUMN=2,SET_VALUE=0,UVALUE='NPARAMS')
method = widget_base(descrip2,/align_center,/column)
widget.method = CW_BGROUP(method,['DoG','M.I'],/EXCLUSIVE,COLUMN=2,SET_VALUE=1,UVALUE='METHOD')
;
;
infoblock3= WIDGET_BASE(base3,/ROW,/FRAME)
descrip3 = WIDGET_BASE(infoblock3,/COLUMN,/BASE_ALIGN_LEFT)
fields3 = WIDGET_BASE(infoblock3,/COLUMN,/BASE_ALIGN_RIGHT)
; button to set the Region of Interest for the rectangular Gaussian mask
item= WIDGET_LABEL(descrip3,VALUE='Define Rectangular Gaussian Mask',font=fontstrsmall,YSIZE=30,XSIZE=300,/ALIGN_LEFT)
buttons3 = WIDGET_BASE(fields3,/ALIGN_CENTER,/ROW)
widget.nummask= WIDGET_TEXT(fields3,VALUE=string(data.nummask,FORMAT='(I3)'),XSIZE=5,UVALUE='NUMMASK',/EDITABLE)
widget.mkmask = WIDGET_BUTTON(fields3,/ALIGN_CENTER,VALUE='Go',UVALUE='MKMASK',/FRAME,SENSITIVE=0)
;
infoblock3b= WIDGET_BASE(base3,/ROW,/FRAME)
descrip3b = WIDGET_BASE(infoblock3b,/COLUMN,/BASE_ALIGN_LEFT)
fields3b = WIDGET_BASE(infoblock3b,/COLUMN,/BASE_ALIGN_RIGHT)
; overall shift parameters to correct for a poorly positioned first image (will be included in all alignment steps...)
item= WIDGET_LABEL(descrip3b,VALUE='Overall image shift parameters ',font=fontstrsmall,YSIZE=30,XSIZE=300,/ALIGN_LEFT)
buttons3 = WIDGET_BASE(fields3b,/ALIGN_CENTER,/ROW)
widget.shx= WIDGET_TEXT(fields3b,VALUE=string(data.shx,FORMAT='(I5)'),XSIZE=8,UVALUE='SHX',/EDITABLE)
widget.shy= WIDGET_TEXT(fields3b,VALUE=string(data.shy,FORMAT='(I5)'),XSIZE=8,UVALUE='SHY',/EDITABLE)
;
infoblock4= WIDGET_BASE(base3,/ROW,/ALIGN_CENTER)
; buttons to start the alingment and loop through after alingment
widget.doseries= WIDGET_BUTTON(infoblock4,/ALIGN_CENTER,VALUE='DoSeries',UVALUE='DOSERIES',/FRAME,SENSITIVE=0)
widget.showseries = WIDGET_BUTTON(infoblock4,/ALIGN_CENTER,VALUE='ShowalignedSeries',UVALUE='SHOWALIGNEDSERIES',/FRAME,SENSITIVE=1)
widget.goalign= WIDGET_BUTTON(infoblock4,/ALIGN_CENTER,VALUE='Apply Align',UVALUE='GOALIGN',/FRAME,SENSITIVE=1)
widget.manalign = WIDGET_BUTTON(infoblock4,/ALIGN_CENTER,VALUE='Manually Align',UVALUE='MANALIGN',/FRAME,SENSITIVE=0)

;widget.goloop = WIDGET_BUTTON(infoblock4,/ALIGN_CENTER,VALUE='Loop',UVALUE='GOLOOP',/FRAME,SENSITIVE=0)
;
infoblock5= WIDGET_BASE(base3,/ROW,/ALIGN_CENTER)
; buttons to start the alingment and loop through after alingment
widget.dotierecon= WIDGET_BUTTON(infoblock5,/ALIGN_CENTER,VALUE='DoTIErecon',UVALUE='DOTIERECON',/FRAME,SENSITIVE=1)

;
;infoblock5= WIDGET_BASE(base3,/ROW,/FRAME)
;descrip5 = WIDGET_BASE(infoblock5,/COLUMN,/BASE_ALIGN_LEFT)
;fields5 = WIDGET_BASE(infoblock5,/COLUMN,/BASE_ALIGN_RIGHT)
; 
;item= WIDGET_LABEL(descrip5,VALUE='Taper Angle [deg.] ',font=fontstrsmall,YSIZE=30,XSIZE=270,/ALIGN_LEFT)
;widget.tpa= WIDGET_TEXT(fields5,VALUE=string(data.tpa,FORMAT='(F6.2)'),XSIZE=8,UVALUE='TPA')
;buttons4 = WIDGET_BASE(fields5,/ALIGN_CENTER,/ROW)
;widget.gettpa = WIDGET_BUTTON(fields5,/ALIGN_CENTER,VALUE='Get Taper Angle',UVALUE='GETTPA',/FRAME,SENSITIVE=0)
;
;infoblock6= WIDGET_BASE(base3,/ROW,/FRAME)
; buttons for starting the entire alignment series
;buttons6 = WIDGET_BASE(infoblock6,/ALIGN_CENTER,/ROW)
;readalign= WIDGET_BUTTON(buttons6,/ALIGN_CENTER,VALUE='Read Alignment',UVALUE='READALIGN',/FRAME,SENSITIVE=1)
;
buttons5 = WIDGET_BASE(base3,/ALIGN_CENTER,/ROW)
quit  = WIDGET_BUTTON(buttons5,/ALIGN_CENTER,VALUE='Quit',UVALUE='QUIT',/FRAME)
;
WIDGET_CONTROL, /REALIZE, widget.base
WIDGET_CONTROL, widget.draw, GET_VALUE=drawID
data.drawID = drawID

; Make the draw widget the current IDL drawable area.
WSET, data.drawID

XMANAGER,"align_series",widget.base,EVENT_HANDLER="align_series_ev",/NO_BLOCK

;
return
end

