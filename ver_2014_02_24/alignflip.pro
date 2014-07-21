@joint_entropy
@joint_histogram
@mdg_filter_2d
@interp2d
;+
;NAME: alignflip.pro
;
;CALLING SEQUENCE: alignflip,unflip=unflip_ptr,flip=flip_ptr
;
;PURPOSE: This program implements a GUI version of the alignflip code to align the unflipped and flipped images
;	        and saves the alignment files in the flipped file location.
;
;INPUT: unflip_ptr - pointer with file information about unflipped images generated either using rd_fls.pro or rdstack.pro
;	      flip_ptr - pointer with file information about flipped images.
;
;REQUIRES: joint_entropy.pro, joint_histrogram.pro, mdg_filter_2d.pro, interp2d.pro
;
;VERSION: 2.0 - GUI version.
;
;AUTHOR: CD Phatak, ANL, Feb. 24th, 2014.
;-
;
;############################################
FUNCTION getmin_fl,param
;############################################
compile_opt idl2

common data_common, data

line = findgen(data.dim)-float(data.dim/2)
bintx=replicate(1.0,data.dim)##line
binty=replicate(1.0,data.dim)#line

;check params
if (data.nparams eq 4) then begin
	c = rot(data.fl_img_al,param[2],param[3],missing=0.0,cubic=-0.5)
	c = interp2d(c,bintx,binty,bintx+param[0],binty+param[1],missing=0.0,cubic=-0.5)
endif else begin
	c = rot(data.fl_img_al,param[2],missing=0.0,cubic=-0.5)
	c = interp2d(c,bintx,binty,bintx*param[3]+binty*param[5]+param[0],bintx*param[6]+binty*param[4]+param[1],cubic=-0.5,missing=0.0)
endelse

;cc = smooth(c,3)
;aa = smooth(a,3)
a = data.un_img

;check method
if (data.method eq 0) then begin
	c = mdg_Filter_2d(c,'LoG',sigma=data.fl_sigma)
	if (data.gtlt eq 1) then c=float(c GT data.fl_thrs)*data.filter else c=float(c LT data.fl_thrs)*data.filter
	a = mdg_Filter_2d(a,'LoG',sigma=data.un_sigma)
	if (data.gtlt eq 1) then a=float(a GT data.un_thrs)*data.filter else a=float(a LT data.un_thrs)*data.filter
	;SSD Norm to be minimized
	cm = -total((c-a)^2)
endif else begin
	imgs=bytarr(2,data.dim,data.dim)
	imgs[0,*,*] = bytscl(a*data.filter)
	imgs[1,*,*] = bytscl(c*data.filter)
	res = joint_entropy(imgs,mutual=cm,/nozero)
endelse
tvscl,rebin(0.5*(a+c)*data.filter,data.dim/2,data.dim/2)
return,-cm
end

;############################################
PRO get_mask_corners_event, ev
;############################################

common data_common,data
common cornercount,ccnt,crnum

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
      if (ccnt eq crnum) then WIDGET_CONTROL, ev.TOP, /DESTROY
    endelse
  endif

END


;############################################
pro get_mask_corners,dummy
;############################################

common data_common,data
common cornercount,ccnt,crnum

  IF(XRegistered("get_mask_corners") NE 0) THEN RETURN
;
  base= WIDGET_BASE(TITLE='click on lower left and upper right mask corners for '+string(crnum,format="(I4)")+' masks', /COLUMN)

  draw = WIDGET_DRAW(base, XSIZE=data.dim/4, YSIZE=data.dim/4, /BUTTON_EVENTS)

  ccnt = 0
  ; Realize the widget hierarchy.
  WIDGET_CONTROL, base, /REALIZE

  ; Retrieve the widget ID of the draw widget. Note that the widget
  ; hierarchy must be realized before you can retrieve this value.
  WIDGET_CONTROL, draw, GET_VALUE=drawID

  ; Make the draw widget the current IDL drawable area.
  WSET, drawID
  ;setup the flipped image as per the settings
  if (data.revset eq 1) then begin
    imb=reverse(rot(data.fl_img,data.rotang,cubic=-0.5,missing=0.0),2)
    mskb=reverse(rot(data.flmask,data.rotang,missing=0.0),2)
  endif else begin
    imb = rot(data.fl_img,data.rotang,cubic=-0.5,missing=0.0)
    mskb = rot(data.flmask,data.rotang,missing=0.0)
  endelse
  imb = interpolate(imb,findgen(data.dim)-data.shx,findgen(data.dim)-data.shy,missing=0.0,cubic=-0.5,/grid)
  mskb = interpolate(mskb,findgen(data.dim)-data.shx,findgen(data.dim)-data.shy,missing=0.0,/grid)
  im = (data.un_img+imb)*data.unmask*mskb
  tvscl,rebin(im,data.dim/4,data.dim/4)
  ; Call XMANAGER to manage the widgets.
  XMANAGER, 'get_mask_corners', base
end

;---------------------------------------------------------------------------------------------------------------------------------------------
;Program that handles the main operations....
pro alignflip_ev, event
;
common fd_data, fid
common widget_common, widget
common data_common, data
common cornercount,ccnt,crn_num

WIDGET_CONTROL, event.id, GET_UVALUE= eventval ; Get the user value that generated the event...
if (N_ELEMENTS(eventval) eq 0) then return
case eventval of

 'NUMMASK'   : begin 
               WIDGET_CONTROL, GET_VALUE=val, widget.nummask
               data.nummask = fix(val(0))
               WIDGET_CONTROL, SET_VALUE=string(val,format='(I3)'), widget.nummask
               endcase
  'IMNM'     : begin
               WIDGET_CONTROL,GET_VALUE=val, widget.imnm
               if (val eq 0) then begin   ;This is unflip image
                 wset,data.drawid
                 tvscl,rebin(data.un_img,widget.wsi,widget.wsi)
                 WIDGET_CONTROL,SET_VALUE=string(data.un_sigma,format='(F6.2)'),widget.sigma
                 WIDGET_CONTROL,SET_VALUE=string(data.un_thrs,format='(F6.2)'),widget.thrs
                 im=data.un_img
                 im=mdg_filter_2d(im,'LoG',sigma=data.un_sigma)
                 range='['+string(min(im))+'-'+string(max(im))+']'
                 WIDGET_CONTROL,SET_VALUE=range,widget.range
               endif else begin     ;Then it is flip image
                 wset,data.drawid
                 tvscl,rebin(data.fl_img,widget.wsi,widget.wsi)
                 WIDGET_CONTROL,SET_VALUE=string(data.fl_sigma,format='(F6.2)'),widget.sigma
                 WIDGET_CONTROL,SET_VALUE=string(data.fl_thrs,format='(F6.2)'),widget.thrs
                 im=data.fl_img
                 im=mdg_filter_2d(im,'LoG',sigma=data.fl_sigma)
                 range='['+string(min(im))+'-'+string(max(im))+']'
                 WIDGET_CONTROL,SET_VALUE=range,widget.range
               endelse
               endcase
  'SIGMA'   :  begin
               WIDGET_CONTROL,GET_VALUE=val,widget.sigma
               ;Check which sigma is being updated
               WIDGET_CONTROL,GET_VALUE=val2,widget.imnm
               if (val2 eq 0) then data.un_sigma=float(val[0]) else data.fl_sigma=float(val[0])
               WIDGET_CONTROL,SET_VALUE=string(val,format='(F6.2)'),widget.sigma
               ;update the intensity range after filtering
               if (val2 eq 0) then im=data.un_img else im=data.fl_img
               imf = mdg_filter_2d(im,'LoG',sigma=float(val[0]))
               range='['+string(min(imf))+'-'+string(max(imf))+']'
               WIDGET_CONTROL,SET_VALUE=range,widget.range
               endcase
  'THRS'    :  begin
               WIDGET_CONTROL,GET_VALUE=val,widget.thrs
               ;Check which sigma is being updated
               WIDGET_CONTROL,GET_VALUE=val2,widget.imnm
               if (val2 eq 0) then data.un_thrs=float(val[0]) else data.fl_thrs=float(val[0])
               WIDGET_CONTROL,SET_VALUE=string(val,format='(F6.2)'),widget.thrs
               endcase
'DOGDISPLAYGT' : begin
                ;Check which image is being displayed
                WIDGET_CONTROL,GET_VALUE=val,widget.imnm
                if (val eq 0) then begin  ;Unflip image
                  im=data.un_img
                  imf = mdg_filter_2d(im,'LoG',sigma=data.un_sigma)*data.unmask*data.filter
                  data.gtlt = 1
                  wset,data.drawid
                  tvscl,rebin(imf GT data.un_thrs,widget.wsi,widget.wsi)
                endif else begin        ;Flip image
                  im=data.fl_img
                  imf = mdg_filter_2d(im,'LoG',sigma=data.fl_sigma)*data.flmask*data.filter
                  data.gtlt = 1
                  wset,data.drawid
                  tvscl,rebin(imf GT data.fl_thrs,widget.wsi,widget.wsi)
                endelse
                endcase
'DOGDISPLAYLT' : begin
                ;Check which image is being displayed
                WIDGET_CONTROL,GET_VALUE=val,widget.imnm
                if (val eq 0) then begin  ;Unflip image
                  im=data.un_img
                  imf = mdg_filter_2d(im,'LoG',sigma=data.un_sigma)*data.unmask*data.filter
                  data.gtlt = -1
                  wset,data.drawid
                  tvscl,rebin(imf LT data.un_thrs,widget.wsi,widget.wsi)
                endif else begin        ;Flip image
                  im=data.fl_img
                  imf = mdg_filter_2d(im,'LoG',sigma=data.fl_sigma)*data.flmask*data.filter
                  data.gtlt = -1
                  wset,data.drawid
                  tvscl,rebin(imf LT data.fl_thrs,widget.wsi,widget.wsi)
                endelse
                endcase
 'NPARAMS'    : begin
                WIDGET_CONTROL, GET_VALUE=val,widget.nparams
                if (val eq 0) then data.nparams=4 else data.nparams=7
                endcase
 'METHOD'     : begin
	              WIDGET_CONTROL,GET_VALUE=val,widget.method
	              if (val eq 0) then data.method=0 else data.method=1
	              endcase
 'SHX'        : begin
                WIDGET_CONTROL,GET_VALUE=val,widget.shx
                data.shx=fix(val[0])
                if (data.revset eq 1) then begin
                  imb=reverse(rot(data.fl_img,data.rotang,cubic=-0.5,missing=0.0),2)
                  mskb=reverse(rot(data.flmask,data.rotang,missing=0.0),2)
                endif else begin
                  imb = rot(data.fl_img,data.rotang,cubic=-0.5,missing=0.0)
                  mskb = rot(data.flmask,data.rotang,missing=0.0)
                endelse
                imb = interpolate(imb,findgen(data.dim)-data.shx,findgen(data.dim)-data.shy,missing=0.0,cubic=-0.5,/grid)
                mskb = interpolate(mskb,findgen(data.dim)-data.shx,findgen(data.dim)-data.shy,missing=0.0,/grid)
                wset,data.drawid
                tvscl,rebin(data.un_img+imb,widget.wsi,widget.wsi)
                WIDGET_CONTROL,SET_VALUE=string(val,format='(I5)'),widget.shx
                endcase
 'SHY'        : begin
                WIDGET_CONTROL,GET_VALUE=val,widget.shy
                data.shy=fix(val[0])
                if (data.revset eq 1) then begin
                  imb=reverse(rot(data.fl_img,data.rotang,cubic=-0.5,missing=0.0),2)
                  mskb=reverse(rot(data.flmask,data.rotang,missing=0.0),2)
                endif else begin
                  imb = rot(data.fl_img,data.rotang,cubic=-0.5,missing=0.0)
                  mskb = rot(data.flmask,data.rotang,missing=0.0)
                endelse
                imb = interpolate(imb,findgen(data.dim)-data.shx,findgen(data.dim)-data.shy,missing=0.0,cubic=-0.5,/grid)
                mskb = interpolate(mskb,findgen(data.dim)-data.shx,findgen(data.dim)-data.shy,missing=0.0,/grid)
                wset,data.drawid
                tvscl,rebin(data.un_img+imb,widget.wsi,widget.wsi)
                WIDGET_CONTROL,SET_VALUE=string(val,format='(I5)'),widget.shy
                endcase
 'ROTANG'        : begin
                WIDGET_CONTROL,GET_VALUE=val,widget.rotang
                data.rotang=fix(val[0])
                if (data.revset eq 1) then begin
                  imb=reverse(rot(data.fl_img,data.rotang,cubic=-0.5,missing=0.0),2)
                  mskb=reverse(rot(data.flmask,data.rotang,missing=0.0),2)
                endif else begin
                  imb = rot(data.fl_img,data.rotang,cubic=-0.5,missing=0.0)
                  mskb = rot(data.flmask,data.rotang,missing=0.0)
                endelse
                imb = interpolate(imb,findgen(data.dim)-data.shx,findgen(data.dim)-data.shy,missing=0.0,cubic=-0.5,/grid)
                mskb = interpolate(mskb,findgen(data.dim)-data.shx,findgen(data.dim)-data.shy,missing=0.0,/grid)
                wset,data.drawid
                tvscl,rebin(data.un_img+imb,widget.wsi,widget.wsi)
                WIDGET_CONTROL,SET_VALUE=string(val,format='(I5)'),widget.rotang
                endcase
 'REVSET'        : begin
                WIDGET_CONTROL,GET_VALUE=val,widget.revset
                if (val eq 0) then data.revset=0 else data.revset=1
                if (data.revset eq 1) then begin
                  imb=reverse(rot(data.fl_img,data.rotang,cubic=-0.5,missing=0.0),2)
                  mskb=reverse(rot(data.flmask,data.rotang,missing=0.0),2)
                endif else begin
                  imb = rot(data.fl_img,data.rotang,cubic=-0.5,missing=0.0)
                  mskb = rot(data.flmask,data.rotang,missing=0.0)
                endelse
                imb = interpolate(imb,findgen(data.dim)-data.shx,findgen(data.dim)-data.shy,missing=0.0,cubic=-0.5,/grid)
                mskb = interpolate(mskb,findgen(data.dim)-data.shx,findgen(data.dim)-data.shy,missing=0.0,/grid)
                wset,data.drawid
                tvscl,rebin(data.un_img+imb,widget.wsi,widget.wsi)
                endcase
 'SHOWSUM'        : begin
                if (data.revset eq 1) then begin
                  imb=reverse(rot(data.fl_img,data.rotang,cubic=-0.5,missing=0.0),2)
                  mskb=reverse(rot(data.flmask,data.rotang,missing=0.0),2)
                endif else begin
                  imb = rot(data.fl_img,data.rotang,cubic=-0.5,missing=0.0)
                  mskb = rot(data.flmask,data.rotang,missing=0.0)
                endelse
                imb = interpolate(imb,findgen(data.dim)-data.shx,findgen(data.dim)-data.shy,missing=0.0,cubic=-0.5,/grid)
                mskb = interpolate(mskb,findgen(data.dim)-data.shx,findgen(data.dim)-data.shy,missing=0.0,/grid)
                wset,data.drawid
                tvscl,rebin(data.un_img+imb,widget.wsi,widget.wsi)
                endcase
 'MKMASK'     : begin
                crn_num = data.nummask
	              get_mask_corners

	              ; create the mask 
	              wexp = 200.D0
	              q = fltarr(data.dim,data.dim)
	              for ii=0,data.nummask-1 do q[4*data.click[ii,0,0]:4*data.click[ii,1,0],4*data.click[ii,0,1]:4*data.click[ii,1,1]] = 1.0

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
                if (data.revset eq 1) then begin
                  imb=reverse(rot(data.fl_img,data.rotang,cubic=-0.5,missing=0.0),2)
                  mskb=reverse(rot(data.flmask,data.rotang,missing=0.0),2)
                endif else begin
                  imb = rot(data.fl_img,data.rotang,cubic=-0.5,missing=0.0)
                  mskb = rot(data.flmask,data.rotang,missing=0.0)
                endelse
                imb = interpolate(imb,findgen(data.dim)-data.shx,findgen(data.dim)-data.shy,missing=0.0,cubic=-0.5,/grid)
                mskb = interpolate(mskb,findgen(data.dim)-data.shx,findgen(data.dim)-data.shy,missing=0.0,/grid)
                data.filter *= (data.unmask*mskb)
	              im = (data.un_img+imb)*data.unmask*mskb*data.filter
	              tvscl,rebin(im,widget.wsi,widget.wsi)
	              ; once the filter is defined, we can proceed with the alignment
	              WIDGET_CONTROL, SENSITIVE=1, widget.doseries
	              endcase
  'DOALIGN'   : begin
              ;Now to begin the actual alignment process
                ;Unflipped image
                aa = smooth(data.un_img,3)
                aa = float(mdg_Filter_2d(aa,'LoG',sigma=data.un_sigma))
		if (data.gtlt eq 1) then aa = aa GT data.un_thrs else aa = aa LT data.un_thrs
		tvscl,rebin(aa,widget.wsi,widget.wsi)
		fa = fft(aa*data.filter,-1)
		
                ;Flipped image
                ;correct the orientation of the flipped image
                if (data.revset eq 1) then begin
                  imb=reverse(rot(data.fl_img,data.rotang,cubic=-0.5,missing=0.0),2)
                endif else begin
                  imb = rot(data.fl_img,data.rotang,cubic=-0.5,missing=0.0)
                endelse
                imb = interpolate(imb,findgen(data.dim)-data.shx,findgen(data.dim)-data.shy,missing=0.0,cubic=-0.5,/grid)
                data.fl_img_al = imb
                imb = smooth(imb,3)
                imb = float(mdg_filter_2d(imb,'LoG',sigma=data.fl_sigma))
                if (data.gtlt eq 1) then imb = imb GT data.fl_thrs else imb = imb LT data.fl_thrs
                tvscl,rebin(imb,widget.wsi,widget.wsi)
		fb = conj(fft(imb*data.filter,-1))
		ab = smooth(abs(fft(fa*fb,1)),5)
		tvscl,rebin(shift(alog(ab+0.01),data.dim/2,data.dim/2),widget.wsi,widget.wsi)
		;determine the max
		qmax=where(ab eq max(ab))
      		roughx = qmax[0] mod data.dim
        	roughy = qmax[0]/data.dim
        	if (roughx gt data.dim/2) then roughx = roughx - data.dim
        	if (roughy gt data.dim/2) then roughy = roughy - data.dim
        	plots,(data.dim/2+[roughx-3,roughx+3])/2,(data.dim/2+[roughy,roughy])/2,/dev,color=230
        	plots,(data.dim/2+[roughx,roughx])/2,(data.dim/2+[roughy-3,roughy+3])/2,/dev,color=230
		;start sub-pixel
		ftol=1.0e-4
		stat_updt = 'Rough shifts= '+string(-roughx)+string(-roughy)
		WIDGET_CONTROL,SET_VALUE=stat_updt,widget.status
		if (data.nparams eq 4) then begin
			pa = [-roughx,-roughy,0.0,1.0]
			scl = [5.0,5.0,5.0,0.05]
		endif else begin
			pa = [-roughx,-roughy,0.0,1.0,1.0,0.0,0.0]
			scl = [1.0,1.0,2.0,0.1,0.1,0.1,0.1]
		endelse
		iter=0
		itmax=3000
		fmin=amoeba(ftol,function_name='getmin_fl',function_value=fv,ncalls=iter,nmax=itmax,p0=pa,scale=scl)
		print,'Iterations=',iter
		print,'Fits=',fmin
		;resetting the file
		openw,1,data.alignfile
		close,1
		openw,2,data.alignfile,/append
		printf,2,data.nparams
		printf,2,[data.shx,data.shy,data.rotang,data.revset]
		printf,2,fmin
		close,2
		;update the params in data structure
		if (data.nparams eq 4) then data.params[0:3]=fmin else data.params[0:6]=fmin
		stat_updt = 'Alignment complete. Params saved in '+data.alignfile+'. Press Show Align to see the alignment.'
		;Apply the alignment and save the image into data structure
		line = findgen(data.dim)-float(data.dim/2)
		bintx = replicate(1.0,data.dim)##line
		binty = replicate(1.0,data.dim)#line
		if (data.nparams eq 4) then begin
			imb = rot(data.fl_img_al,data.params[2],data.params[3],cubic=-0.5,missing=0.0)
			imb = interp2d(imb,bintx,binty,bintx+data.params[0],binty+params[1],cubic=-0.5,missing=0.0)
		endif else begin
			imb = rot(data.fl_img_al,data.params[2],cubic=-0.5,missing=0.0)
			imb = interp2d(imb,bintx,binty,bintx*data.params[3]+binty*data.params[5]+data.params[0],$
						bintx*data.params[6]+binty*data.params[4]+data.params[1],missing=0.0,cubic=-0.5)
		endelse
		data.fl_img_al = imb
		WIDGET_CONTROL,SET_VALUE=stat_updt,widget.status
		WIDGET_CONTROL,SENSITIVE=1,widget.showseries
		WIDGET_CONTROL,SENSITIVE=1,widget.goalign
		endcase
  'SHOWALIGN' : begin
		wset,data.drawid
		for cnt = 0,10 do begin
			tvscl,rebin(data.un_img*data.filter,widget.wsi,widget.wsi)
			wait,0.1
			tvscl,rebin(data.fl_img_al*data.filter,widget.wsi,widget.wsi)
			wait,0.1
		endfor
		endcase
  'GOALIGN'   : begin
		stat_updt = 'Applying Alignment to Flipped Images...'
		print,data.params
		WIDGET_CONTROL,SET_VALUE=stat_updt,widget.status
		;Input data
		openr,10,fid[1].path+fid[1].prefix+'.aligned'
		numim = fid[1].num
		q10 = assoc(10,fltarr(data.dim,data.dim),12)
		;Output data
		openw,20,fid[1].path+fid[1].prefix+'_final.aligned'
		q20=assoc(20,intarr(3))
		q20[0]=[data.dim,data.dim,numim]
		q20=assoc(20,fltarr(data.dim,data.dim),12)
		;Mask from flipped alignment
		mask = replicate(1.0,data.dim,data.dim)*data.flmask
		line = findgen(data.dim)-float(data.dim/2)
		bintx = replicate(1.0,data.dim)##line
		binty = replicate(1.0,data.dim)#line
		for ii=0,numim-1 do begin
			print,'Image #',ii
			imb = q10[ii]
			;Correct the flip...
			 if (data.revset eq 1) then begin
                  		imb=reverse(rot(imb,data.rotang,cubic=-0.5,missing=-1.0),2)
                	 endif else begin
                  		imb = rot(imb,data.rotang,cubic=-0.5,missing=-1.0)
                	 endelse
                	 imb = interpolate(imb,findgen(data.dim)-data.shx,findgen(data.dim)-data.shy,missing=-1.0,cubic=-0.5,/grid)
			 ;Apply alignment
			 if (data.nparams eq 4) then begin
				imb = rot(imb,data.params[2],data.params[3],cubic=-0.5,missing=-1.0)
				imb = interp2d(imb,bintx,binty,bintx+data.params[0],binty+params[1],cubic=-0.5,missing=-1.0)
		         endif else begin
				imb = rot(imb,data.params[2],cubic=-0.5,missing=-1.0)
				imb = interp2d(imb,bintx,binty,bintx*data.params[3]+binty*data.params[5]+data.params[0],$
						bintx*data.params[6]+binty*data.params[4]+data.params[1],missing=-1.0,cubic=-0.5)
		 	 endelse
			 ;modify the mask
			 qq=where(imb eq -1.0)
			 if (qq[0] ne -1) then mask[qq]=0.0
			 ;store the image
			 q20[ii]=imb
			 ;Show the image
			 wset,data.drawid
			 tvscl,rebin(imb,widget.wsi,widget.wsi)
		endfor
		print,'----> Done..'
		;save mask
		q20[numim]=mask
		;close files
		close,10
		close,20
		endcase
  'QUIT'      : begin
                WIDGET_CONTROL,widget.base,/DESTROY
                endcase
  else: MESSAGE,"Event User Value Not Found.."
endcase
end

;---------------------------------------------------------------------------------------------------------------------------------------------
;Program that defines the GUI window and related structures....

pro alignflip,unflip=unflip,flip=flip
compile_opt idl2

common fd_data,fid
common widget_common, widget
common data_common, data

wsi=1024
dim=unflip.dim
fid=[unflip,flip]

widget = {widgetstruct, base:long(0), draw:long(0), status:long(0),rotang:long(0), stan:long(0), imnm:long(0), sigma:long(0), range:long(0), $
           thrs:long(0), goalign:long(0), manalign:long(0), goloop:long(0), mkmask:long(0), goDoGGT:long(0), goDoGLT:long(0), gettpa:long(0), tpa:long(0), doseries:long(0),  $
	   wsi:long(wsi), showsum:long(0), wshedb:long(0), revset:long(0), nummask:long(0), showseries:long(0), shx:long(0), shy:long(0), dotierecon:long(0), nparams:long(0), method:long(0)}

data = {datastruct, dim:fix(dim), status:'waiting for input', revset:fix(0), fl_sigma:float(1.0), rotang:fix(0), imnm:fix(0), un_sigma:float(1.0), $
	range:'This will contain the range', un_thrs:float(128), fl_thrs:float(128), alignfile:flip.path+flip.prefix+'_flipunflip.align', un_img:fltarr(dim,dim), fl_img:fltarr(dim,dim), $
	click:intarr(25,2,2), unmask:fltarr(dim,dim), flmask:fltarr(dim,dim),curcorner:fix(1), drawID:long(0), params:fltarr(7), gtlt:fix(1), fl_img_al:fltarr(dim,dim),  $
	dirname:strarr(4), nummask:fix(1), shx:fix(0), shy:fix(0), filter:replicate(1.0,dim,dim), nparams:fix(4), method:fix(2) $
	}

;Load the in-focus images from each set into the data structure
openr,1,unflip.path+unflip.prefix+'.aligned'
q=assoc(1,fltarr(dim,dim),12)
data.unmask=q[unflip.num]
data.un_img=q[0]*data.unmask
close,1
openr,1,flip.path+flip.prefix+'.aligned'
q=assoc(1,fltarr(dim,dim),12)
data.flmask=q[flip.num]
data.fl_img=q[0]*data.flmask
close,1

IF (XRegistered("alignflip") NE 0) then return

DEVICE, GET_SCREEN_SIZE=scr
ydim = fix(scr[1] * 0.9) < wsi; set the viewing box to be 80% of the screen height
xdim = fix(ydim * 1.4); make it a rectangular widget with aspect ratio 1.45

widget.base= WIDGET_BASE(TITLE='Alignment of Unflip/Flip images', /COLUMN, SCR_XSIZE=xdim, SCR_YSIZE=ydim)
base2 = WIDGET_BASE(widget.base,/ROW)
base3 = WIDGET_BASE(base2,/column)
base4 = WIDGET_BASE(base2,/column)
;info blocks

; create the drawing area
widget.draw= WIDGET_DRAW(base4,COLOR_MODEL=2,RETAIN=2,X_SCROLL_SIZE=ydim*0.9,Y_SCROLL_SIZE=ydim*0.9,/FRAME,/SCROLL,XSIZE=wsi, YSIZE=wsi)
widget.status= WIDGET_TEXT(base4,VALUE=string(" ",FORMAT='(A)'),XSIZE=20,UVALUE='STAT')

infoblock2= WIDGET_BASE(base3,/ROW,/FRAME)
descrip2 = WIDGET_BASE(infoblock2,/COLUMN,/BASE_ALIGN_LEFT)
fields2 = WIDGET_BASE(infoblock2,/COLUMN,/BASE_ALIGN_RIGHT)
; image number to consider for DoG display
item= WIDGET_LABEL(descrip2,VALUE='Image',font=fontstrsmall,YSIZE=30,XSIZE=200,/ALIGN_LEFT)
widget.imnm = CW_BGROUP(fields2,['Unflip','Flip'],/EXCLUSIVE,COLUMN=2,SET_VALUE=0,UVALUE='IMNM')
; sigma
item= WIDGET_LABEL(descrip2,VALUE='Sigma',font=fontstrsmall,YSIZE=30)
widget.sigma= WIDGET_TEXT(fields2,VALUE=string(data.un_sigma,FORMAT='(F8.2)'),XSIZE=10,UVALUE='SIGMA',/EDITABLE)
; threshold parameter
item=WIDGET_LABEL(descrip2,VALUE='Range',font=fontstrsmall,ysize=30)
widget.range= WIDGET_TEXT(fields2,VALUE=string(data.range,FORMAT='(A)'),XSIZE=30,UVALUE='RANGE')
item=WIDGET_LABEL(descrip2,VALUE='Threshold',font=fontstrsmall,ysize=30)
widget.thrs= WIDGET_TEXT(fields2,VALUE=string(data.un_thrs,FORMAT='(F8.4)'),XSIZE=8,UVALUE='THRS',/EDITABLE)
; button to display the DoG filtered image
buttons2 = WIDGET_BASE(fields2,/ALIGN_RIGHT,/ROW)
widget.goDoGGT = WIDGET_BUTTON(buttons2,/ALIGN_CENTER,VALUE='DoG > THR',UVALUE='DOGDISPLAYGT',/FRAME,SENSITIVE=1)
widget.goDoGLT = WIDGET_BUTTON(buttons2,/ALIGN_CENTER,VALUE='DoG < THR',UVALUE='DOGDISPLAYLT',/FRAME,SENSITIVE=1)
;
;
infoblock3b= WIDGET_BASE(base3,/ROW,/FRAME)
descrip3b = WIDGET_BASE(infoblock3b,/COLUMN,/BASE_ALIGN_LEFT)
fields3b = WIDGET_BASE(infoblock3b,/COLUMN,/BASE_ALIGN_RIGHT)
; overall shift parameters to correct for a poorly positioned image
item= WIDGET_LABEL(descrip3b,VALUE='Overall image shift parameters ',font=fontstrsmall,YSIZE=30,XSIZE=300,/ALIGN_LEFT)
buttons3 = WIDGET_BASE(fields3b,/ALIGN_CENTER,/ROW)
widget.shx= WIDGET_TEXT(fields3b,VALUE=string(data.shx,FORMAT='(I5)'),XSIZE=8,UVALUE='SHX',/EDITABLE)
widget.shy= WIDGET_TEXT(fields3b,VALUE=string(data.shy,FORMAT='(I5)'),XSIZE=8,UVALUE='SHY',/EDITABLE)
; overall rotate and reverse setting parameters
infoblock3c= WIDGET_BASE(base3,/ROW,/FRAME)
descrip3c = WIDGET_BASE(infoblock3c,/COLUMN,/BASE_ALIGN_LEFT)
fields3c = WIDGET_BASE(infoblock3c,/COLUMN,/BASE_ALIGN_RIGHT)
item= WIDGET_LABEL(descrip3c,VALUE='Overall image rotation angle',font=fontstrsmall,YSIZE=30,XSIZE=300,/ALIGN_LEFT)
buttons3 = WIDGET_BASE(fields3c,/ALIGN_CENTER,/ROW)
widget.rotang= WIDGET_TEXT(fields3c,VALUE=string(data.shx,FORMAT='(I5)'),XSIZE=8,UVALUE='ROTANG',/EDITABLE)
item = WIDGET_LABEL(descrip3c,VALUE='Reverse setting for images',font=fontstrsmall,YSIZE=30,XSIZE=300,/ALIGN_LEFT)
rev_set = WIDGET_BASE(fields3c,/align_left,/column)
widget.revset= CW_BGROUP(rev_set,['0','1'],/EXCLUSIVE,COLUMN=2,SET_VALUE=1,UVALUE='REVSET')
but = WIDGET_BASE(descrip3c,/align_center,/column)
widget.showsum = WIDGET_BUTTON(but,/ALIGN_CENTER,VALUE='Show Sum',UVALUE='SHOWSUM',/FRAME)
;
;Alignment settings
infoblock3d= WIDGET_BASE(base3,/ROW,/FRAME)
descrip3d = WIDGET_BASE(infoblock3d,/COLUMN,/BASE_ALIGN_LEFT)
fields3d = WIDGET_BASE(infoblock3d,/COLUMN,/BASE_ALIGN_RIGHT)
item=widget_label(descrip3d,VALUE='No. of parameters',font=fontstrsmall,YSIZE=30,XSIZE=300,/ALIGN_LEFT)
nparams = widget_base(fields3d,/align_center,/column)
widget.nparams = CW_BGROUP(nparams,['4','7'],/EXCLUSIVE,COLUMN=2,SET_VALUE=0,UVALUE='NPARAMS')
method = widget_base(fields3d,/align_center,/column)
item=widget_label(descrip3d,VALUE='Method for Alignement',font=fontstrsmall,YSIZE=30,XSIZE=300,/ALIGN_LEFT)
widget.method = CW_BGROUP(method,['DoG','M.I'],/EXCLUSIVE,COLUMN=2,SET_VALUE=1,UVALUE='METHOD')
;
infoblock3= WIDGET_BASE(base3,/ROW,/FRAME)
descrip3 = WIDGET_BASE(infoblock3,/COLUMN,/BASE_ALIGN_LEFT)
fields3 = WIDGET_BASE(infoblock3,/COLUMN,/BASE_ALIGN_RIGHT)
; button to set the Region of Interest for the rectangular Gaussian mask
item= WIDGET_LABEL(descrip3,VALUE='Define Rectangular Gaussian Mask',font=fontstrsmall,YSIZE=30,XSIZE=300,/ALIGN_LEFT)
buttons3 = WIDGET_BASE(fields3,/ALIGN_RIGHT,/ROW)
widget.nummask= WIDGET_TEXT(fields3,VALUE=string(data.nummask,FORMAT='(I3)'),XSIZE=5,UVALUE='NUMMASK',/EDITABLE)
widget.mkmask = WIDGET_BUTTON(fields3,/ALIGN_CENTER,VALUE='Go',UVALUE='MKMASK',/FRAME,SENSITIVE=1)
;
;Action Buttons
infoblock4= WIDGET_BASE(base3,/ROW,/ALIGN_CENTER)
; buttons to start the alingment and loop through after alingment
widget.doseries= WIDGET_BUTTON(infoblock4,/ALIGN_CENTER,VALUE='Do Align',UVALUE='DOALIGN',/FRAME,SENSITIVE=0)
widget.showseries = WIDGET_BUTTON(infoblock4,/ALIGN_CENTER,VALUE='Show Align',UVALUE='SHOWALIGN',/FRAME,SENSITIVE=0)
widget.goalign= WIDGET_BUTTON(infoblock4,/ALIGN_CENTER,VALUE='Apply Align',UVALUE='GOALIGN',/FRAME,SENSITIVE=0)
;
;End
buttons5 = WIDGET_BASE(base3,/ALIGN_CENTER,/ROW)
quit  = WIDGET_BUTTON(buttons5,/ALIGN_CENTER,VALUE='Quit',UVALUE='QUIT',/FRAME)
;
WIDGET_CONTROL, /REALIZE, widget.base
WIDGET_CONTROL, widget.draw, GET_VALUE=drawID
data.drawID = drawID
;
;check if alignment file exists and update numbers accordingly
fi=file_info(data.alignfile)
if (fi.exists eq 1) then begin
  openr,1,data.alignfile
  tmp=0
  readf,1,tmp
  if (tmp eq 4) then begin
    params=fltarr(4)
    WIDGET_CONTROL,SET_VALUE=0,widget.nparams
    data.nparams=4
  endif else begin
    params=fltarr(7)
    WIDGET_CONTROL,SET_VALUE=1,widget.nparams
    data.nparams=7
  endelse
  al_set=fltarr(4)
  readf,1,al_set
  readf,1,params
  close,1
  if (data.nparams eq 4) then data.params[0:3]=params else data.params[0:6]=params
  data.shx=al_set[0]
  data.shy=al_set[1]
  data.rotang=al_set[2]
  data.revset=al_set[3]
  WIDGET_CONTROL,SET_VALUE=string(al_set[0],format='(I5)'),widget.shx
  WIDGET_CONTROL,SET_VALUE=string(al_set[1],format='(I5)'),widget.shy
  WIDGET_CONTROL,SET_VALUE=string(al_set[2],format='(I5)'),widget.rotang
  WIDGET_CONTROL,SET_VALUE=fix(al_set[3]),widget.revset
  WIDGET_CONTROL,SENSITIVE=1,widget.showseries
  WIDGET_CONTROL,SENSITIVE=1,widget.goalign
  WIDGET_CONTROL,SET_VALUE='Alignment file exists. Pre-loaded the values. Press Show Align to check.',widget.status
  ;Apply the alignment parameters
  if (data.revset eq 1) then begin
	 imb=reverse(rot(data.fl_img,data.rotang,cubic=-0.5,missing=0.0),2)
  endif else begin
         imb = rot(data.fl_img,data.rotang,cubic=-0.5,missing=0.0)
  endelse
  imb = interpolate(imb,findgen(data.dim)-data.shx,findgen(data.dim)-data.shy,missing=0.0,cubic=-0.5,/grid)
  line = findgen(data.dim)-float(data.dim/2)
  bintx = replicate(1.0,data.dim)##line
  binty = replicate(1.0,data.dim)#line
  if (data.nparams eq 4) then begin
      imb = rot(imb,data.params[2],data.params[3],cubic=-0.5,missing=0.0)
      imb = interp2d(imb,bintx,binty,bintx+data.params[0],binty+params[1],cubic=-0.5,missing=0.0)
  endif else begin
      imb = rot(imb,data.params[2],cubic=-0.5,missing=0.0)
      imb = interp2d(imb,bintx,binty,bintx*data.params[3]+binty*data.params[5]+data.params[0],$
			bintx*data.params[6]+binty*data.params[4]+data.params[1],missing=0.0,cubic=-0.5)
  endelse
  data.fl_img_al = imb
endif
;
; Make the draw widget the current IDL drawable area.
WSET, data.drawID
tvscl,rebin(data.un_img,widget.wsi,widget.wsi)
;
XMANAGER,"alignflip",widget.base,EVENT_HANDLER="alignflip_ev",/NO_BLOCK
;
;
return
end

