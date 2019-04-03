function joint_histogram,images,verbose=verbose
;
; this function returns the joint histogram, also known as the feature space, 
; for two, three, or four images or stacks.
;
; arguments:
; images               : bytarr(number_of_images,#cols,#rows[,#slices])
;
; options:
; verbose              : spit out some information
;
; This code combines the individual 2-D, 3-D, and 4-D joint histogram codes, and
; employs Dave Rowenhorst's suggestion to use the reverse indices option to histogram.
;
; Written by MDG, 9/27/07
;

; first make sure images is a byte array with 2, 3, or 4 images or stacks
sz = size(images)
if (sz[sz[0]+1] ne 1) then begin
	print,'This function requires a byte array as input'
	return,-1
endif

if ((sz[1] lt 2) or (sz[1] gt 4)) then begin
	print,'This function is only implemented for 2, 3, or 4 images or stacks'
	return,-1
endif

dim = 256
d2 = dim/2

case sz[1] of
  2 : begin
	  fs = lonarr(dim,dim)
	  if (sz[0] eq 3) then begin
	  	  hista = histogram(reform(images[0,*,*]), min=0, max=255, binsize=1, reverse_in=ri)
		  imageB = reform(images[1,*,*])
		  if keyword_set(verbose) then print,'two images detected'
	  end else begin
	  	  hista = histogram(reform(images[0,*,*,*]), min=0, max=255, binsize=1, reverse_in=ri)
		  imageB = reform(images[1,*,*,*])
		  if keyword_set(verbose) then print,'two data stacks detected'
	  endelse
	  for i=0, dim-1 do begin
	        if (hista[i] gt 0) then begin
			q = ri[ ri[i]:ri[i+1]-1 ]
			h = histogram(imageB[q], min=0, max=255, binsize=1)
			fs[i, 0:*] = h[0:*]
	 	endif
	  endfor
      end

  3 : begin
	  fs = lonarr(dim,dim,dim)
	  if (sz[0] eq 3) then begin
	  	  hista = histogram(reform(images[0,*,*]), min=0, max=255, binsize=1, reverse_in=ria)
		  imageB = reform(images[1,*,*])
		  imageC = reform(images[2,*,*])
		  if keyword_set(verbose) then print,'three images detected'
	  end else begin
	  	  hista = histogram(reform(images[0,*,*,*]), min=0, max=255, binsize=1, reverse_in=ria)
		  imageB = reform(images[1,*,*,*])
		  imageC = reform(images[2,*,*,*])
		  if keyword_set(verbose) then print,'three data stacks detected'
	  endelse
	  for i=0,dim-1 do begin
		  if (hista[i] gt 0) then begin
			q  = ria[ ria[i]:ria[i+1]-1 ]
			histb = histogram(imageB[q], min=0, max=255, binsize=1, reverse_in=rib)
			Cq = imageC[q]
			for j=0,dim-1 do begin
		  		if (histb[j] gt 0) then begin
					r = rib[ rib[j]:rib[j+1]-1 ]
					h = histogram(Cq[r],min=0,max=255,binsize=1)
					fs[i,j,0:*] = h[0:*]
				endif
			endfor
		  endif
	  endfor
      end

  4 : begin
	  fs = intarr(d2,d2,d2,d2)    ; note the intarr instead of lonarr, which would be too large !!!
	  if (sz[0] eq 3) then begin
	  	  hista = histogram(reform(images[0,*,*]), min=0, max=255, binsize=2, reverse_in=ria)
		  imageB = reform(images[1,*,*])
		  imageC = reform(images[2,*,*])
		  imageD = reform(images[3,*,*])
		  if keyword_set(verbose) then print,'four images detected'
	  end else begin
	  	  hista = histogram(reform(images[0,*,*,*]), min=0, max=255, binsize=2, reverse_in=ria)
		  imageB = reform(images[1,*,*,*])
		  imageC = reform(images[2,*,*,*])
		  imageD = reform(images[3,*,*,*])
		  if keyword_set(verbose) then print,'four data stacks detected'
	  endelse
	  for i=0,d2-1 do begin
		  if (hista[i] gt 0) then begin
			q = ria[ ria[i]:ria[i+1]-1 ]
			histb = histogram(imageB[q], min=0, max=255, binsize=2, reverse_in=rib)
			Cq = imageC[q]
			Dq = imageD[q]
			for j=0,d2-1 do begin
		  		if (histb[j] gt 0) then begin
					r = rib[ rib[j]:rib[j+1]-1 ]
					histc = histogram(Cq[r], min=0, max=255, binsize=2, reverse_in=ric)
					Dr = Dq[r]
					for k=0,d2-1 do begin
		  				if (histc[k] gt 0) then begin
							s = ric[ ric[k]:ric[k+1]-1 ]
							h = histogram(Dr[s],min=0,max=255,binsize=2)
							fs[i,j,k,0:*] = h[0:*]
						endif
					endfor
				endif
			endfor
		  endif
	  endfor
      end

  else: begin
	  print,'option ',sz[0],' is not implemented '
	  return,-1
	end
endcase

; spit out some information
if keyword_set(verbose) then begin
	help,fs
	mfs = max(fs)
	q = array_indices(fs,where(fs eq mfs))
	print,'Maximum value in joint histogram = ',mfs
	print,'at location ',q
endif

; that's it...
return,fs
end