function  mdg_filter_2D,image,ftype,sigma=sigma,maxsize=maxsize
;
;+
; NAME:
;     mdg_filter_2D
; PURPOSE:
;     Compute one of several discrete filters in 2D; return it as a normalized array for
;     convolution operations;  the size of the filter is three times the standard 
;     deviation (times two).  If maxsize is specified, then it is used instead of 
;     3*sigma.
;
; EXPLANATION:
;
; CALLING SEQUENCE:
;     filter= mdg_filter_2D(image,ftype,sigma[,maxsize=maxsize])
;
; INPUT PARAMETERS
;     image = 2D input image
;     ftype = one off the following ['Mean','Median','Gaussian','ConservativeSmoothing', 
;             'LoG','Unsharp']
;     sigma = standard deviation for Gaussian and related filters
;     maxsize (if present) = maximum array dimension for filter (in pixels)
;
; OUTPUT PARAMETERS
;     RESULT - filtered image as a floating point image
;
; NOTES:

; REVISION HISTORY
;     Written by MDG, 5/19/06
;-

; compute the size of the filter unless maxsize is defined (maxsize overrides sigma)
if arg_present(maxsize) then begin $
  dim = fix(maxsize)
end else begin $
  dim = round(6.0*sigma)
endelse
; make sure filter has odd size and is at least 3 pixels
if (dim mod 2 eq 0) then dim=dim+1
if (dim lt 3) then dim=3
d2 = (dim-1)/2

; some useful auxiliary variables
kernel = fltarr(dim,dim)
line = findgen(dim)-float(dim/2)
x = line#replicate(1.0,dim)
y = rotate(x,3)
r = x^2+y^2

fimage = float(image)
si = size(image)

case (ftype) of
 'Mean'    :  begin $
                kernel = replicate(1.0/float(dim)^2,dim,dim)
                newimage = convol(fimage,kernel,/edge_truncate)
              end
 'Median'  :  begin $
                newimage = median(fimage,dim,/even)
              end
 'Gaussian':  begin $
                kernel = exp(-r*0.5/sigma^2)/(2.0*!pi*sigma^2)
                newimage = convol(fimage,kernel,/edge_truncate)
              end
 'ConservativeSmoothing' : begin $
                newimage = fimage
                for i=d2,si[1]-1-d2 do begin $
                  for j=d2,si[2]-1-d2 do begin $
                    w = reform(fimage[i-d2:i+d2,j-d2:j+d2],dim,dim)
                    v = w[d2,d2]
                    w[d2,d2]=(total(w)-v)/(float(dim)^2-1.0)
                    miw = min(w)
                    maw = max(w)
                    if (v lt miw) then newimage[i,j] = miw
                    if (v gt maw) then newimage[i,j] = maw
                  endfor
                endfor
              end
 'LoG'     :  begin $
                kernel = -(1.0-r*0.5/sigma^2)*exp(-r*0.5/sigma^2)/(!pi*sigma^4)
                newimage = convol(fimage,kernel,/edge_truncate)
              end
 'Unsharp' :  begin $
                newimage = unsharp_mask(fimage,radius=sigma)
              end

 else : begin $
          print,'filter type ',ftype,' unknown; check spelling?'
          return,-1
        end
endcase


return,newimage
end
