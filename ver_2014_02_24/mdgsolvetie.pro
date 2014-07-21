;
; mdgsolvetie.pro, a general routine for solving the Transport-of-Intensity Equation
;
;
; routines are based on several versions written since 1999
;
; this version, 05/08/2007, written at Argonne National Laboratory
;
;
; improvements to be made:
;
;	- currently, the derivatives are three-point derivatives, computed
;	  via convolutions with an appropriate kernel.  Perhaps we could 
;	  replace those by actual derivatives using Fourier transforms...
;	  In that case, we must be careful with the multiplicative prefactors !!!
;
;--------------------------------


;--------------------------------
function symmetrize,image	; performs the TIE image symmetrization
;--------------------------------

sz = size(image,/dimensions)
dim = 2*sz[0]

;imi = dblarr(dim,dim)
;imi[0,dim/2] = image  
;imi[dim/2,dim/2] = -reverse(image,1)
;imi[0,0] = -reverse(image,2)  
;imi[dim/2,0] = reverse(reverse(image,1),2)
;imi = shift(imi,0,dim/2)

imi=dblarr(dim,dim)
imi[0:dim/2-1,0:dim/2-1] = image
imi[0:dim/2-1,dim/2:dim-1] = -reverse(image,2)
imi[dim/2:dim-1,0:dim-1] = reverse(imi[0:dim/2-1,0:dim-1],1)


return,imi
end
;--------------------------------


;--------------------------------
pro cdreconstruct,ptie,method,symmetrize=symmetrize
compile_opt idl2
;--------------------------------

; get the array size for proper scaling of ffts
sz = size((*ptie).ima_dIdz,/dimensions)
dim = sz[0]

if keyword_set(symmetrize) then dim = 2*dim

; kernels for gradient 
kx=double([[0.0,0.0,0.0],[-1.0,0.0,1.0],[0.0,0.0,0.0]])
ky=double([[0.0,1.0,0.0],[0.0,0.0,0.0],[0.0,-1.0,0.0]])

if keyword_set(symmetrize) then begin 
; symmetrize both the in-focus image and the intensity derivative
  dIdz= symmetrize((*ptie).ima_dIdz)
  imh = 1.D0/symmetrize((*ptie).ima_if)
  if (method eq 'laplacian') then begin
	  qi = (*ptie).qis
  end else begin
	  gamx = (*ptie).gamxs
	  gamy = (*ptie).gamys
  endelse
end else begin
  dIdz= (*ptie).ima_dIdz
  imh = 1.D0/(*ptie).ima_if
  if (method eq 'laplacian') then begin
	  qi = (*ptie).qi
  end else begin
	  gamx = (*ptie).gamx
	  gamy = (*ptie).gamy
  endelse
endelse

; compute Fourier transform of longitudinal intensity derivative
fd=fft(dIdz,-1,/double)

if (method eq 'laplacian') then begin
; apply first inverse Laplacian operator
  tmp=-float(fft(fd*qi,1,/double))

; apply gradient operator and divide by in-focus image
  dtmpdx=convol(tmp,kx,/center,/edge_wrap) * imh
  dtmpdy=convol(tmp,ky,/center,/edge_wrap) * imh

; at this point we have the integrated induction as an intermediate result
  if keyword_set(symmetrize) then begin
  	(*ptie).byt=(*ptie).preb * dtmpdx[0:dim/2-1,0:dim/2-1] * float(dim)^2
  	(*ptie).bxt=(*ptie).preb * dtmpdy[0:dim/2-1,0:dim/2-1] * float(dim)^2
  end else begin
  	(*ptie).byt=(*ptie).preb * dtmpdx * float(dim)^2
  	(*ptie).bxt=(*ptie).preb * dtmpdy * float(dim)^2
  endelse

; apply second gradient operator
  dtmpdx=convol(dtmpdx,kx,/center,/edge_wrap)
  dtmpdy=convol(dtmpdy,ky,/center,/edge_wrap)

; combine gradient components and apply second inverse Laplacian
  fd=fft(dtmpdx+dtmpdy,-1,/double)
  tmp=-fft(fd*qi,1,/double)

; scale and put minimum to zero
  if keyword_set(symmetrize) then begin
  	(*ptie).phase = float(tmp[0:dim/2-1,0:dim/2-1]) * (*ptie).pre * float(dim)^4 
  end else begin
  	(*ptie).phase = float(tmp) * (*ptie).pre * float(dim)^4 
  endelse
  ;(*ptie).phase -= min((*ptie).phase)

end else begin ; reconstruct the phase using the inverse gradient method

;apply the inverse gradient, divide by in-focus image, and inverse gradient again
  dtmpdx = fft(gamx*fd,1,/double)*imh
  dtmpdy = fft(gamy*fd,1,/double)*imh
  phix = double(fft(fft(dtmpdx,-1,/double)*gamx,1,/double))
  phiy = double(fft(fft(dtmpdy,-1,/double)*gamy,1,/double))

; scale and put min(phase) at zero, since origin of phase is not defined
  dtmpdx = phix+phiy
  if keyword_set(symmetrize) then begin
  	(*ptie).phase = -(*ptie).preg * dtmpdx[0:dim/2-1,0:dim/2-1] * float(dim)^2
  end else begin
  	(*ptie).phase = -(*ptie).preg * dtmpdx * float(dim)^2
  endelse
  ;(*ptie).phase -= min((*ptie).phase)

; compute integrated induction from phase and use the proper scaling factors to make sure the units are T nm 
  (*ptie).byt=convol((*ptie).phase,kx,/center,/edge_wrap)*(*ptie).prebg
  (*ptie).bxt=convol((*ptie).phase,ky,/center,/edge_wrap)*(*ptie).prebg
endelse
end
;--------------------------------


;+
;program to recover the phase from an aligned through focus series
;based on the TIE formalism.
;@author charudatta phatak, cmu, 09/28/06
; modified by MDG 5/7/2007
;-

;--------------------------------
pro mdgsolvetie,method,ptie,skip=skip,symmetrize=symmetrize,tikhonov=tikhonov
compile_opt idl2
;--------------------------------

; the first time this routine is called, all auxiliary arrays will be computed;
; all subsequent calls should use the keyword /skip to skip this portion so that 
; the computation will be a little faster.

if keyword_set(skip) then goto,skipall

;get the image dimensions
sz = size((*ptie).phase,/dimensions)
dim = sz[0]
if keyword_set(symmetrize) then dim = 2*dim 

; define arrays for inverse gradient reconstruction 
line = (dindgen(dim)-double(dim/2))
gamy = rebin(reform(line,1,dim),dim,dim)
gamx = rotate(gamy,3)
g = gamx^2+gamy^2
g[dim/2,dim/2] = 1.D0
g = 1.D0/g
g[dim/2,dim/2] = 0.D0
if keyword_set(symmetrize) then begin
	(*ptie).gamxs = shift(gamx*g,dim/2,dim/2)
	(*ptie).gamys = shift(gamy*g,dim/2,dim/2)
end else begin 
	(*ptie).gamx = shift(gamx*g,dim/2,dim/2)
	(*ptie).gamy = shift(gamy*g,dim/2,dim/2)
endelse

; define array for inverse laplacian reconstruction 
q=double(dist(dim))
q[0,0]=1.D0
if keyword_set(tikhonov) then begin	
	gc = double(tikhonov);1.5D0
	print,gc
	qi = q^2/(q^2+gc^2)^2
end else begin
	qi=1.D0/q^2
endelse

qi[0,0]=0.D0
if keyword_set(symmetrize) then begin
	(*ptie).qis=qi
end else begin 
	(*ptie).qi=qi
endelse

; remove variables that are no longer needed
  gamx = 0 & gamy = 0 & q =0 & q = 0 & line = 0 & qi = 0

skipall:

; make sure the in-focus image has no zeroes in it
s=where((*ptie).ima_if lt 0.001)
if (s[0] ne -1) then (*ptie).ima_if[s]=0.001

; reconstruct the phase using inverse gradient or inverse laplacian methods, symmetrized or not
if keyword_set(symmetrize) then cdreconstruct,ptie,method,/symmetrize else cdreconstruct,ptie,method
end
;--------------------------------


