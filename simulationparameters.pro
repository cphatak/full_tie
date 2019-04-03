;+
;function that returns all the parameters required in the
;plate tfs,phase,tie,recons routines
;@author charudatta phatak,cmu,09/28/06
;modified: 12/21/06: all types of configs now possible. No. of sides, type of
;magnetization i.e, closure or uniform
;
; renamed to simulationparameters.pro by MDG 5/8/07
;modified by CD, ANL. This version has the correct prefactors for phase reconstruction.
;-

; this will require a bit more work to make it more generally useful.  
; the basic procedure should be such that this routine provides default
; values for all entries;  default values can be changed when passed as a named
; argument....

function simulationparameters,help=help,dim=dim,delta=delta,EE=EE,cs=cs,thetac=thetac,defstep=defstep,spread=spread

;--------------------------
; image-related parameters
image = { dim:512, $		; image dimension (dimxdim)
	  delta:1.D0} 		; pixel size [nm]
if keyword_set(dim) then begin
	image.dim = dim
	print,'image dimension redefined to ',image.dim
endif
if keyword_set(delta) then begin
	image.delta= delta
	print,'image pixel size redefined to ',image.delta
endif
if keyword_set(help) then begin
	print,'The following parameters are the default settings for the simulationparameters routine;'
	print,' Quantities preceded by a star (*) can be changed as optional arguments'
	print,' '
	print,'*image dimension                 : ',image.dim
	print,'*pixel size                [nm]  : ',image.delta
end 
pimage=ptr_new(image)

; use this value to define the tie structure below
dim = long((*pimage).dim)

;--------------------------
;physical constants
physcon={e:1.602177D-19, $	; electron charge [C]
	 c:2.99792458D+8, $	; velocity of light [m/s]
	 h:6.626075D-34, $	; Planck's constant [Js]
	 k:1.380658D-23, $	; Boltzmann's constant [J/K]
	 mu0:1.2566371D-06, $	; vacuum permeability [-]
	 mue:9.284770D-24, $	; electron magnetic moment [J/T]
	 eps0:8.854188D-12, $	; vacuum permittivity [F/m]
	 phi0:2.067834D-15, $	; flux quantum [Wb]
	 m:9.109389D-31}  	; electron rest mass [kg]
if keyword_set(help) then begin
	print,'electron charge            [C]   : ',physcon.e
	print,'velocity of light          [m/s] : ',physcon.c
	print,'Planck''s constant          [Js]  : ',physcon.h
	print,'Boltzmann''s constant       [J/K] : ',physcon.k
	print,'vacuum permeability        [-]   : ',physcon.mu0
	print,'electron magnetic moment   [J/T] : ',physcon.mue
	print,'vacuum permittivity        [F/m] : ',physcon.eps0
	print,'flux quantum               [Wb]  : ',physcon.phi0
	print,'electron rest mass         [kg]  : ',physcon.m
end 
pphyscon=ptr_new(physcon)

;--------------------------
;microscope parameters
scope={EE:200000.0D, $ 	; accelerating voltage in Volt
       cs:1.0D+6, $	; spherical aberration [nm]
       thetac:6.0D-4, $	; beam divergence [radian]
       defstep:200.0D,$	; defocus step size [nm]
       spread:8.0D, $   ; defocus spread [nm]
       lambda:0.0D,$	; electron wavelength [nm, to be computed]
       omega:0.0D, $	; relativistic factor [to be computed]
       gamma:0.0D, $	; relativistic factor [to be computed]
       sigma:0.0D}	; interaction constant [to be computed]

if keyword_set(EE) then scope.EE= EE
if keyword_set(cs) then scope.cs= cs
if keyword_set(thetac) then scope.thetac= thetac
if keyword_set(defstep) then scope.defstep= defstep
if keyword_set(spread) then scope.spread= spread
epsilon=0.5D0*physcon.e/physcon.m/physcon.c^2
scope.lambda=physcon.h*1.0D+9/sqrt(2.0D*physcon.m*physcon.e)/sqrt(scope.EE+epsilon*scope.EE^2)
scope.omega=physcon.e*scope.EE/physcon.m/physcon.c^2
gama=1.D0+scope.omega
scope.gamma = gama
scope.sigma=2.0*!pi*physcon.m*scope.gamma*physcon.e*scope.lambda*1.0D-18/physcon.h^2
if keyword_set(help) then begin
	print,'*accelerating voltage      [V]   : ',scope.EE
	print,'*spherical aberratin       [nm]  : ',scope.cs
	print,'*beam divergence           [rad] : ',scope.thetac
	print,'*defocus step size         [nm]  : ',scope.defstep
	print,'*defocus spread            [nm]  : ',scope.spread
	print,'electron wavelength        [nm]  : ',scope.lambda
	print,'relativistic factor gamma  [-]   : ',scope.gamma
	print,'interaction constant    [1/V/nm] : ',scope.sigma
end 
pscope=ptr_new(scope)

if keyword_set(help) then goto,skiptie
;--------------------------
; define a new structure to hold all the variables relevant to the TIE reconstruction
;tie={   ima_dIdz:dblarr(dim,dim),  $    ; longitudinal intensity derivative (overfocus - underfocus)
;	ima_if:dblarr(dim,dim),$        ; in-focus image
;	phase:dblarr(dim,dim),$         ; reconstructed phase
;	qi:dblarr(dim,dim),$            ; inverse Laplacian array (precomputed)
;	gamx:dblarr(dim,dim),$          ; inverse gradient x-component (precomputed)
;	gamy:dblarr(dim,dim),$          ; inverse gradient y-component (precomputed)
;	qis:dblarr(2*dim,2*dim),$       ; symmetrized inverse Laplacian array (precomputed)
;	gamxs:dblarr(2*dim,2*dim),$     ; symmetrized inverse gradient x-component (precomputed)
;	gamys:dblarr(2*dim,2*dim),$     ; symmetrized inverse gradient y-component (precomputed)
;	pre:0.0D0,$                     ; pre-factor for Laplacian reconstruction
;	preb:0.0D0,$                    ; pre-factor for Bx, By reconstruction
;	preg:0.0D0,$                    ; pre-factor for gradient reconstruction
;	bxt:dblarr(dim,dim),$           ; reconstructed bx integrated induction (intermediate result)
;	byt:dblarr(dim,dim)}            ; reconstructed by integrated induction (intermediate result)
tie={   ima_dIdz:fltarr(dim,dim),  $    ; longitudinal intensity derivative (overfocus - underfocus)
	ima_if:fltarr(dim,dim),$        ; in-focus image
	phase:fltarr(dim,dim),$         ; reconstructed phase
	qi:fltarr(dim,dim),$            ; inverse Laplacian array (precomputed)
	gamx:fltarr(dim,dim),$          ; inverse gradient x-component (precomputed)
	gamy:fltarr(dim,dim),$          ; inverse gradient y-component (precomputed)
	qis:fltarr(2*dim,2*dim),$       ; symmetrized inverse Laplacian array (precomputed)
	gamxs:fltarr(2*dim,2*dim),$     ; symmetrized inverse gradient x-component (precomputed)
	gamys:fltarr(2*dim,2*dim),$     ; symmetrized inverse gradient y-component (precomputed)
	pre:0.0D0,$                     ; pre-factor for Laplacian reconstruction
	preb:0.0D0,$                    ; pre-factor for Bx, By reconstruction
	preg:0.0D0,$                    ; pre-factor for gradient reconstruction
	prebg:0.0D0,$			; pre-factor for Bx,By in gradient method reconstruction
	bxt:fltarr(dim,dim),$           ; reconstructed bx integrated induction (intermediate result)
	byt:fltarr(dim,dim)}            ; reconstructed by integrated induction (intermediate result)

skiptie:
if (not keyword_set(help)) then begin
	ptie=ptr_new(tie)
; define the scaling factors for the TIE reconstructions
	;(*ptie).pre=-1.0/(*pimage).delta^2/(64.D0*!dpi^3*(*pscope).lambda*(*pscope).defstep)
	;(*ptie).preg=-1.0/(*pimage).delta^2/(4.D0*!dpi*(*pscope).lambda*(*pscope).defstep) /4.D0
	;(*ptie).preb=(*ptie).preg * (*pphyscon).h/(2.D0*!dpi*(*pphyscon).e)*1.0D18 *(*pimage).delta   ; to get T nm for integrated induction
	;(*ptie).preb=-1.0/((*pimage).delta*8.D0*!dpi*(*pscope).lambda*(*pscope).defstep)
	siconstant = (*pphyscon).h/(2.D0*!dpi*(*pphyscon).e)*1.0D18 ; to get induction in T nm.
	(*ptie).pre = -1.0*(*pimage).delta^2/(64.D0*!dpi^3*(*pscope).lambda*(*pscope).defstep)
	(*ptie).preb = -1.0*(*pimage).delta/(8.D0*!dpi*(*pscope).lambda*(*pscope).defstep) * siconstant
	(*ptie).preg = -1.0*(*pimage).delta^2/(4.D0*!dpi*(*pscope).lambda*(*pscope).defstep)
	(*ptie).prebg = 1.0/(2.D0*(*pimage).delta)*siconstant
	return, [pimage,pscope,pphyscon,ptie]
end else begin
	; all pointers should be released, since we do not need them in the calling program
	ptr_free, ptr_valid()
	return,1
endelse
end
