;NAME:
;  k_arr2d
;PURPOSE:
;  Generate an array of wavenumber values (by default, in units of cycles per 
;  sample interval), arranged corresponding to the FFT of an Nx x Ny-element array.
;CALLING SEQUENCE:
;  k = k_arr2d(Nx, Ny [, dx=dx] [, dy=dy] [, /radians] [, kx=kx] [, ky=ky])
;OUTPUT:
;  k = an Nx x Ny array of scalar wavenumbers. These are arranged just as they would
;     be for a 2D FFT.
;INPUTS:
;  Nx, Ny = dimensions of a hypothetical image for which to calculate the wavenumber.
;OPTIONAL KEYWORD INPUTS:
;  dx = sample interval. Default=1. If included, then k will be converted into units
;     of cycles (or radians, if that keyword is set) per unit distance (or time),
;     corresponding to the units of dx.
;  dy = sample interval along the vertical direction (i.e., along the second
;     array index). Default: dy=dx.
;  radians = if set, then let k be in radians (rather than cycles).
;OPTIONAL KEYWORD OUTPUTS:
;  kx = Nx-dimensional array of horizontal wavenumbers
;  ky = Ny-dimensional array of vertical wavenumbers
;  big_kx, big_ky = Nx x Ny arrays of horizontal and vertical wavenumbers, respectively.
;     This is handy for constructing digital filters as a function of the vector
;     wavenumber, to be applied to the FFT of an image.
;HISTORY:
;  2013-Jun-26  C. Kankelborg
;
function k_arr2d, Nx, Ny, dx=dx, dy=dy, radians=radians, $
   kx=kx, ky=ky, big_kx=big_kx, big_ky=big_ky

;plate scale
if not keyword_set(dx) then dx=1.0
if not keyword_set(dy) then dy=dx

kx = k_arr(Nx, dx=dx, radians=radians) 
ky = k_arr(Ny, dx=dy, radians=radians)

big_kx = kx # replicate(1.0, Ny)
big_ky = replicate(1.0, Nx) # ky

k = sqrt(big_kx^2 + big_ky^2)

return, k
end