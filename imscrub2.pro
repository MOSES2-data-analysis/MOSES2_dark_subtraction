;+
;NAME:
; IMSCRUB2
;PURPOSE:
; Remove spikes, lines, cosmic ray hits, etc. from an otherwise
; smooth image. Bad pixels are identified by contrast between
; the image and a median-smoothed image. Then the bad pixels are
; replaced from a smoothed image created by ck_convol (which handles
; bad pixels more flexibly than IDL's convol).
;CALLING SEQUENCE:
; result = imscrub2(image [, mindiff=mindiff] [, ndev=ndev]  $
;   [, thresh=thresh] [, medwidth=medwidth]                 $
;   [, image_marked=image_marked] )
;INPUT PARAMETERS:
; image --- 2d float or double array. If there are negative pixels,
;   imscrub may not behave well.
;OPTIONAL KEYWORD INPUTS:
; mindiff --- fixed threshold for declaring a bad pixel. The
;   default is 1.0.
; ndev --- sqrt(smoothed signal) scaled threshold for declaring
;   a bad pixel. If the signal is in counts, and the noise is
;   mainly Poisson, and the number of counts is large, then ndev
;   is approximately the number of standard deviations. Default = 3.0.
; thresh --- intensity-normalized threshold for declaring a
;   bad pixel.  Default = 0.5.
; medwidth --- width of the median neighborhood. Default = 5 pixels.
; expand_bad --- if set, expand the width of the bad pixel regions
;   by as many pixels as the keyword is set to. Default = 1.
;OPTIONAL KEYWORD OUTPUTS:
; image_marked --- copy of image with bad pixels set to NaN.
;PROCEDURE:
; First the median-smoothed image, median, is calculated.
; Three significance criteria are used for identifying bad pixels
; from the absolute difference map, difference = abs(image - median).
; All three criteria must evaluate true to deem the image pixel bad:
;   1. difference greater than a fixed threshold, mindiff.
;   2. difference greater than shot noise, ndev * sqrt(median)
;   3. difference greater than scaled threshold, thresh * median.
;MODIFICATION HISTORY:
; 2006-Jun-19  C. Kankelborg  based n imscrub.pro. Added tripartite
;   bad pixel test.
;-
function imscrub2, image, thresh=thresh, medwidth=medwidth, $
    image_marked=image_marked, expand_bad=expand_bad, ndev=ndev, $
    mindiff=mindiff
    
  isize = size(image)
  Nx = isize[1]
  Ny = isize[2]
  NaN = !values.f_nan ;used for marking bad data.
  maxit = 16 ;maximum number of smoothing iterations
  kernel = [[1,2,1],[2,4,2],[1,2,1]]/16.0 ;fairly minimal smoothing kernel
  
  ;Defaults for keywords that control identification of bad pixels:
  if (n_elements(mindiff)    eq 0) then    mindiff = 1.0
  if (n_elements(ndev)       eq 0) then       ndev = 3.0
  if (n_elements(thresh)     eq 0) then     thresh = 0.5
  if (n_elements(medwidth)   eq 0) then   medwidth = 5
  if (n_elements(expand_bad) eq 0) then expand_bad = 1
  
  ;Create median smoothed image:
  med = median(image, medwidth)
  
  ;Find and mark bad pixels.
  difference = abs(image - med)
  badpix = float( ( difference gt mindiff        )  and   $
    ( difference gt ndev*sqrt(med) )  and   $
    ( difference gt thresh*med     )        )
  ;Expand bad pixel map if desired.
  if (expand_bad ne 0) then begin
    for i=1,expand_bad do begin
      badpix = convol(badpix, kernel, /edge_truncate)
    endfor
  endif
  ss = where(badpix)
  ss_good = where(~badpix)
  image_marked = image
  if (ss[0] eq -1) then begin
    message,'No bad pixels found.',/informational
    return, image
  endif else begin
    print,'Found ',n_elements(ss),' bad pixels.'
  endelse
  image_marked[ss] = NaN
  
  ;Replace bad pixels as appropriate
  smoothed = image_marked
  for i = 1L,maxit do begin
    print,'   Smoothing iteration ',i
    smoothed = ck_convol(smoothed, kernel, method='redeem_taint', /edge_truncate)
    smoothed[ss_good] = image[ss_good] ;no need to smooth good pixels!
    if min(finite(smoothed[ss])) then break
  endfor
  result = image_marked
  result[ss] = smoothed[ss]
  
  return,result
end
