;+
;NAME:
;  MOSES_PREP
;PURPOSE:
;  Routine MOSES data reduction. Extract and dark subtract
;  all the data, and apply flat fields. Save as 'mosesLevelOne.sav'.
;	It is assumed that imageindex.xml is in the current working
;	directory.
;		Saturated pixels are flagged as missing (missing_data=NaN).
;	Pixels that are negative after background subtraction are mapped
;	into the domain [0, 0.01] using the POSITIVITY function from the
;	IDL Astronomy User's Library.
;		There is a script embedded at bottom that lets MOSES_PREP run
;	as a batch job under ssw_batch.
;CALLING SEQUENCE:
;	From IDL,
;  	IDL> moses_prep [, /byteorder]
;	or from the shell,
;		% unlimit
;		% ssw_batch moses_prep moses_prep.log /date
;OPTIONAL INPUT KEYWORDS:
;	byteorder --- passed straight through to moses_read. This should be
;		set for big-endian machines such as PPC, IRIX and Sparc. It should
;		not be set for x86 or Dec.
;HISTORY:
;  2006-May-20 C. Kankelborg  based on movie0.pro.
;	2006-May-22 C. Kankelborg  fixed byteorder problem.
;	2006-May-23 C. Kankelborg  flag saturated pixels as missing data;
;		map negative pixels into [0, 0.01].
;	2006-May-25 C. Kankelborg  fixed a bug that swapped the plus
;		and zero orders.
;	2006-Jun-17 C. Kankelborg  Fixed sun orientation (flip E-W).
;		Added (best available, old) flat fields for minus and zero orders.
;		Switched to NaN for missing data (formerly used 0.0).
;	2006-Jun-19 C. Kankelborg  Added filtering (IMSCRUB) after application
;		of flat fields to detect and correct spikes, dust, etc.
;	2006-Jun-21 C. Kankelborg  Removed the /all from the save statement.
;		It creates incompatibilities with other IDL versions! Repaired
;		existing save file by loading and resaving.
;  2012-Feb-02 C. Kankelborg  Added detector 3 (m=+1) flat field, finally!
;  2012-Feb-13 C. Kankelborg  Commented out /byteorder flag from the
;     ssw_batch script at the bottom.
;-
pro moses_prep, byteorder=byteorder

  print, systime()+' MOSES_PREP starting.'
  
  
  dt = 250000.0 ;exposure time resolution (microsec)
  forward_function hist_equal_long
  
  ;Set up dataspace
  Nx = 2048
  Ny = 1024
  
  index = mxml2('mdata/imageindex.xml',curdir())
  
  N_images = n_elements(index.filename)
  
  for i=0,N_images-1 do begin
  
    ;TALLY DATA IMAGES
    if (strpos(index.seqname[i],'data') ne -1) then begin
      print,'Found DATA image ',index.seqname[i]
      add_element, data_list, i
    endif
    
    ;TALLY DARK IMAGES
    if (strpos(index.seqname[i],'dark') ne -1) then begin
      print,'Found DARK image ',index.seqname[i]
      add_element, dark_list, i
    endif
    
    ;CREATE/UPDATE EXPOSURE LIST, exp_list
    if n_elements(exp_list) eq 0 then exp_list = [index.exptime[i]] else begin
      ss = where(abs(exp_list - index.exptime[i]) lt index.exptime[i]/10.0)
      if ss eq -1 then add_element, exp_list, $
        dt*round(index.exptime[i]/dt)
    endelse
    
  endfor
  
  ;PROCESS DARK IMAGES
  Nexptimes = n_elements(exp_list)
  darks_minus = fltarr(Nx, Ny, Nexptimes) ;array of darks (DN)
  darks_zero  = fltarr(Nx, Ny, Nexptimes) ;array of darks (DN)
  darks_plus  = fltarr(Nx, Ny, Nexptimes) ;array of darks (DN)
  
  darkRepeats = intarr(Nexptimes) ;initialized to zero
  print,'Processing dark images...'
  Ndark = n_elements(dark_list)
  for j = 0, Ndark-1 do begin
    i = dark_list[j]
    error=0L
    moses2_read, index.filename[i],minus,zero,plus,noise,          $
      directory=curdir(), error=error, $
      byteorder=byteorder
    if error eq 0L then begin
      print,"Got ",index.filename[i]
      ;locate exposure time in table
      s = where(abs(exp_list - index.exptime[i]) lt index.exptime[i]/10.0)
      if s[0] eq -1 or n_elements(s) ne 1 then message,'Exposure list error.'
      s = s[0] ;make it into a scalar rather than a 1-element array.
      ;incorporate zero order image into average dark (DN).
      darks_minus[*,*,s] = ( darkRepeats[s]*darks_minus[*,*,s] + $
        minus ) / (1 + darkRepeats[s])
      darks_zero[*,*,s] = ( darkRepeats[s]*darks_zero[*,*,s] + $
        zero ) / (1 + darkRepeats[s])
      darks_plus[*,*,s] = ( darkRepeats[s]*darks_plus[*,*,s] + $
        plus ) / (1 + darkRepeats[s])
      darkRepeats[s] = darkRepeats[s] + 1
    endif
  endfor
  
  
  ;PROCESS DATA IMAGES
  Ndata = n_elements(data_list)
  cube_minus = fltarr(Nx,Ny,Ndata)
  cube_zero  = fltarr(Nx,Ny,Ndata)
  cube_plus  = fltarr(Nx,Ny,Ndata)
  for j = 0, Ndata-1 do begin
    i = data_list[j]
    error=0L
    moses2_read, index.filename[i],minus,zero,plus,noise,          $
      directory=curdir(), error=error, $
      byteorder=byteorder
    if error eq 0L then begin
      print,"Got ",index.filename[i]
      print,"Data ranges for minus, zero, plus:"
      pmm, minus
      pmm, zero
      pmm, plus
      
      ;find saturated pixels
      saturated = 2^14 - 1  ;corresponds to saturation of 14-bit ADC
      sat_minus = where(minus eq saturated)
      sat_zero  = where(zero  eq saturated)
      sat_plus  = where(plus  eq saturated)
      
      ;locate exposure time in table
      s = where(abs(exp_list - index.exptime[i]) lt index.exptime[i]/10.0)
      if s[0] eq -1 or n_elements(s) ne 1 then message,'Exposure list error.'
      s = s[0] ;make it into a scalar rather than a 1-element array.
      
      ;Dark subtract and map to all positive values
      zero_value = 0.1 ;the value zero will map to
      minus = positivity(minus - darks_minus[*,*,s], epsilon=4*zero_value^2)
      zero  = positivity(zero  - darks_zero[*,*,s],  epsilon=4*zero_value^2)
      plus  = positivity(plus  - darks_plus[*,*,s],  epsilon=4*zero_value^2)
      
      ;Mark saturated pixels as missing data
      if (sat_minus[0] ne -1) then minus[sat_minus] = missing_data
      if (sat_zero[0]  ne -1) then zero[sat_zero]   = missing_data
      if (sat_plus[0]  ne -1) then plus[sat_plus]   = missing_data
      
      ;Store in cube, flipped E-W to view sun correctly.
      cube_minus[*,*,j] = rotate(minus, 5)
      cube_zero[*,*,j]  = rotate(zero, 5)
      cube_plus[*,*,j]  = rotate(plus, 5)
    endif
  endfor
  
  
  ;Flat fields
  print, systime()+' Applying flat fields to all spectral orders....'
  restore, "flats/det1ff_042605-2nd.sav"
  flat_minus  = rotate(q/max(q), 7)
  ;normalize to unit maximum and orient to the Sun.
  
  restore, "flats/det2ff_042705-2nd.sav"
  flat_zero = rotate(q/max(q), 7)
  ;normalize to unit maximum and orient to the Sun.
  
  restore, "flats/det3ff_111606-2nd.sav
  flat_plus = rotate(q/max(q), 7)
  ;normalize to unit maximum and orient to the Sun.
  
  for i=0, Ndata-1 do begin
    cube_minus[*,*,i] /= flat_minus
    cube_zero[*,*,i]  /= flat_zero
    cube_plus[*,*,i]  /= flat_plus
  endfor
  
  
  ;Post-filtering to remove dust, particle hits, etc.
  print, systime()+' Post-filtering with IMSCRUB....'
  
  ;Set up some storage to track bad pixel maps for all 3 data cubes:
  cube_minus_badpix = bytarr(Nx, Ny, Ndata)
  cube_zero_badpix  = bytarr(Nx, Ny, Ndata)
  cube_plus_badpix  = bytarr(Nx, Ny, Ndata)
  
  ;Parameters for IMSCRUB2:
  mindiff    = 1.0
  ndev       = 3.0
  thresh     = 0.25
  medwidth   = 5
  expand_bad = 1
  print,"For the record: mindiff=",mindiff,", ndev=", ndev,", thresh=",thresh
  print,"    medwidth=",medwidth,",  expand_bad=",expand_bad
  
  for i=0, Ndata-1 do begin
    print,'Exposure ',i,' order m = -1'
    cube_minus[*,*,i] = imscrub2(cube_minus[*,*,i], image_marked=image_marked, $
      medwidth = medwidth, expand_bad=expand_bad,                             $
      mindiff=mindiff, ndev=ndev, thresh=thresh)
    cube_minus_badpix[*,*,i] = ~finite(image_marked)
    
    print,'Exposure ',i,' order m =  0'
    cube_zero[*,*,i]  = imscrub2(cube_zero[*,*,i], image_marked=image_marked, $
      medwidth = medwidth, expand_bad=expand_bad,                            $
      mindiff=mindiff, ndev=ndev, thresh=thresh)
    cube_zero_badpix[*,*,i] = ~finite(image_marked)
    
    print,'Exposure ',i,' order m = +1'
    cube_plus[*,*,i]  = imscrub2(cube_plus[*,*,i], image_marked=image_marked, $
      medwidth = medwidth, expand_bad=expand_bad,                            $
      mindiff=mindiff, ndev=ndev, thresh=thresh)
    cube_plus_badpix[*,*,i] = ~finite(image_marked)
  endfor
  
  
  print, systime()+' MOSES_PREP saving to disk.'
  save, filename='mosesLevelOne.sav'
  print, systime()+' MOSES_PREP completed.'
end



;*** Here is the script that allows execution in batch mode: ***
set_plot,'z'

print, !version
moses_prep    ;,/byteorder
end

