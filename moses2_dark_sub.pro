; NAME: MOSES2_dark_subtraction
; PURPOSE: Subtract dark current from the MOSES2 images
; OPTIONAL ARGUMENTS:
;	- index_dir = the location of imageindex.xml in the filesystem
;

pro moses2_dark_sub, index_dir=index_dir

  !EXCEPT=2 ; report line numbers of arithmetic errors

  ; If the optional command line argument is not defined
  ; Set the director to the defaulf flight directory
  if ~keyword_set(index_dir) then begin
    index_dir = 'IDLWorkspace/MOSES2_LevelZero_Dataset'
  endif
  
  print, systime() + ' MOSES2 Dark Subtraction Starting...'
  
  ; Initialize variable used for missing data
  NaN = !values.f_nan
  missing_data = NaN
  
  ; Define Dataspace
  Nx = 2048
  Ny = 1024
  
  ; Open index file to find filenames
  index = mxml2('imageindex.xml',index_dir)
  N_images = n_elements(index.filename)
  
  ; Catalog darks and data images
  for i = 0, N_images - 1 do begin
  
    ; Locate dark images
    if(strpos(index.seqname[i], 'dark') ne -1) then begin
      print, 'Found new DARK image: ', index.seqname[i], ', ', index.filename[i], ', ', index.exptime[i]
      add_element, dark_list, i
    endif
    
    ; Locate data images
    if(strpos(index.seqname[i], 'data') ne -1) then begin
      print, 'Found new DATA image: ', index.seqname[i], ', ', index.filename[i] , ', ', index.exptime[i]
      add_element, data_list, i
    endif
    
  endfor
  
  
  
  ; PROCESS DARK IMAGES
  Ndark = 7 ; Not all of the dark images are good in this dataset
  
  ; Store the average dark image in DN/usec
  ; dark_minus = fltarr(Nx, Ny) eliminated since data is lost here
  dark_zero = intarr(Nx, Ny, Ndark)
  dark_plus = intarr(Nx, Ny, Ndark)
  
  ; Find the average value of the dark images
  
  for j = 0, Ndark - 1 do begin
    i = dark_list[j]	; find index of next dark image
    
    ; Copy next image into memory
    moses2_read, index.filename[i],minus,zero,plus,noise, directory=index_dir, error=error, byteorder=byteorder
    
    dark_zero[*,*,j] = zero
    dark_plus[*,*,j] = plus
    
  endfor
  
  ; Calculate the median of the dark cube
  med_zero = median(dark_zero, DIMENSION=3)
  med_plus = median(dark_plus, DIMENSION=3)
  
  
  
  ; Print statistics of median dark
  print, 'ZERO-ORDER MEDIAN DARK'
  print, 'Max = ',max(med_zero)
  print, 'Mean = ',mean(med_zero)
  print, 'Median = ',median(med_zero)
  print, 'Min = ',min(med_zero)
  print
  print, 'POSITIVE-ORDER MEDIAN DARK'
  print, 'Max = ', max(med_plus)
  print, 'Mean = ', mean(med_plus)
  print, 'Median = ', median(med_plus)
  print, 'Min = ',min(med_plus)
  
  ;  xtv, med_zero, screenwidth = 1800, screenheight = 900
  ;  xtv, med_plus, screenwidth = 1800, screenheight = 900
  
  ; DATA PROCESSING
  ; Cubes for storing dark-subtracted data
  Ndata = n_elements(data_list)
  cube_plus = fltarr(Nx, Ny, Ndata)
  cube_zero = fltarr(Nx, Ny, Ndata)
  
  ; Subtract average value of dark images from data
  for j = 0, Ndata - 1 do begin
    i = data_list[j]	; find index of next data image
    
    ; Copy next image into memory
    moses2_read, index.filename[i],minus,zero,plus,noise, directory=index_dir, error=error, byteorder=byteorder

    ; Write images as tiff for testing
    output_dir = 'IDLWorkspace/MOSES2_LevelOne_Dataset/MOSES2_tiff_images_darksub'
    ;output_dir = 'IDLWorkspace/temp'
    zero_tiff = bytscl(sqrt(rotate(zero,5)))
    plus_tiff = bytscl(sqrt(rotate(plus,5)))
    write_tiff, output_dir+'/zero/'+strmid(index.filename[i],7,12)+string((index.exptime[i]*1e-6), FORMAT='(f20.2)')+'  darksub.tif', zero_tiff
    write_tiff, output_dir+'/plus/'+strmid(index.filename[i],7,12)+string((index.exptime[i]*1e-6), FORMAT='(f20.2)')+'  darksub.tif', plus_tiff

    
    ; Print image statistics
    print,"Got ",index.filename[i]
    print,"Data ranges for minus, zero, plus:"
    pmm, minus
    pmm, zero
    pmm, plus
    
    ;find saturated pixels
    saturated = 2^14 - 1  ;corresponds to saturation of 14-bit ADC
    sat_zero  = where(zero  eq saturated)
    sat_plus  = where(plus  eq saturated)
    
    ;Dark subtract and map to all positive values
    zero_value = 0.1 ;the value zero will map to
    zero  = positivity(zero - med_zero, epsilon = 4*zero_value^2)
    plus  = positivity(plus - med_plus, epsilon = 4*zero_value^2)
    
    ;Mark saturated pixels as missing data
    if (sat_zero[0]  ne -1) then zero[sat_zero]   = missing_data
    if (sat_plus[0]  ne -1) then plus[sat_plus]   = missing_data
    
    ;Store in cube, flipped E-W to view sun correctly.)
    cube_zero[*,*,j]  = rotate(zero, 5)
    cube_plus[*,*,j]  = rotate(plus, 5)
    
    ; Write images as tiff for testing
    output_dir = 'IDLWorkspace/MOSES2_LevelZero_Dataset/MOSES2_tiff_images_raw'
    ;output_dir = 'IDLWorkspace/temp'
    zero_tiff = bytscl(sqrt(cube_zero[*,*,j]))
    plus_tiff = bytscl(sqrt(cube_plus[*,*,j]))
    write_tiff, output_dir+'/zero/'+strmid(index.filename[i],7,12)+string((index.exptime[i]*1e-6), FORMAT='(f20.2)')+'  raw.tif', zero_tiff
    write_tiff, output_dir+'/plus/'+strmid(index.filename[i],7,12)+string((index.exptime[i]*1e-6), FORMAT='(f20.2)')+'  raw.tif', plus_tiff

    
  endfor
  
  ;Flat fields
  print, systime()+' Applying flat fields to all spectral orders....'
  restore, "IDLWorkspace/MOSES2_dark_subtraction/flats/det2ff_042705-2nd.sav"
  flat_zero = rotate(q/max(q), 7)
  ;normalize to unit maximum and orient to the Sun.
  
  restore, "IDLWorkspace/MOSES2_dark_subtraction/flats/det3ff_111606-2nd.sav"
  flat_plus = rotate(q/max(q), 7)
  ;normalize to unit maximum and orient to the Sun.
  
  for i=0, Ndata-1 do begin
    cube_zero[*,*,i]  /= flat_zero
    cube_plus[*,*,i]  /= flat_plus
  endfor
  
  
 
  
  
  ;Post-filtering to remove dust, particle hits, etc.
  print, systime()+' Post-filtering with IMSCRUB....'
  
  ;Set up some storage to track bad pixel maps for all 3 data cubes:
  cube_zero_badpix  = bytarr(Nx, Ny, Ndata)
  cube_plus_badpix  = bytarr(Nx, Ny, Ndata)
  
  ;Parameters for IMSCRUB2:
  mindiff    = 1.0
  ndev       = 5.0
  thresh     = 0.25
  medwidth   = 5
  expand_bad = 1
  print,"For the record: mindiff=",mindiff,", ndev=", ndev,", thresh=",thresh
  print,"    medwidth=",medwidth,",  expand_bad=",expand_bad
  
  for i=0, Ndata-1 do begin
    
    ; imscrub does not work for images with low S/
    j = data_list[i] ; find index of next data imageN
    if (index.exptime[j] gt 2e6) then begin
      print, index.exptime[j]
      
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
    
    endif
    
  endfor
  
  ; Save images in TIFF format for easy viewing
  print, 'Saving images to TIFF format...'
  output_dir = 'IDLWorkspace/MOSES2_LevelOne_Dataset/MOSES2_tiff_images_level1'
  for i = 0, Ndata-1 do begin

      ;if index.filename[j] > 5e6 then begin
      j = data_list[i] ; find index of next data image
    
      ; Scale into 8 bit channel
      zero_tiff = bytscl(sqrt(cube_zero[*,*,i]))
      plus_tiff = bytscl(sqrt(cube_plus[*,*,i]))
    
      write_tiff, output_dir+'/zero/'+strmid(index.filename[j],7,12)+string((index.exptime[j]*1e-6), FORMAT='(f20.2)')+'  level1.tif', zero_tiff
      write_tiff, output_dir+'/plus/'+strmid(index.filename[j],7,12)+string((index.exptime[j]*1e-6), FORMAT='(f20.2)')+'  level1.tif', plus_tiff

    ;endif
  endfor

  print, 'images written successfully'
  

  
;  ; Double image removal
;  for i=0, Ndata-1 do begin
;  
;  	;xtv, sqrt(cube_zero[*,*,i]), screenwidth=1800, screenheight=900
;  	cube_zero[*,*,i]=moses2_deconvolve_fft(cube_zero[*,*,i],7,1,0.4) 
;  	;xtv, sqrt(cube_zero[*,*,i]), screenwidth=1800, screenheight=900
;  	
;  endfor
;  
;  ; Save images in TIFF format for easy viewing
;  print, 'Saving images to TIFF format...'
;  output_dir = 'IDLWorkspace/MOSES2_LevelOne_Dataset/MOSES2_tiff_images_level1'
;  for i = 0, Ndata-1 do begin
;  
;    j = data_list[i] ; find index of next data image
;    
;    ; Scale into 8 bit channel
;    zero_tiff = bytscl(cube_zero[*,*,i])
;    
;    write_tiff, output_dir+'/zero/'+strmid(index.filename[j],6,18)+'.corrected.tif', zero_tiff
;    
;  endfor
  
  print, 'images written successfully'
  
  print, systime()+' MOSES_PREP saving to disk.'
  save, filename='IDLWorkspace/MOSES2_LevelOne_Dataset/mosesLevelOne.sav'
  print, systime()+' MOSES_PREP completed.'
  
 
  
end
