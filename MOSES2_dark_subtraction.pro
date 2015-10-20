; NAME: MOSES2_dark_subtraction
; PURPOSE: Subtract dark current from the MOSES2 images
; 

pro MOSES2_dark_subtraction, index_dir

JOURNAL, 'MOSES2_dark_subtraction' + systime()  ; Store output in log file

print, systime() + 'MOSES2 Dark Subtraction Starting...'

; Initialize variable used for missing data
NaN = !values.f_nan

; Define Dataspace
Nx = 2048
Ny = 1024

; Open index file to find filenames
index = mxml2('imageindex.xml',index_dir)
N_images = n_elements(index.filename)

; Catalog darks and data images
for i = 0, N_images - 1 do begin

  ; Darks
  if(strpos(index.seqname[i], 'dark') ne -1) then begin
    print, 'Found new DARK image ', index.seqname[i]
    
    
endfor