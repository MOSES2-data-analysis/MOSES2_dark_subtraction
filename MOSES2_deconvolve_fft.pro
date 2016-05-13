; Deconvolve using Fourier Translation theorem
; If the measured data is represented as
; f(t) = g(t) + b*g(t + a)
; where g(t) is the single image (the target), 
; f(t) is the measured image
; b is a scalar representing the relative intensity of the peaks
; a is a vector representing the spatial distance between the peaks

function moses2_deconvolve_fft, img, ax, ay, intensity

	sz = size(img)	
	
	Nx = sz[1]
	Ny=sz[2]

	;ax = -6
	;ay = -3
	;intensity = 0.8

	; Load MOSES2 level one dataset into memory
	;restore, '/nfs/home/roysmart/mosesLevelOne.sav'
	
	K = k_arr2d(Nx, Ny, /radians, big_kx=kx, big_ky=ky)
	
	
	; temp
	;i = 16
	J = complex(0,1) ; imaginary unit
	
	
	
	; Loop through data cube to deconvolve each exposure
	;for i = 0, Ndata-1 do begin
	
		;xtv, cube_zero[*,*,i], screenwidth=1800, screenheight=900
	
		zero_ft = fft(img)
		
		
			
			zero_ft /= (1 + intensity * exp(-J*(kx*ax+ky*ay)))
			

		
		img = real_part(fft(zero_ft, /INVERSE))
		
	 	xtv, img, screenwidth=1800, screenheight=900
	 
	;endfor

return, img

end
