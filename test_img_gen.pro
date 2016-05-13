function test_img_gen, ax, ay, intensity
	
	Nx = 2048
	Ny = 1024
	
	;ax = -8
	;ay = 3
	;intensity = 0.5
	base = 1024

	
	x_low = fix(Nx / 3)
	x_up = fix(2 * x_low)
	
	y_low = fix(Ny / 3)
	y_up = fix(2 * y_low)
	
	image = intarr(Nx, Ny)
	
	for i = x_low, x_up do begin
		for j = y_low, y_up do begin
		
			image[i,j] = base
		
		endfor
	endfor
	
	for i = x_low+ax, x_up+ax do begin
		for j = y_low+ay, y_up+ay do begin
		
			image[i,j] += base * intensity
		
		endfor
	endfor
	
return, image	
end
