pro slicemap
	;step 1: find slice pixels
	ind = 1
	seg = 13
	margin = 5.
	smsize = 1.2
	fitsfile = '../../../N/N.fits'
	rmsfile = '../../../N/13CO_rms.fits'
	sklfile = '../../../width/allwidth/skl.fits'
	unit = '10!U22!N cm!U-2!N'
	file_name = file_search("../slices/coordinate/idl/" + strtrim(string(ind), 2) + '/' + $
		strtrim(string(seg), 2) + '/*')
	num = n_elements(file_name)
	picfile1 = '../slicemap/slicemap.eps'
	;picfile2 = 'slicemap_modified.eps'
	readcol, '../contalist/contalist' + strtrim(string(ind), 2) + '.cat', seg_num, category, slice, $
		format = 'A, I, A'
	remain = fltarr(5)
	for i = 0, 4 do remain[i] = strmid(slice[seg], i, 1)
	print, remain
	allslice = [1, 1, 1, 1, 1]
	x = fltarr(num, 21)
	y = fltarr(num, 21)
	for i = 0, num - 1 do begin
		readcol, '../../width/allwidth/slices/coordinate/idl/' + strtrim(string(ind), 2) + '/' + $
			strtrim(string(seg), 2) + '/slice' + strtrim(string(i), 2) + '.cat', $
				slice_x, slice_y, format = 'f, f'
		x[i, *] = slice_x
		y[i, *] = slice_y
	endfor
	;step 2: identify map sizes
	readcol, '../../width/allwidth/segment/' + strtrim(string(ind), 2) + '/segment' + strtrim(string(seg), 2) + '.cat', $
		segment_x, segment_y, format = 'f, f'
	x_min = min(segment_x) - margin
	y_min = min(segment_y) - margin
	x_max = max(segment_x) + margin
	y_max = max(segment_y) + margin
	mid_x = (x_min + x_max)/2.
	mid_y = (y_min + y_max)/2.
	x_size = x_max - x_min
	y_size = y_max - y_min
	l = max(x_size, y_size)/2.
	xmin = round(mid_x - l); - 0.5
	xmax = round(mid_x + l); + 0.5
	ymin = round(mid_y - l); - 0.5
	ymax = round(mid_y + l); + 0.5
	;step 3: draw
	fits_read, fitsfile, data, hdr
	fits_read, rmsfile, rms
	cgloadct, 71, /reverse; 39
	stretch, 0, 255, 0.7
	newdata = data[xmin:xmax, ymin:ymax]
	max = max((newdata[where(finite(newdata))]))
	min = min(newdata[where(finite(newdata))])
	;max = 9e22;max(data[where(data gt 0)])*0.8
	min = min(newdata[where(newdata gt 0)])
	xticks = 3
	s = size(newdata, /dim)
	pic_aspect_ratio = double(s[0])/s[1]
	fits_read, sklfile, skl
	sk = skl[xmin:xmax, ymin:ymax]
	indices = array_indices(skl, where(skl gt 0))
	for i = 0, n_elements(indices[0, *]) - 1 do begin
		for j = 0, n_elements(x[*, 10]) - 1 do begin
			if indices[0, i] eq x[j, 10] and indices[1, i] eq y[j, 10] then begin
				indices[0, i] = -1
				indices[1, i] = -1
				break
			endif
		endfor
	endfor
	sk_x = indices[0, *]
	sk_y = indices[1, *]
	cgps_open, picfile1, xsize = 9, ysize = 9, /encapsulated, /portrait  
	!p.font = -1
	!p.thick = 2
	!p.charthick = 3
	!p.CHARSIZE = 2.
	;print, img
	pos = [0.2,0.2,0.9,0.7]
	cgimage, newdata, position=pos, /save, /KEEP_ASPECT_RATIO
	cgplot, [0], [0], xrange = [xmin-0.5, xmax+0.5], yrange = [ymin-0.5, ymax+0.5], /nodata, psym = 4, position = pos, /noerase, $
		xtitle = '!17 x (pixel)', ytitle = '!17 y (pixel)', xtickinterval = 10,  ytickinterval = 10, xminor = 5, yminor = 5
	cgcolorbar, position = [pos[0], pos[3]+0.01, pos[2], pos[3]+0.03],range = [min/1e22, max/1e22],$
    	title = unit, charsize = !p.CHARSIZE, textthick = !P.CharThick, ncolors = 254, $
    	AnnotateColor='black', /top
    cgplot, sk_x, sk_y, psym = 7, symsize = smsize,thick = 5, /overplot, color = 'navy'
    cgplot, x[*, 10], y[*, 10], psym = 16, symsize = smsize*0.8, /overplot, color = 'navy'
    for i = 0, num - 1 do begin
    	if remain[i] eq 1 then cgplot, x[i, *], y[i, *], psym = 6, thick = 5, symsize = smsize, /overplot, color = 'dark_green'
    	if remain[i] eq 0 then cgplot, x[i, *], y[i, *], psym = 6, thick = 5, symsize = smsize, /overplot, color = 'black'
    endfor
    cgps_close

end
