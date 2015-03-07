;This is designed to work at 3.2 sec using VCF reprojected from 250m MODIS and 3.2sec ALOS.


PRO mask_line, hlorey_block, vcf_block, hv_block, work_line
	;use work_line, even though it is part of hlorey_block because it may not be in the center of the block
	;also because it will be modified for output

	;we can obtain win_size from the y-dimension of the blocks
	sizeinfo = size(hlorey_block)
	win_size = sizeinfo[2]
	xdim = sizeinfo[1]

	;set hlorey to water (-1) where it is NaN
	index = where((finite(work_line) eq 0), count)
	if (count gt 0) then work_line[index] = -1

	;set hv = 0 to water
	index = where(hv_block[*,win_size/2] eq 0, count)
	if (count gt 0) then work_line[index] = -1

	;set hv < 0.005 to 0
	index = where(hv_block[*,win_size/2] lt 0.005, count)
	if (count gt 0) then work_line[index] = 0

	;work on first win_size/2 columns
	index = where((vcf_block[win_size,*] eq 200), count)
	if (count gt 0) then begin
		for i=0ULL, win_size/2 do begin
			if (hv_block[i,win_size/2] lt 0.01) then work_line[i] = -1
		endfor
	endif
	
	;work on middle portion
	for i=ulong(win_size)/2+1, xdim-win_size/2-1 do begin
		index = where((vcf_block[i-(win_size/2):i+(win_size/2),*] eq 200), count)
		if (count gt 0) then begin
			if (hv_block[i,win_size/2] lt 0.01) then work_line[i] = -1
		endif
	endfor

	;work on remaining columns
	index = where((vcf_block[xdim-win_size:xdim-1,*] eq 200), count)
	if (count gt 0) then begin	
		for i=xdim-win_size/2, xdim-1 do begin
			if (hv_block[i,win_size/2] lt 0.01) then work_line[i] = -1
		endfor
	endif	

End


PRO multi_mask, hlorey_file, vcf_file, hv_file, xdim, ydim, out_file


	;Settings ----------
	win_size = 5ULL    ;For 5x5 window   win_size should be an odd number greater than or equal to 3

	hlorey_block = fltarr(xdim, win_size)
	vcf_block = bytarr(xdim, win_size)
	hv_block = fltarr(xdim, win_size)

	hlorey_line = fltarr(xdim)
	vcf_line = bytarr(xdim)
	hv_line = fltarr(xdim)
	work_line = fltarr(xdim)
	;-------------------

	;print configurations....
	print, 'hlorey_file:', hlorey_file
	print, 'vcf_file:', vcf_file
	print, 'hv_file:', hv_file
	print, 'xdim, ydim:', xdim, ydim
	print, 'out_file:', out_file
	print, 'masking file...', systime()

	;Open files for processing
	openr, hlorey_lun, hlorey_file
	openr, vcf_lun, vcf_file
	openr, hv_lun, hv_file

	openw, out_lun, out_file	


	;Read and work on first win_size lines
	readu, hlorey_lun, hlorey_block
	readu, vcf_lun, vcf_block
	readu, hv_lun, hv_block
	
	for j=0ULL, 	win_size/2 do begin
		work_line = hlorey_block[*,j]
		mask_line, hlorey_block, vcf_block, hv_block, work_line
		writeu, out_lun, work_line
	endfor

	;middle part

	for j=ulong(win_size)/2+1, ydim-(win_size/2)-1 do begin

		if (j mod 1000 eq 0) then print, j
		readu, hlorey_lun, hlorey_line
		readu, vcf_lun, vcf_line
		readu, hv_lun, hv_line

		for i=0, win_size-2 do begin
			hlorey_block[*,i] = hlorey_block[*,i+1]
			vcf_block[*,i] = vcf_block[*,i+1]
			hv_block[*,i] = hv_block[*,i+1]
		endfor
		hlorey_block[*,win_size-1] = hlorey_line
		vcf_block[*,win_size-1] = vcf_line
		hv_block[*,win_size-1] = hv_line

		work_line[*] = hlorey_block[*,win_size/2]
		mask_line, hlorey_block, vcf_block, hv_block, work_line
		writeu, out_lun, work_line
	endfor



	;Read and work on last win_size lines
	;we have already read to the end of file, so we just need to keep the blocks the same, and work through remaining lines
	
	for j=0ULL, win_size/2-1  do begin
		work_line[*] = hlorey_block[*,win_size/2+j+1]
		mask_line, hlorey_block, vcf_block, hv_block, work_line
		writeu, out_lun, work_line
	endfor


	free_lun, hlorey_lun, vcf_lun, hv_lun, out_lun

	print, 'Done!', systime()

End
