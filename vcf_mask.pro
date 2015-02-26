;Apply MODIS vcf mask to input file.
;in_file and mask_file should have the same dimension.  mask is assumed to be byte and infile float


PRO vcf_mask, in_file, mask_file, out_file, threshold

	print, 'Performing vcf mask:', systime()
	print, 'Input file:', in_file
	print, 'Output file:', out_file
	print, 'Mask file:', mask_file
	print, 'Threshold:', threshold

	openr, in_lun, in_file, /get_lun
	openr, mask_lun, mask_file, /get_lun

	openw, out_lun, out_file, /get_lun

	infile_info = file_info(in_file)
	infile_size = infile_info.size

	tot_pix = infile_size/4

	read_pix = tot_pix / 100

	remainder_pix = tot_pix mod 100

	in_line = fltarr(read_pix)
	in_mask = bytarr(read_pix)
	out_line = fltarr(read_pix)

	for i=0, 99 do begin
		print, i
		readu, in_lun, in_line
		readu, mask_lun, in_mask
		

		index = where(in_mask lt threshold, count)
		if (count gt 0) then in_line[index] = 0

		index = where(in_mask eq 200, count)  ;water
		if (count gt 0) then in_line[index] = -1

		writeu, out_lun, in_line
	endfor

	if (remainder_pix gt 0) then begin

		in_line = fltarr(remainder_pix)
		in_mask = bytarr(remainder_pix)
		out_line = fltarr(remainder_pix)

		readu, in_lun, in_line
		readu, mask_lun, in_mask
			

		index = where(in_mask lt threshold, count)
		if (count gt 0) then in_line[index] = 0

		index = where(in_mask eq 200, count)  ;water
		if (count gt 0) then in_line[index] = -1

		writeu, out_lun, in_line
	endif

	free_lun, in_lun
	free_lun, mask_lun
	free_lun, out_lun	

	print, 'Done!', systime()
End
