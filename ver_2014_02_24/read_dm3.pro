;
;
;modified to work with 3D DM3 stacks
;modified to read the dimensional info and scaling factor

PRO read_dm3_datatype,field_code,bytelength,dm3_data_name,$
                      idl_data_name,idl_variable
  if (field_code eq 2) then begin
     bytelength = 2
     dm3_data_name = 'short'
     idl_data_name = 'int'
     idl_variable = fix(0)
  endif else if (field_code eq 3) then begin
     bytelength = 4
     dm3_data_name = 'long'
     idl_data_name = 'long'
     idl_variable = long(0)
  endif else if (field_code eq 4) then begin
     bytelength = 2
     dm3_data_name = 'ushort'
     idl_data_name = 'uint'
     idl_variable = uint(0)
  endif else if (field_code eq 5) then begin
     bytelength = 4
     dm3_data_name = 'ulong'
     idl_data_name = 'ulong'
     idl_variable = ulong(0)
  endif else if (field_code eq 6) then begin
     bytelength = 4
     dm3_data_name = 'float'
     idl_data_name = 'float'
     idl_variable = float(0)
  endif else if (field_code eq 7) then begin
     bytelength = 8
     dm3_data_name = 'double'
     idl_data_name = 'double'
     idl_variable = double(0)
  endif else if (field_code eq 8) then begin
     bytelength = 1
     dm3_data_name = 'boolean'
     idl_data_name = 'byte'
     idl_variable = byte(0)
  endif else if (field_code eq 9) then begin
     bytelength = 1
     dm3_data_name = 'char'
     idl_data_name = 'byte'
     idl_variable = byte(0)
  endif else if (field_code eq 10) then begin
     bytelength = 1
     dm3_data_name = 'octet'
     idl_data_name = 'byte'
     idl_variable = byte(0)
  endif else begin
     bytelength = 0
     dm3_data_name = 'unknown'
     idl_data_name = 'unknown'
     idl_variable = fix(0)
  endelse

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
PRO read_dm3, filename, imginfo, verbose=verbose


if (n_params() lt 2) then begin
   print,'Usage: read_dm3,filename,image,[/verbose or verbose=2]'
   return
endif

svec = size(filename)
if ((svec[0] ne 0) or (svec[1] ne 7)) then begin
   print,'Usage: read_dm3,filename,image with "filename" an ascii string'
   return
endif

svec = size(verbose)
if (svec[svec[0]+1] eq 0) then verbose = 0

if (verbose ne 0) then print,'Modified version of read_dm3 for 2D/3D images.'

  ;; Establish our own endianness.  Network order is big endian.
endian_test = fix(1)
byteorder,endian_test,/ntohs
if (endian_test eq 1) then begin
   i_am_little_endian = 0
endif else begin
   i_am_little_endian = 1
endelse

;; This is used as a marker in the file
percent_long = long(byte('%%%%'),0)

get_lun,lun
on_ioerror,read_dm3_bailout
openr,lun,filename

dm3_version = long(0)
n_bytes = long(0)
tag_is_little_endian = long(0)

readu,lun,dm3_version,n_bytes,tag_is_little_endian
dm3_version = swap_endian(dm3_version,/swap_if_little_endian)
n_bytes = swap_endian(n_bytes,/swap_if_little_endian)
tag_is_little_endian = swap_endian(tag_is_little_endian,$
                                   /swap_if_little_endian)
if (tag_is_little_endian eq i_am_little_endian) then begin
   swap_tags = 0
endif else begin
   swap_tags = 1
endelse

if (verbose gt 1) then begin
   if (i_am_little_endian eq 0) then begin
      endian_string = 'This computer is big endian '
      if (tag_is_little_endian eq 0) then begin
         endian_string = endian_string+$
                         'and tags are big endian (swap='+$
                         strtrim(string(swap_tags),2)+')'
      endif else begin
         endian_string = endian_string+$
                         'and tags are little endian (swap='+$
                         strtrim(string(swap_tags),2)+')'
      endelse
   endif else begin
      endian_string = 'This computer is little endian '
      if (tag_is_little_endian eq 0) then begin
         endian_string = endian_string+$
                         'and tags are big endian (swap='+$
                         strtrim(string(swap_tags),2)+')'
      endif else begin
         endian_string = endian_string+$
                         'and tags are little endian (swap='+$
                         strtrim(string(swap_tags),2)+')'
      endelse
   endelse
   print,endian_string
endif

tag_groups_are_sorted = byte(0)
tag_groups_are_open = byte(0)
n_tags = long(0)
readu,lun,tag_groups_are_sorted,tag_groups_are_open,n_tags
n_tags = swap_endian(n_tags,/swap_if_little_endian)

if (verbose gt 1) then begin
   print,'Gatan Digital Micrograph file version='+$
         strtrim(string(dm3_version),2)+', size='+$
         strtrim(string(n_bytes),2)+' bytes, with '+$
         strtrim(string(n_tags),2)+' tags.'
endif


;; "Dimensions" holds dimensions of the data

;; Loop through all TagEntries
previous_tag_label = ''
got_nx = 0
got_ny = 0
nx = 0
ny = 0
nz = 0
sc_x = 0.0
sc_y = 0.0
sc_z = 0
got_scx = 0
got_scy = 0
if_first=0
keep_going = 1
i_tag = 0
while (keep_going eq 1) do begin
   tag_group_or_data = byte(0)
   tag_label_length = fix(0)
   
   if (verbose gt 1) then begin
      point_lun,-lun,at_tag
      print,'  Next tag data starts at file offset '+$
            strtrim(string(at_tag),2)
   endif
 ; stop

   readu,lun,tag_group_or_data,tag_label_length
   tag_label_length = swap_endian(tag_label_length,/swap_if_little_endian)
   tag_label = ''
   if (tag_label_length ne 0) then begin
      tag_label = bytarr(tag_label_length)
      readu,lun,tag_label
      tag_label = string(tag_label)
   endif
   if (tag_group_or_data eq 0) then begin
      if (verbose ne 0) then begin
         print,'End of file'
      endif
      goto,read_dm3_close_file
   endif else if (tag_group_or_data eq 20) then begin
      tag_group_is_sorted = byte(0)
      tag_group_is_open = byte(0)
      tag_group_n_tags = long(0)
      readu,lun,tag_group_is_sorted,tag_group_is_open,tag_group_n_tags
      if (verbose ne 0) then begin
         print,'Group tag '+strtrim(string(i_tag+1),2)+' "'+tag_label+$
               '" has '+strtrim(string(tag_group_n_tags),2)+' tags'
      endif
   endif else if (tag_group_or_data eq 21) then begin
      if (verbose ne 0) then begin
         print,'Data tag '+strtrim(string(i_tag+1),2)+' "'+tag_label+'":'
      endif
      percent_marker = bytarr(4)
      readu,lun,percent_marker
      if (long(percent_marker,0) ne percent_long) then begin
         print,'  Instead of "%%%%", got "'+string(percent_marker)+'"'
         goto,read_dm3_bailout
      endif
      tag_type_info_length = long(0)
      readu,lun,tag_type_info_length
      tag_type_info_length = swap_endian(tag_type_info_length,$
                                         /swap_if_little_endian)
      if (tag_type_info_length eq 0) then begin
         ;; Do nothing
      endif else if (tag_type_info_length eq 1) then begin
         ;; We have a single variable
         tag_type = long(0)
         readu,lun,tag_type
         tag_type = swap_endian(tag_type,/swap_if_little_endian)
         read_dm3_datatype,tag_type,bytelength,dm3_data_name,$
                           idl_data_name,idl_variable
         if (bytelength eq 0) then begin
            print,'  Tag type '+strtrim(string(tag_type),2)+' is unknown'
            goto,read_dm3_bailout
         endif else begin
            readu,lun,idl_variable
            if (swap_tags eq 1) then begin
               idl_variable = swap_endian(idl_variable)
            endif
            if (verbose ne 0) then begin
               print,'  Single '+dm3_data_name+' tag value=',idl_variable
            endif
            if (strcmp(tag_label,'Scale') eq 1) then begin
                ;This is either scx,scy,or scz
                if_first+=1
                ;print,if_first
                if (got_scx eq 0 and got_scy eq 0 and if_first gt 1) then begin
                scx = idl_variable
                got_scx=1
                ;print,'Got Dimx = ',scx
                endif else if (got_scy eq 0 and if_first gt 1) then begin
                scy = idl_variable
                got_scy=1
                ;print,'Got Dimy =',scy
                endif else if (if_first gt 1) then begin
                scz=idl_variable
                got_scx=0
                got_scy=0
                ;print,'Got Dimz =',scz
                endif
            endif
            if ((strcmp(previous_tag_label,'Dimensions') eq 1) and $
                (strlen(tag_label) eq 0) and $
                (previous_tag_value eq 2)) then begin
               ;; This is either nx or ny or nz
               if (got_nx eq 0 and got_ny eq 0) then begin
                  nx = idl_variable
                  got_nx = 1
                  if (verbose ne 0) then begin
                     point_lun,-lun,this_file_position
                     print,'    At file offset '+$
                           strtrim(string(this_file_position),2)+$
                           ', got nx='+strtrim(string(nx),2)
                  endif
               endif else if (got_ny eq 0) then begin
                  ny = idl_variable
                  got_ny = 1
                  if (verbose ne 0) then begin
                     point_lun,-lun,this_file_position
                     print,'    At file offset '+$
                           strtrim(string(this_file_position),2)+$
                           ', got ny='+strtrim(string(ny),2)
                  endif
               endif else begin
                 nz = idl_variable
                 got_nx = 0
                 got_ny = 0
                  if (verbose ne 0) then begin
                     point_lun,-lun,this_file_position
                     print,'    At file offset '+$
                           strtrim(string(this_file_position),2)+$
                           ', got nz='+strtrim(string(nz),2)
                  endif
                endelse
            endif
         endelse
      endif else begin
         ;; tag_type_info_length greater than 1
         tag_type = long(0)
         readu,lun,tag_type
         tag_type = swap_endian(tag_type,/swap_if_little_endian)
         
         if (tag_type eq 15) then begin
            ;; Struct
            struct_group_length = long(0)
            struct_n_entries = long(0)
            readu,lun,struct_group_length,struct_n_entries
            struct_group_length = swap_endian(struct_group_length,$
                                              /swap_if_little_endian)
            struct_n_entries = swap_endian(struct_n_entries,$
                                           /swap_if_little_endian)
            if (verbose ne 0) then begin
               print,'  Tag type is struct, with '+$
                     strtrim(string(struct_n_entries),2)+' entries'
            endif
            point_lun,-lun,this_entry_address
            ;; This is where we will have to jump to get the data values
            this_entry_address = this_entry_address+8*struct_n_entries
            for struct_entry=0,(struct_n_entries-1) do begin
               this_field_name_length = long(0)
               this_field_type = long(0)
               readu,lun,this_field_name_length,this_field_type
               this_field_name_length = swap_endian(this_field_name_length,$
                                                    /swap_if_little_endian)
               this_field_type = swap_endian(this_field_type,$
                                             /swap_if_little_endian)
               read_dm3_datatype,this_field_type,bytelength,dm3_data_name,$
                                 idl_data_name,idl_variable               
               ;; Record the present file position
               point_lun,-lun,present_position
               ;; Jump to get the data value
               point_lun,lun,this_entry_address
               readu,lun,idl_variable
               if (swap_tags eq 1) then begin
                  idl_variable = swap_endian(idl_variable)
               endif
               ;; This is where the next data value will be
               this_entry_address = this_entry_address+bytelength
               ;; Go back to where we were in reading entry info
               point_lun,lun,present_position
               if (verbose gt 1) then begin
                  print,'    Entry '+strtrim(string(struct_entry+1),2)+$
                        ': '+dm3_data_name+'='+strtrim(string(idl_variable),2)
               endif
            endfor
            ;; At the end, we need to go beyond the last entry value to
            ;; where the next whole new entry starts
            point_lun,lun,this_entry_address
         endif else if (tag_type eq 18) then begin
            ;; string
            print,'  Tag type is string.  Not yet written'
            goto,read_dm3_bailout
         endif else if (tag_type eq 20) then begin
            ;; Array
            array_type = long(0)
            readu,lun,array_type
            array_type = swap_endian(array_type,$
                                     /swap_if_little_endian)
            if (array_type eq 15) then begin
               ;; Array of groups
               array_group_length = long(0)
               n_dims = long(0)
               readu,lun,array_group_length,n_dims
               array_group_length = swap_endian(array_group_length,$
                                                /swap_if_little_endian)
               n_dims = swap_endian(n_dims,/swap_if_little_endian)
               if (verbose gt 1) then begin
                  print,'  Tag type is group array, with '+$
                        strtrim(string(n_dims),2)+' dimensions'
               endif
               dimensions_agree = 1
               output_string = 'Dimension types: '
               for i_dim=0,(n_dims-1) do begin
                  this_dim_name_length = long(0)
                  this_dim_type = long(0)
                  readu,lun,this_dim_name_length,this_dim_type
                  this_field_name_length = swap_endian(this_field_name_length,$
                                                       /swap_if_little_endian)
                  this_dim_type = swap_endian(this_dim_type,$
                                              /swap_if_little_endian)
                  read_dm3_datatype,this_dim_type,bytelength,dm3_data_name,$
                                    idl_data_name,idl_variable
                  if (i_dim eq 0) then begin
                     dimension_types = [this_dim_type]
                     output_string = output_string+dm3_data_name+', '
                     first_dim_type = this_dim_type
                  endif else begin
                     dimension_types = [dimension_types,this_dim_type]
                     if (i_dim eq (n_dims-1)) then begin
                        output_string = output_string+dm3_data_name
                     endif else begin
                        output_string = output_string+dm3_data_name+', '
                     endelse
                     if (this_dim_type ne first_dim_type) then begin
                        dimensions_agree = 0
                     endif
                  endelse
               endfor
               if (dimensions_agree eq 0) then begin
                  output_string = output_string+' do NOT agree'
                  print,'Data tag '+strtrim(string(i_tag+1),2)+' "'+$
                        tag_label+'":'
                  print,output_string
               endif else begin
                  output_string = output_string+' agree'
                  if (verbose gt 1) then begin
                     print,output_string
                  endif
               endelse
               array_size = long(0)
               readu,lun,array_size
               array_size = swap_endian(array_size,/swap_if_little_endian)
               if (verbose ne 0) then begin
                  print,'    Array dimensions: ['+strtrim(string(n_dims),2)+$
                        ','+strtrim(string(array_size),2)+']'
               endif
               this_array = replicate(idl_variable,n_dims,array_size)
               readu,lun,this_array
               if (swap_tags eq 1) then begin
                  this_array = swap_endian(this_array)
               endif
               if (strcmp(tag_label,'Data') eq 1) then begin
                  print,'2D data ['+strtrim(string(n_dims),2)+','+$
                        strtrim(string(array_size),2)+']'
                  data_2d = data_array
               endif
            endif else begin
               ;; 1D array
               array_size = long(0)
               readu,lun,array_size
               array_size = swap_endian(array_size,$
                                        /swap_if_little_endian)
               read_dm3_datatype,array_type,bytelength,dm3_data_name,$
                                 idl_data_name,idl_variable
               if (verbose ne 0) then begin
                  print,'    1D array with '+$
                        strtrim(string(array_size),2)+' '+$
                        dm3_data_name+'-type values (IDL '+idl_data_name+')'
               endif
               if (array_size ne 0) then begin
                  if (bytelength eq 0) then begin
                     print,'  Array data type is unknown'
                     goto,read_dm3_bailout
                  endif else begin
                     data_array = replicate(idl_variable,array_size)
                     readu,lun,data_array
                     if (swap_tags eq 1) then begin
                        data_array = swap_endian(data_array)
                     endif
                  endelse
               endif
               if (strcmp(tag_label,'Data') eq 1) then begin
                  if (verbose ne 0) then begin
                     print,'1D data ['+strtrim(string(array_size),2)+']'
                  endif
                  data_1d = data_array
               endif
               if (verbose gt 1) then begin
                  point_lun,-lun,this_file_position
                  print,'    Now at byte position '+$
                        strtrim(string(this_file_position),2)+' in the file'
               endif
            endelse
         endif
      endelse
   endif 
   if (strlen(tag_label) ne 0) then begin
      previous_tag_label = tag_label
      previous_tag_value = idl_variable
   endif
   i_tag = i_tag+1
endwhile


read_dm3_close_file:
if (nz eq 0) then begin
image=reform(data_1d,nx,ny)
  dims=[nx,ny]
  scale=[scx,scy]
  imginfo={img2d,dims:dims,scale:scale,image:image}
endif else begin 
image = reform(data_1d,nx,ny,nz)
dims=[nx,ny,nz]
scale=[scx,scy,scz]
imginfo={img3d,dims:dims,scale:scale,image:image}
endelse
read_dm3_bailout:
close,lun
free_lun,lun
on_ioerror,null

end
