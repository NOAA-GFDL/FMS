
integer :: parse
integer :: is, ie, id, k, last

  parse = 0

  last = len_trim(text)
  is = index( trim(text), trim(label) )
  if (is == 0) return

! position of initial starting pointer
  is = is + len_trim(label)            ! move starting pointer after label
  is = is + scan( text(is:last), '=' ) ! move starting pointer after "="

  do k = 1, size(values(:))

   ! position of ending pointer
     id = scan( text(is:last), ',' )
     if (id == 0) then
         ie = last     ! no trailing comma, use end of text string
     else
         ie = is+id-2
     endif

   ! decode value
   ! print *, 'k,is,ie,id,last=',k,is,ie,id,last
   ! print *, 'DECODE: ', text(is:ie)
     read ( text(is:ie), *, err=99 ) values(k)
     parse = parse+1   ! parse is the number of values decoded

     if (ie == last) exit  ! end of text string ... DONE
     is = is+id            ! move starting pointer after ","
     if (is  > last) exit  ! end of text string ... DONE
  enddo
  return
  
99 continue
   call mpp_error (FATAL,'in parse, error decoding data') 

