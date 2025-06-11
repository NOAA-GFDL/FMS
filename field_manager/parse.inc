!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************


  character(*), parameter :: SPACE = ' '
  character(*), parameter :: DELIM = SPACE//','
  integer :: parse
  integer :: is, ie, id, k
  integer :: ts, last, i

  parse = 0; ts = 1; last=len_trim(text)
  do
    i=scan(text(ts:last),'=') ! location of the next equal sign in the input test
    if (i == 0) return
    ! find the last non-space character before the equal sign
    do ie = ts+i-2,ts,-1
       if (scan(text(ie:ie),SPACE)==0) exit
    enddo
    ! find the last delimeter preceding spaces and equal sign
    do is = ie,ts,-1
       if (scan(text(is:is),DELIM)>0) exit
    enddo
    if (trim(label)==text(is+1:ie)) exit ! from outer loop: found the label
    ! for the next iteration of the loop
    ts = ts+i+1 ! shift the beginning of the line
  enddo

  is = ts+i
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
     read ( text(is:ie), *, err=99, end=99 ) values(k)
     parse = parse+1   ! parse is the number of values decoded

     if (ie == last) exit  ! end of text string ... DONE
     is = is+id            ! move starting pointer after ","
     if (is  > last) exit  ! end of text string ... DONE
  enddo
  return

99 continue
   call mpp_error (FATAL,'in parse, error decoding "'//trim(label)//'" in text "'//text//'"')
