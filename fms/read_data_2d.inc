
   if (do_init) call fms_init ( )

   if (.not.associated(Domain)) call error_mesg &
        ('read_data in fms_mod', 'set_domain not called', FATAL)

       if (present(end)) end = .false. 

       if (do_read()) read (unit,end=10) gdata
       if (.not.read_all_pe) then
         ! did not read on all PEs need to broadcast data
           len = size(gdata,1)*size(gdata,2)
           call mpp_transmit ( gdata(isg,jsg), len, ALL_PES, &
                               gdata(isg,jsg), len, 0        )       
       endif   
     ! return data for compute domain
       data(is:ie,js:je) = gdata(is:ie,js:je)
       return  

   10  call read_eof (end)

