
   if (.not.module_is_initialized) call fms_io_init ( )

   if (.not.associated(Current_domain)) call mpp_error &
        (FATAL,'set_domain not called')

       if (present(end)) end = .false. 

       if (do_read()) read (unit,end=10) gdata
       if (.not.read_all_pe) then
         ! did not read on all PEs need to broadcast data
           len = size(gdata,1)*size(gdata,2)
!           call mpp_transmit ( gdata, len, ALL_PES, &
!                               gdata, len, 0        ) 
           call mpp_broadcast ( gdata, len, mpp_root_pe() )      
       endif   
     ! return data for compute domain
       data(is:ie,js:je) = gdata(is:ie,js:je)
       return  

   10  call read_eof (end)

