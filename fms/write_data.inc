
   if (do_init) call fms_init ( )

   if (.not.associated(Domain)) call error_mesg &
        ('write_data in fms_mod', 'set_domain not called', FATAL)

 ! get the global data and write only on root pe
   call mpp_global_field ( Domain, data, gdata ) 
   if ( mpp_pe() == mpp_root_pe() ) write (unit) gdata

