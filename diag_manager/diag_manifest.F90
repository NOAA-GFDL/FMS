MODULE diag_manifest_mod

  USE diag_data_mod, ONLY: files,&  ! TYPE(file_type) --- diagnostic files
       & output_fields,& ! TYPE(output_field_type) --- field  in diagnostic file
       & input_fields ! TYPE(input_field_type) --- field from diag_table
  USE mpp_io_mod, ONLY: mpp_open,&
       & MPP_OVERWR,&
       & MPP_ASCII,&
       & MPP_SEQUENTIAL,&
       & MPP_SINGLE, &
       & mpp_close
  USE mpp_mod, ONLY: mpp_pe,&
       & mpp_root_pe
  USE fms_mod, ONLY: error_mesg,&
       & WARNING
  
  ! Some type to hold data for manifest
  TYPE manifest_field_type
     CHARACTER(len=128) :: output_name !< output field name in diagnostic file (from diag_table)
     CHARACTER(len=128) :: module_name !< model module that has this field
     CHARACTER(len=128) :: input_name !< field name in model land
     CHARACTER(len=50) :: time_method !< string to hold the time redux method.  If static, the .false.
     INTEGER :: packing !< packing value
     INTEGER :: nDim !< number of dimensions
  END TYPE manifest_field_type

  PRIVATE
  PUBLIC :: write_diag_manifest
  
CONTAINS

  SUBROUTINE write_diag_manifest(file)
    INTEGER, INTENT(in) :: file

    INTEGER :: i, j, o
    INTEGER :: file_unit, ios
    TYPE(manifest_field_type) :: maniField
    CHARACTER(len=128) :: maniFileName

    ! This entire routine should only be called by the rootPE
    IF ( mpp_pe() .EQ. mpp_root_pe() ) THEN
       ! Get the file name.
       !
       ! May need to worry about tile count files(file)%tile_count may have that
       ! information Also need to verify ensembles.  May be better to not use
       ! tile/ensemble as all should have the same data.
       maniFileName = TRIM(files(file)%name)//".mfst"

       ! Open the file for writing
       !
       ! Not using mpp_open, as this routine forces to only write from the root
       ! PE, and each root PE should have its own set of files to write.
       OPEN(UNIT=file_unit, FILE=TRIM(maniFileName), ACCESS='SEQUENTIAL', FORM='FORMATTED',&
            & ACTION='WRITE', POSITION='REWIND', IOSTAT=ios)
       IF ( ios .NE. 0 ) THEN
          CALL error_mesg('diag_manifest_mod::write_diag_manifest',&
               & 'Unable to open file "'//TRIM(maniFileName)//'".  No manifest file will be created.',&
               & WARNING)
       ELSE
          ! Loop over all fields in file
          DO j = 1, files(file)%num_fields
             o = files(file)%fields(j) ! Position of this field in output_fields array
             i = output_fields(o)%input_field ! Position of the input fields associated with this output_field
             
             ! this is information I currently know we want to save, and where it is:
             maniField%output_name = output_fields(o)%output_name
             maniField%module_name = input_fields(i)%module_name
             maniField%input_name = input_fields(i)%field_name
             IF ( output_fields(o)%static ) THEN
                maniField%time_method = ".false."
             ELSE
                maniField%time_method = output_fields(o)%time_method
             END IF
             maniField%packing = output_fields(o)%pack
             maniField%nDim = output_fields(o)%num_axes
             ! Write the data to the manifest file.  Cannot use mpp_write as it expects
             ! time dependent fields
             IF ( output_fields(o)%written_once ) THEN
                WRITE(file_unit,'(A128,",",A128,",",A128,",",A50,",",i2,",",i2)') maniField%output_name, manifield%module_name,&
                     & maniField%input_name, maniField%time_method, maniField%packing, maniField%nDim
             END IF
          END DO
          
          ! Close the file
          CLOSE(file_unit)
       END IF
    END IF
  END SUBROUTINE write_diag_manifest
END MODULE diag_manifest_mod
