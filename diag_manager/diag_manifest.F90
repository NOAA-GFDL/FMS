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

  IMPLICIT NONE

  INTERFACE ASSIGNMENT(=)
     MODULE PROCEDURE manifest_field_type_assign
  END INTERFACE ASSIGNMENT(=)
  
  ! Some type to hold data for manifest
  TYPE manifest_field_type
     CHARACTER(len=128) :: output_name !< output field name in diagnostic file (from diag_table)
     CHARACTER(len=128) :: module_name !< model module that has this field
     CHARACTER(len=128) :: input_name !< field name in model land
     CHARACTER(len=50) :: time_method !< string to hold the time redux method.  If static, the .false.
     INTEGER :: packing !< packing value
     INTEGER :: nDim !< number of dimensions
  END TYPE manifest_field_type

  TYPE manifest_fields_type
     INTEGER :: num_1d = 0
     INTEGER :: num_2d = 0
     INTEGER :: num_3d = 0
     INTEGER :: num_4d = 0
     TYPE(manifest_field_type), DIMENSION(:), ALLOCATABLE :: fields_1d
     TYPE(manifest_field_type), DIMENSION(:), ALLOCATABLE :: fields_2d
     TYPE(manifest_field_type), DIMENSION(:), ALLOCATABLE :: fields_3d
     TYPE(manifest_field_type), DIMENSION(:), ALLOCATABLE :: fields_4d
  END TYPE manifest_fields_type
  
  PRIVATE
  PUBLIC :: write_diag_manifest
  
CONTAINS

  ! PUBLIC routines
  SUBROUTINE write_diag_manifest(file)
    INTEGER, INTENT(in) :: file

    INTEGER :: file_unit, ios
    INTEGER :: num_static, num_temporal
    TYPE(manifest_fields_type) :: static_fields
    TYPE(manifest_fields_type) :: temporal_fields

    CHARACTER(len=128) :: maniFileName

    ! This entire routine should only be called by the rootPE
    IF ( mpp_pe() .EQ. mpp_root_pe() ) THEN
       ! Get the file name.
       !
       ! May need to worry about tile count files(file)%tile_count may have that
       ! information Also need to verify ensembles.  May be better to not use
       ! tile/ensemble as all should have the same data.
       maniFileName = TRIM(files(file)%name)//".mfst"

       static_fields = get_diagnostic_fields(file, static=.TRUE.)
       temporal_fields = get_diagnostic_fields(file, static=.FALSE.)

       ! Get the number of fields to write to manifest file
       num_static = static_fields%num_1d + static_fields%num_2d + static_fields%num_3d + static_fields%num_4d
       num_temporal = temporal_fields%num_1d + temporal_fields%num_2d + temporal_fields%num_3d + temporal_fields%num_4d
       
       ! Open the file for writing, but only if we have something to write
       IF ( num_static + num_temporal .GT. 0 ) THEN
          ! Not using mpp_open, as this routine forces to only write from the root
          ! PE, and each root PE should have its own set of files to write.
          OPEN(UNIT=file_unit, FILE=TRIM(maniFileName), ACCESS='SEQUENTIAL', FORM='FORMATTED',&
               & ACTION='WRITE', POSITION='REWIND', IOSTAT=ios)
          IF ( ios .NE. 0 ) THEN
             CALL error_mesg('diag_manifest_mod::write_diag_manifest',&
                  & 'Unable to open file "'//TRIM(maniFileName)//'".  No manifest file will be created.',&
                  & WARNING)
          ELSE
             ! Open JSON
             write(file_unit,'(A1)') '{'
             ! Fill in other data
             CALL write_manifest(file_unit, static_fields, static=.TRUE.)
             CALL write_manifest(file_unit, temporal_fields, static=.FALSE.)
             ! Close JSON
             write(file_unit,'(A1)') '}'
             !!WRITE(file_unit,'(A128,",",A128,",",A128,",",A50,",",i2,",",i2)') maniField%output_name, manifield%module_name,&
             !!        & maniField%input_name, maniField%time_method, maniField%packing, maniField%nDim
             ! Close the file
             CLOSE(file_unit)
          END IF
       END IF
    END IF
  END SUBROUTINE write_diag_manifest

  ! PRIVATE routines
  SUBROUTINE manifest_field_type_assign(lhs,rhs)
    TYPE(manifest_field_type), INTENT(out) :: lhs
    TYPE(manifest_field_type), INTENT(in) :: rhs

    lhs%output_name = rhs%output_name
    lhs%module_name = rhs%module_name
    lhs%input_name = rhs%input_name
    lhs%time_method = rhs%time_method
    lhs%packing = rhs%packing
    lhs%nDim = rhs%nDim
  END SUBROUTINE manifest_field_type_assign

  SUBROUTINE write_fields(unit, fields)
    INTEGER, INTENT(in) :: unit
    TYPE(manifest_field_type), DIMENSION(:), INTENT(in) :: fields

    INTEGER :: i
    CHARACTER(LEN=*), PARAMETER :: FMT_FLD = "(12X,'""',A,'""',': {')"
    CHARACTER(LEN=*), PARAMETER :: FMT_MOF = "(16X,'""model_field"":','""',A,'"",')"
    CHARACTER(LEN=*), PARAMETER :: FMT_MOD = "(16X,'""module"":','""',A,'"",')"
    CHARACTER(LEN=*), PARAMETER :: FMT_PAK = "(16X,'""packing"":',I1,',')"
    CHARACTER(LEN=*), PARAMETER :: FMT_TAV = "(16X,'""time_averaging"":','""',A,'""')"
    
    DO i=1, SIZE(fields)
       WRITE (unit,FMT_FLD) TRIM(fields(i)%output_name)
       WRITE (unit,FMT_MOF) TRIM(fields(i)%input_name)
       WRITE (unit,FMT_MOD) TRIM(fields(i)%module_name)
       WRITE (unit,FMT_PAK) fields(i)%packing
       WRITE (unit,FMT_TAV) TRIM(fields(i)%time_method)
       IF ( i.EQ.SIZE(fields) ) THEN
          WRITE (unit,'(12X,A1)') '}'
       ELSE
          WRITE (unit,'(12X,A2)') '},'
       END IF
    END DO
  END SUBROUTINE write_fields

  SUBROUTINE write_manifest(unit, fields, static)
    INTEGER, INTENT(in) :: unit
    TYPE(manifest_fields_type), INTENT(in) :: fields
    LOGICAL, INTENT(in) :: static

    CHARACTER(len=*), PARAMETER :: FMT_DIM = "(8X,'""',A2,'""',': {')"
    CHARACTER(len=*), PARAMETER :: FMT_STA = "(4X,'""',A6,'""',': {')"
    CHARACTER(len=*), PARAMETER :: FMT_TEM = "(4X,'""',A8,'""',': {')"
    
    ! Static / Temporal
    IF ( static ) THEN
       WRITE (unit,FMT_STA) 'Static'
    ELSE
       WRITE (unit,FMT_TEM) 'Temporal'
    END IF

    ! 1D fields
    WRITE (unit,FMT_DIM) '1D'
    CALL write_fields(unit, fields%fields_1d(1:fields%num_1d))
    WRITE (unit,'(8X,A2)') '},'

    ! 2D fields
    WRITE (unit,FMT_DIM) '2D'
    CALL write_fields(unit, fields%fields_2d(1:fields%num_2d))
    WRITE (unit,'(8X,A2)') '},'

    ! 3D fields
    WRITE (unit,FMT_DIM) '3D'
    CALL write_fields(unit, fields%fields_3d(1:fields%num_3d))
    WRITE (unit,'(8X,A2)') '},'

    ! 4D fields
    WRITE (unit,FMT_DIM) '4D'
    CALL write_fields(unit, fields%fields_4d(1:fields%num_4d))
    WRITE (unit,'(8X,A1)') '}'

    ! Static / Temporal
    IF ( static ) THEN
       WRITE (unit,'(4X,A2)') '},'
    ELSE
       WRITE (unit,'(4X,A1)') '}'
    END IF
  END SUBROUTINE write_manifest
  
  TYPE(manifest_fields_type) FUNCTION get_diagnostic_fields(file, static)
    INTEGER, INTENT(in) :: file !< diagnostic file, as defined by diag_manager_mod
    LOGICAL, INTENT(in) :: static !< Indicates if looking for static or non-static
                                  ! fields.  .TRUE. indicates looking only for
                                  ! static files.  .FALSE. indicates looking only
                                  ! for non-static fields.
    
    INTEGER :: i, j, o
    INTEGER :: istat
    TYPE(manifest_field_type) :: maniField
    CHARACTER(len=128) :: maniFileName

    maniFileName = TRIM(files(file)%name)//".mfst"
    
    DO j=1, files(file)%num_fields
       o = files(file)%fields(j) ! Position of this field in output_fields array
       IF ( output_fields(o)%written_once .AND. (static.EQV.output_fields(o)%static) ) THEN
          ! output field was written to file, and is static/non-static, whichever was requested
          ! Gather the information to record it.
          i = output_fields(o)%input_field ! Position of the input fields associated with this output_field
             
          ! this is information I currently know we want to save, and where it is:
          maniField%output_name = output_fields(o)%output_name
          maniField%module_name = input_fields(i)%module_name
          maniField%input_name = input_fields(i)%field_name
          IF ( output_fields(o)%static ) THEN
             ! Static fields MUST have a time_method of .false.
             maniField%time_method = ".false."
          ELSE
             maniField%time_method = output_fields(o)%time_method
          END IF
          maniField%packing = output_fields(o)%pack
          maniField%nDim = output_fields(o)%num_axes

          ! Now that we have the information about the field, add to type based on dimensions of field
          SELECT CASE (maniField%nDim)
          CASE (1)
             get_diagnostic_fields%num_1d = get_diagnostic_fields%num_1d + 1
             IF ( .NOT.ALLOCATED(get_diagnostic_fields%fields_1d) ) THEN
                ! Allocate to the max number of fields
                ALLOCATE(get_diagnostic_fields%fields_1d(files(file)%num_fields), STAT=istat)
                IF ( istat.NE.0 ) THEN
                   CALL error_mesg('diag_manifest_mod::get_diagnostic_fields',&
                        & 'Unable to allocate 1d array for manifest file "'//TRIM(maniFileName)//'".  Manifest incomplete.',&
                        & WARNING)
                   ! Resetting count to 0 to keep from writing out
                   get_diagnostic_fields%num_1d = 0
                   CYCLE
                END IF
             END IF
             IF ( ALLOCATED(get_diagnostic_fields%fields_1d) ) THEN
                get_diagnostic_fields%fields_1d(get_diagnostic_fields%num_1d) = maniField
             END IF
          CASE (2)
             get_diagnostic_fields%num_2d = get_diagnostic_fields%num_2d + 1
             IF ( .NOT.ALLOCATED(get_diagnostic_fields%fields_2d) ) THEN
                ! Allocate to the max number of fields
                ALLOCATE(get_diagnostic_fields%fields_2d(files(file)%num_fields), STAT=istat)
                IF ( istat.NE.0 ) THEN
                   CALL error_mesg('diag_manifest_mod::get_diagnostic_fields',&
                        & 'Unable to allocate 2d array for manifest file "'//TRIM(maniFileName)//'".  Manifest incomplete.',&
                        & WARNING)
                   ! Resetting count to 0 to keep from writing out
                   get_diagnostic_fields%num_2d = 0
                   CYCLE
                END IF
             END IF
             IF ( ALLOCATED(get_diagnostic_fields%fields_2d) ) THEN
                get_diagnostic_fields%fields_2d(get_diagnostic_fields%num_2d) = maniField
             END IF
          CASE (3)
             get_diagnostic_fields%num_3d = get_diagnostic_fields%num_3d + 1
             IF ( .NOT.ALLOCATED(get_diagnostic_fields%fields_3d) ) THEN
                ! Allocate to the max number of fields
                ALLOCATE(get_diagnostic_fields%fields_3d(files(file)%num_fields), STAT=istat)
                IF ( istat.NE.0 ) THEN
                   CALL error_mesg('diag_manifest_mod::get_diagnostic_fields',&
                        & 'Unable to allocate 3d array for manifest file "'//TRIM(maniFileName)//'".  Manifest incomplete.',&
                        & WARNING)
                   ! Resetting count to 0 to keep from writing out
                   get_diagnostic_fields%num_3d = 0
                   CYCLE
                END IF
             END IF
             IF ( ALLOCATED(get_diagnostic_fields%fields_3d) ) THEN
                get_diagnostic_fields%fields_3d(get_diagnostic_fields%num_3d) = maniField
             END IF
          CASE (4)
             get_diagnostic_fields%num_4d = get_diagnostic_fields%num_4d + 1
             IF ( .NOT.ALLOCATED(get_diagnostic_fields%fields_4d) ) THEN
                ! Allocate to the max number of fields
                ALLOCATE(get_diagnostic_fields%fields_4d(files(file)%num_fields), STAT=istat)
                IF ( istat.NE.0 ) THEN
                   CALL error_mesg('diag_manifest_mod::get_diagnostic_fields',&
                        & 'Unable to allocate 4d array for manifest file "'//TRIM(maniFileName)//'".  Manifest incomplete.',&
                        & WARNING)
                   ! Resetting count to 0 to keep from writing out
                   get_diagnostic_fields%num_4d = 0
                   CYCLE
                END IF
             END IF
             IF ( ALLOCATED(get_diagnostic_fields%fields_4d) ) THEN
                get_diagnostic_fields%fields_4d(get_diagnostic_fields%num_4d) = maniField
             END IF
          END SELECT
       END IF
    END DO
  END FUNCTION get_diagnostic_fields
END MODULE diag_manifest_mod
