!***********************************************************************
!*                   GNU Lesser General Public License
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* FMS is free software: you can redistribute it and/or modify it under
!* the terms of the GNU Lesser General Public License as published by
!* the Free Software Foundation, either version 3 of the License, or (at
!* your option) any later version.
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
!* FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
!* for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with FMS.  If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************

!> \author Seth Underwood <Seth.Underwood@noaa.gov>
!!
!! \brief diag_manifest_mod writes out a manifest file for each diagnostic output
!!        file defined in the diag_table file.
!!
!!  diag_manifest_mod writes a JSON formatted manifest file for each diagnostic
!!  file defined in the diag_table file.  The manifest file contains basic
!!  information about each field.  Although, this manifest file is for use in the
!!  future Chaco release of the FMS Runtime Environment (FRE), others may find the
!!  information in this file useful.
!!
!!  Although some FMS components write diagnostic files separated by tiles
!!  (Cubed-sphere atmosphere), and some models are run with multiple ensembles the
!!  only one manifest file will be written for each.  That is, although an
!!  atmos_cubed_sphere component may write `atmos_month.tile[1-6].nc`, only one
!!  manifest file `atmos_month.mfst` will be written.  This was done as
!!  diag_manager_mod does not allow a tile or ensemble to write out a different
!!  set of diagnostics.  All tiles, and ensemble members read the same diag_table
!!  file.
MODULE diag_manifest_mod

  USE diag_data_mod, ONLY: files,&  ! TYPE(file_type) --- diagnostic files
       & output_fields,& ! TYPE(output_field_type) --- field  in diagnostic file
       & input_fields,& ! TYPE(input_field_type) --- field from diag_table
       & prepend_date,& ! LOGICAL --- indicates if the date should be prepended to files
       & diag_init_time ! TYPE(time_type) -- model time when diag_manager initialized
  USE mpp_mod, ONLY: mpp_pe,&
       & mpp_root_pe,&
       & get_unit,& ! Get a good file unit value
       & mpp_npes,& ! Get number of PEs in pelist
       & mpp_gather
  USE fms_mod, ONLY: error_mesg,&
       & WARNING
  USE fms_io_mod, ONLY: get_filename_appendix
  USE time_manager_mod, ONLY: get_date

  IMPLICIT NONE

  !> \brief Assignment operator for TYPE(manifest_field_type)
  !!
  !! Allow the TYPE(manifest_field_type) to be assigned properly.  In most cases,
  !! this shouldn't be needed, but it is added here just in case some compiler
  !! just doesn't want to do the correct thing.
  INTERFACE ASSIGNMENT(=)
     MODULE PROCEDURE manifest_field_type_assign
  END INTERFACE ASSIGNMENT(=)

  !> \brief A type to hold the data required for the manifest file.
  !!
  !! The data collected in this type is directly from the other types used in
  !! diag_manager, namely: output_fields and input_fields.
  TYPE manifest_field_type
     CHARACTER(len=128) :: output_name !< output field name in diagnostic file (from diag_table)
     CHARACTER(len=128) :: module_name !< model module that has this field
     CHARACTER(len=128) :: input_name !< field name in model land
     CHARACTER(len=50) :: time_method !< string to hold the time redux method.  If static, the .false.
     INTEGER :: packing !< packing value
     INTEGER :: nDim !< number of dimensions
  END TYPE manifest_field_type

  !> \brief A type to hold all the fields by dimension size
  !!
  !! The fields in the manifest file are separated by the number of axis
  !! dimensions (minus the time dimension).  This type is to facilitate this
  !! separation.
  TYPE manifest_fields_type
     INTEGER :: num_1d = 0 !< Number of 1D fields in fields_1d
     INTEGER :: num_2d = 0 !< Number of 2D fields in fields_2d
     INTEGER :: num_3d = 0 !< Number of 3D fields in fields_3d
     INTEGER :: num_4d = 0 !< Number of 4D fields in fields_4d
     TYPE(manifest_field_type), DIMENSION(:), ALLOCATABLE :: fields_1d !< Array of 1D fields
     TYPE(manifest_field_type), DIMENSION(:), ALLOCATABLE :: fields_2d !< Array of 2D fields
     TYPE(manifest_field_type), DIMENSION(:), ALLOCATABLE :: fields_3d !< Array of 3D fields
     TYPE(manifest_field_type), DIMENSION(:), ALLOCATABLE :: fields_4d !< Array of 4D fields
  END TYPE manifest_fields_type

  PRIVATE
  PUBLIC :: write_diag_manifest

CONTAINS

  ! PUBLIC routines
  !> \brief Public routine that will start the writing of the manifest file.
  !!
  !! This routine is written in such a way that only the root MPI process and the
  !! master OpenMP thread will attempt to write the file.
  SUBROUTINE write_diag_manifest(file)
    INTEGER, INTENT(in) :: file

    INTEGER :: file_unit, ios !< Fortran file unit, and status of file open
    INTEGER :: num_static, num_temporal !< Used to know if any fields are recorded
    INTEGER :: year, month, day, hour, minute, second !< to hold data on current model time.
    TYPE(manifest_fields_type) :: static_fields !< Type to hold all static fields
    TYPE(manifest_fields_type) :: temporal_fields !< Type to hold all non-static fields
    CHARACTER(len=128) :: maniFileName !< Manifest file name
    CHARACTER(len=32) :: filename_appendix !< to hold file name appendix from fms_io
    CHARACTER(len=24) :: start_date !< String to hold init time of diag_manager

    ! Used to determine if the ensemble number.  filename_appendix will contain an
    ! the string ens_ if running with multiple ensembles.  If running only one
    ! ensemble, then filename_appendix will not contain that string.
    CALL get_filename_appendix(filename_appendix)

    ! Get the file name.  Do not need to worry about tiles or ensembles.  Only
    ! writing one manifest file per history file defined in diag_table.
    maniFileName = TRIM(files(file)%name)//".mfst"
    ! prepend the file start date if prepend_date == .TRUE.
    IF ( prepend_date ) THEN
       call get_date(diag_init_time, year, month, day, hour, minute, second)
       write (start_date, '(1I20.4, 2I2.2)') year, month, day

       maniFileName = TRIM(adjustl(start_date))//'.'//TRIM(maniFileName)
    END IF

    ! Extract static and non-static fields data
    static_fields = get_diagnostic_fields(file, static=.TRUE.)
    temporal_fields = get_diagnostic_fields(file, static=.FALSE.)

    ! Get the number of fields to write to manifest file
    ! Need to gather data from all PEs for the component/pelist
    num_static = static_fields%num_1d + static_fields%num_2d + static_fields%num_3d + static_fields%num_4d
    num_temporal = temporal_fields%num_1d + temporal_fields%num_2d + temporal_fields%num_3d + temporal_fields%num_4d

    ! This bulk of this routine should only be called by the rootPE, and only from
    ! ens_01 If running a single ensemble, filename_appendix will not contain the
    ! string ens_

!$OMP MASTER
    IF ( mpp_pe() .EQ. mpp_root_pe() .AND.&
         & (INDEX(filename_appendix,'ens_').EQ.0 .OR. INDEX(filename_appendix,'ens_01').GT.0) ) THEN
       ! Open the file for writing, but only if we have something to write
       IF ( num_static + num_temporal .GT. 0 ) THEN
          ! Get a free Fortran file unit number
          file_unit = get_unit()
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
!$OMP END MASTER
    ! Free up memory used
    CALL destroy_manifest_fields_type(static_fields)
    CALL destroy_manifest_fields_type(temporal_fields)
  END SUBROUTINE write_diag_manifest

  ! PRIVATE routines
  !> \brief De-allocate arrays used in the manifest_fields_type
  SUBROUTINE destroy_manifest_fields_type(manifest_fields)
    TYPE(manifest_fields_type), INTENT(inout) :: manifest_fields

    ! Set all num_?d to 0
    manifest_fields%num_1d = 0
    manifest_fields%num_2d = 0
    manifest_fields%num_3d = 0
    manifest_fields%num_4d = 0
    ! De-allocate the arrays
    IF ( ALLOCATED(manifest_fields%fields_1d) ) DEALLOCATE(manifest_fields%fields_1d)
    IF ( ALLOCATED(manifest_fields%fields_2d) ) DEALLOCATE(manifest_fields%fields_2d)
    IF ( ALLOCATED(manifest_fields%fields_3d) ) DEALLOCATE(manifest_fields%fields_3d)
    IF ( ALLOCATED(manifest_fields%fields_4d) ) DEALLOCATE(manifest_fields%fields_4d)
  END SUBROUTINE destroy_manifest_fields_type

  !> \brief Allow ASSIGNMENT(=) operator to work on TYPE(manifest_field_type)
  !!
  !! Simply assign the type on the rhs to the type on the lhs of the `=`.
  SUBROUTINE manifest_field_type_assign(lhs,rhs)
    TYPE(manifest_field_type), INTENT(out) :: lhs !< lhs, target
    TYPE(manifest_field_type), INTENT(in) :: rhs !< rhs, source

    lhs%output_name = rhs%output_name
    lhs%module_name = rhs%module_name
    lhs%input_name = rhs%input_name
    lhs%time_method = rhs%time_method
    lhs%packing = rhs%packing
    lhs%nDim = rhs%nDim
  END SUBROUTINE manifest_field_type_assign

  !> \brief Write the JSON format of the field object.
  SUBROUTINE write_fields(unit, fields)
    INTEGER, INTENT(in) :: unit !< File unit number.  File should already be opened.
    TYPE(manifest_field_type), DIMENSION(:), INTENT(in) :: fields !< Array of fields to write

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

  !> \brief Write the JSON format of the static/temporal object.
  SUBROUTINE write_manifest(unit, fields, static)
    INTEGER, INTENT(in) :: unit !< File unit number.  File should already be opened.
    TYPE(manifest_fields_type), INTENT(in) :: fields !< All fields to be written to manifest file
    LOGICAL, INTENT(in) :: static !< Indicate if the fields in the fields array
                                  !! are static or non-static fields

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

  !> \brief Extract the diagnostic fields, and collect the information about the
  !! fields.
  TYPE(manifest_fields_type) FUNCTION get_diagnostic_fields(file, static)
    INTEGER, INTENT(in) :: file !< diagnostic file, as defined by diag_manager_mod
    LOGICAL, INTENT(in) :: static !< Indicates if looking for static or non-static
                                  !! fields.  .TRUE. indicates looking only for
                                  !! static files.  .FALSE. indicates looking only
                                  !! for non-static fields.

    INTEGER :: i, j, o
    INTEGER :: istat
    TYPE(manifest_field_type) :: maniField
    CHARACTER(len=128) :: maniFileName
    LOGICAL, DIMENSION(:), ALLOCATABLE :: data_written !< Array to indicate if
                                                       !! field was written to file

    ! manifest file name
    maniFileName = TRIM(files(file)%name)//".mfst"

    ALLOCATE(data_written(mpp_npes()), STAT=istat)
    IF ( istat.NE.0 ) THEN
       CALL error_mesg('diag_manifest_mod::get_diagnostic_fields',&
            & 'Unable to allocate array to determine if field written to file.  No manifest file will be created.',&
            & WARNING)
       ! Set all num_?d to 0, to verify they are set
       get_diagnostic_fields%num_1d = 0
       get_diagnostic_fields%num_2d = 0
       get_diagnostic_fields%num_3d = 0
       get_diagnostic_fields%num_4d = 0
    ELSE
       DO j=1, files(file)%num_fields
          o = files(file)%fields(j) ! Position of this field in output_fields array
          ! Determine if any PE has written file

          ! This is a hack for now.  A future version will use a more elaborate
          ! fix.
          IF ( output_fields(o)%local_output ) THEN
             ! Field is only written for specific regions.  Need to mpp_gather to
             ! know if written on any PE other than root_pe -- as only the root_pe
             ! will write the manifest file
             CALL mpp_gather((/output_fields(o)%written_once/), data_written)
          ELSE
             ! Assuming root_pe was involved in writing of the field --- if written
             data_written = output_fields(o)%written_once
          END IF

          IF ( ANY(data_written) .AND. (static.EQV.output_fields(o)%static) ) THEN
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
    END IF
    ! Clean up allocated arrays
    IF (ALLOCATED(data_written)) DEALLOCATE(data_written)
  END FUNCTION get_diagnostic_fields
END MODULE diag_manifest_mod
