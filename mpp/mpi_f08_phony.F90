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

! Phony mpi_f08 module which provides type names and null values that are needed in the non-MPI build.

module mpi_f08_phony
#ifndef use_libMPI
  implicit none

  type :: mpi_comm
    integer :: mpi_val = -1
  end type mpi_comm

  type :: mpi_info
    integer :: mpi_val = -1
  end type mpi_info

  type :: mpi_group
  end type mpi_group

  type :: mpi_request
  end type mpi_request

  type :: mpi_datatype
  end type mpi_datatype

  type(mpi_comm), parameter :: MPI_COMM_NULL = mpi_comm()
  type(mpi_info), parameter :: MPI_INFO_NULL = mpi_info()
  type(mpi_group), parameter :: MPI_GROUP_NULL = mpi_group()
  type(mpi_request), parameter :: MPI_REQUEST_NULL = mpi_request()
  type(mpi_datatype), parameter :: MPI_DATATYPE_NULL = mpi_datatype()

  type(mpi_datatype), parameter :: MPI_REAL4 = mpi_datatype()
  type(mpi_datatype), parameter :: MPI_REAL8 = mpi_datatype()
#endif
end module mpi_f08_phony
