
C ... Without MPI: Replace mpiaux.f with this

      MODULE MPI
       INTEGER, PARAMETER :: MPI_STATUS_SIZE = 11
       INTEGER MPI_INTEGER,MPI_COMM_WORLD,
     + MPI_CHARACTER,MPI_LOGICAL,MPI_REAL8,MPI_SUM,MPI_LOR,
     + MPI_REAL8,MPI_INTEGER8
       REAL, EXTERNAL ::  PMPI_WTIME
      END MODULE MPI

      subroutine mpi_bcast(NAME,i,MPI_CHARACTER,j,MPI_COMM_WORLD,IERR)

      write(*,*) 'MPI problems'
      stop
      return
      end

      subroutine mpi_recv(ICON,NBCS,MPI_INTEGER,i,IPRO,
     +        MPI_COMM_WORLD,STATUS,IERR)

      write(*,*) 'MPI problems'
      stop
      return
      end

      subroutine mpi_initialized(joo,IERR)
      logical joo

      joo = .false.
      return
      end

      subroutine mpi_finalize(RC)

      write(*,*) 'MPI problems'
      stop
      return
      end

      subroutine mpi_send(ZZZ,NTOT,MPI_REAL8,JPRO,IP,MPI_COMM_WORLD,
     +        IERR)

      write(*,*) 'MPI problems'
      stop
      return
      end

      subroutine mpi_ssend(ZZZ,NTOT,MPI_REAL8,JPRO,IP,MPI_COMM_WORLD,
     +        IERR)

      write(*,*) 'MPI problems'
      stop
      return
      end

      subroutine mpi_gather(IAA,i,MPI_INTEGER,IA,j,MPI_INTEGER2,k,
     +     MPI_COMM_WORLD,IERR)

      write(*,*) 'MPI problems'
      stop
      return
      end

      subroutine mpi_init(IERR)

      write(*,*) 'MPI problems'
      stop
      return
      end

      subroutine mpi_comm_rank( MPI_COMM_WORLD, myid, ierr )

      write(*,*) 'MPI problems'
      stop
      return
      end

      subroutine mpi_comm_size( MPI_COMM_WORLD, numprocs, ierr )

      write(*,*) 'MPI problems'
      stop
      return
      end

      subroutine mpi_abort(MPI_COMM_WORLD,ERRORCODE,IERR)

      write(*,*) 'MPI problems'
      stop
      return
      end

      subroutine mpi_wtime()

      write(*,*) 'MPI problems'
      stop
      return
      end

      subroutine mpi_sendrecv_replace(ZZZ,ILNG,
     +     MPI_REAL8,JPRO2,JP,JPRO,IGP,MPI_COMM_WORLD,STATUS,IERR)

      write(*,*) 'MPI problems'
      stop
      return
      end

      subroutine mpi_wtick()

      write(*,*) 'MPI problems'
      stop
      return
      end

      subroutine mpi_attr_get(MPI_COMM_WORLD,MPI_TAG_UB,ITGI,FLAGI,
     +      IERROR)

      write(*,*) 'MPI problems'
      stop
      return
      end

      subroutine mpi_reduce(CLS, CLSM,ntot,MPI_REAL8,MPI_SUM,IPRO,
     +         MPI_COMM_WORLD,IERR)

      write(*,*) 'MPI problems'
      stop
      return
      end

      subroutine mpi_allreduce(CLS, CLSM,ntot,MPI_REAL8,MPI_SUM,
     +         MPI_COMM_WORLD,IERR)

      write(*,*) 'MPI problems'
      stop
      return
      end

      subroutine mpi_barrier(MPI_COMM_WORLD,IERR)

      write(*,*) 'MPI problems'
      stop
      return
      end

      subroutine  PMPI_WTIME()

      write(*,*) 'MPI problems'
      stop
      return
      end
