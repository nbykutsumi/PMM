c      subroutine INGEST_GSMaP_AMSU(fname, first_line, nread, nfound,
c     $ size_A, size_B, year, day, UTC, startDate, endDate,
c     $ latB, lonB, Sfc, SI, LZ, QC, RR)

      IMPLICIT NONE

      integer    AMSU_PIX
      integer    AMSU_MAX_SCANS
      parameter (AMSU_PIX = 96)
      parameter (AMSU_MAX_SCANS = 9000)

      integer mtin
      parameter   (mtin = 10)

      CHARACTER*300 fname

      INTEGER first_line, nread, nfound, status

      INTEGER size_A, size_B
      INTEGER year(AMSU_MAX_SCANS)
      INTEGER day(AMSU_MAX_SCANS)
      INTEGER UTC(AMSU_MAX_SCANS)
      INTEGER startDate(3)
      INTEGER endDate(3)
      REAL latB(AMSU_PIX,AMSU_MAX_SCANS)
      REAL lonB(AMSU_PIX,AMSU_MAX_SCANS)
      REAL Sfc(AMSU_PIX,AMSU_MAX_SCANS)
      REAL SI(AMSU_PIX,AMSU_MAX_SCANS)
      REAL LZ(AMSU_PIX,AMSU_MAX_SCANS)
      REAL RR(AMSU_PIX,AMSU_MAX_SCANS)
      REAL SNOW(AMSU_PIX,AMSU_MAX_SCANS) 
      REAL FLZ(MWS_PIX,MWS_MAX_SCANS)

      INTEGER OROG(MWS_PIX,MWS_MAX_SCANS)
      INTEGER QC(AMSU_PIX,AMSU_MAX_SCANS)
      INTEGER SnowM(AMSU_PIX,AMSU_MAX_SCANS)

      INTEGER i,j,ipix,ict
      integer iargc,n_arguments

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      n_arguments = iargc()
      if (n_arguments .eq. 1) then
        call getarg(1, fname)
      else
        stop
      endif

      write(*,*)'Opening File : '//fname

      Open(mtin,file=fname,FORM='UNFORMATTED',
     $     ACCESS='DIRECT',RECL=4,err=103)

      read(mtin,rec=1,err=105)size_A
      write(*,*)size_A
      read(mtin,rec=2,err=105)size_B
      write(*,*)size_B

      ict=3
      do i=1,size_B
        read(mtin,rec=ict,err=105)year(i)
        ict=ict+1
      enddo
      do i=1,size_B
        read(mtin,rec=ict,err=105)day(i)
        ict=ict+1
      enddo
      do i=1,size_B
        read(mtin,rec=ict,err=105)UTC(i)
        ict=ict+1
      enddo

      do i=1,3
        read(mtin,rec=ict,err=105)startDate(i)
        ict=ict+1
      enddo
      do i=1,3
        read(mtin,rec=ict,err=105)endDate(i)
        ict=ict+1
      enddo

      do j=1,size_B
        do i=1,AMSU_PIX
          read(mtin,rec=ict,err=105)latB(i,j)
          ict=ict+1
        enddo
      enddo
      do j=1,size_B
        do i=1,AMSU_PIX
          read(mtin,rec=ict,err=105)lonB(i,j)
          ict=ict+1
        enddo
      enddo
      do j=1,size_B
        do i=1,AMSU_PIX
          read(mtin,rec=ict,err=105)Sfc(i,j)
          ict=ict+1
        enddo
      enddo
      do j=1,size_B
        do i=1,AMSU_PIX
          read(mtin,rec=ict,err=105)SI(i,j)
          ict=ict+1
        enddo
      enddo
      do j=1,size_B
        do i=1,AMSU_PIX
          read(mtin,rec=ict,err=105)LZ(i,j)
          ict=ict+1
        enddo
      enddo
      do j=1,size_B
        do i=1,AMSU_PIX
          read(mtin,rec=ict,err=105)QC(i,j)
          ict=ict+1
        enddo
      enddo
      do j=1,size_B
        do i=1,AMSU_PIX
          read(mtin,rec=ict,err=105)SnowM(i,j)
          ict=ict+1
        enddo
      enddo
      do j=1,size_B
        do i=1,AMSU_PIX
          read(mtin,rec=ict,err=105)RR(i,j)
          ict=ict+1
        enddo
      enddo
      do j=1,size_B
        do i=1,MWS_PIX
          read(mtin,rec=ict,err=105)FLZ(i,j)
          ict=ict+1
        enddo
      enddo
      do j=1,size_B
        do i=1,MWS_PIX
          read(mtin,rec=ict,err=105)OROG(i,j)
          ict=ict+1
        enddo
      enddo

      close(mtin)

C----------------------------------------
      do j=1,size_B
        do i=1,AMSU_PIX
           if (RR(i,j).gt.0.0) then
              write(*,*) lonB(i,j),latB(i,j),RR(i,j)
           endif
        enddo
      enddo

C-------- (exit code) ---------------------
      call exit(0)

 103  write(*,*) " [E] 103: open error (input) "//fname
      close(mtin)
      call exit(103)

 105  write(*,*) " [E] 105: read error (input) "//fname
      close(mtin)
      call exit(105)

      end

