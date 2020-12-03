      IMPLICIT  NONE

c-----------------------------------------------------------------------
c     宣言文
c-----------------------------------------------------------------------
C---- set parameter
      include 'param_grid.inc'

C---- REGION NO TEIGI
      integer iyy,imm,idd,izz,minute,idseq
      INTEGER*2 year(MAX_SCANS)
      INTEGER*2 month(MAX_SCANS),day(MAX_SCANS)
      INTEGER*2 hour(MAX_SCANS),minu(MAX_SCANS)

      character*300 seqfile
      character*13  etxt1
      character*8   etxt2
      character*5   etxt3

      integer i, j, k, ipix
      integer iscan, icount, ncount, nwrite

C--- rain_spc variables
      integer iseqm
      real zlon(NPIX),zlat(NPIX),zitoil(NPIX)
      real zirflg(NPIX),rainfg(NPIX)
      real zsnowprb(NPIX)

      real slon(NPIX, MAX_SCANS),slat(NPIX, MAX_SCANS)
      real sitoil(NPIX, MAX_SCANS),sirflg(NPIX, MAX_SCANS)
      real srainfg(NPIX, MAX_SCANS)
      real ssnowprb(NPIX, MAX_SCANS)
      real sorogflg(NPIX, MAX_SCANS)

C----------- Function Declarations
      integer iargc,n_arguments

c-----------------------------------------------------------------------
c     MAIN
c-----------------------------------------------------------------------
c--- call input file
      n_arguments = iargc()
      if (n_arguments .eq. 1) then
        call getarg(1, seqfile)
      else
        stop
      endif

c--- init
      etxt2 = 'Scan No:'
      etxt3 = 'ipix:'

      iscan=0

c--- open seqfiles
      open(50, file=seqfile, access='SEQUENTIAL', 
     *     form='UNFORMATTED', status='OLD', err=103)

c--- read ncount
      read(50, err=105) ncount
      write(*,*) 'ncount ', ncount

c--- check scan number
      if ( ( ncount.le.0 ).or.( ncount.gt.MAX_SCANS ) ) go to 202

c--- loop start
      do icount=1,ncount
         do ipix=1,NPIX

c--- read physical quantities
            read(50, err=105) zlon(ipix),zlat(ipix),iseqm,zitoil(ipix),
     &                        zirflg(ipix),rainfg(ipix),zsnowprb(ipix)

c--- check longitude
            if (zlon(ipix).gt.-900.0) then
            if ((zlon(ipix).lt.-180.0).or.(zlon(ipix).gt.360.0)) then
               etxt1 = '[Error: lon]'
               write(*,600) etxt1,etxt2,icount,etxt3,ipix,zlon(ipix)
            endif
            endif

            if (zlat(ipix).gt.-900.0) then
            if (abs(zlat(ipix)).gt.90.0) then
               etxt1 = '[Error: lat]'
               write(*,600) etxt1,etxt2,icount,etxt3,ipix,zlat(ipix)
            endif
            endif

c--- check zitoil
            if ((rainfg(ipix).ge.-10.0).and.(rainfg(ipix).lt.0.0).and.
     1          (zitoil(ipix).ge.0)) then
               zitoil(ipix)=rainfg(ipix)
            endif

            if ((rainfg(ipix).lt.-10.0).and.(zitoil(ipix).ge.0)) then
               zitoil(ipix)=UNDEF
            endif

         end do                 !--- ipix

c--- read date
         idseq=iseqm

         call symdhm(iyy,imm,idd,izz,minute,idseq)
c         write(*,*) 'date check ',iyy,imm,idd,izz,minute

c--- kakunou
         iscan=iscan+1

         year(iscan)=iyy
         month(iscan)=imm
         day(iscan)=idd
         hour(iscan)=izz           
         minu(iscan)=minute

         do ipix=1,NPIX
            slon(ipix,iscan)=zlon(ipix)
            slat(ipix,iscan)=zlat(ipix)
            sitoil(ipix,iscan)=zitoil(ipix)
            sirflg(ipix,iscan)=zirflg(ipix)
            srainfg(ipix,iscan)=rainfg(ipix)
            ssnowprb(ipix,iscan)=zsnowprb(ipix)
            sorogflg(ipix,iscan)=UNDEF
         enddo

 9999 enddo                     !--- icount

      close(50)

c--- 
      nwrite=iscan

      do icount=iscan+1,MAX_SCANS
         do ipix=1,NPIX
            slon(ipix,icount)=UNDEF
            slat(ipix,icount)=UNDEF
            sitoil(ipix,icount)=UNDEF
            sirflg(ipix,icount)=UNDEF
            srainfg(ipix,icount)=UNDEF
            ssnowprb(ipix,icount)=UNDEF
            sorogflg(ipix,iscan)=UNDEF
         enddo
         year(icount)=IUNDEF
         month(icount)=IUNDEF
         day(icount)=IUNDEF
         hour(icount)=IUNDEF
         minu(icount)=IUNDEF
      enddo

C------------------------------------------
      write(*,*)'Scan No.(nwrite) : ',nwrite

      do j=1,nwrite
         do i=1,NPIX
            if (srainfg(i,j).gt.0.0) then
               write(*,*) slon(i,j),slat(i,j),
     &         srainfg(i,j),ssnowprb(i,j),sorogflg(i,j)
            endif
         enddo
      enddo

      call exit(0)

C-------- (exit code) ---------------------
 103  write(*,*) " [E] 103: open error (input) "//seqfile
      close(50)
      call exit(103)

 105  write(*,*) " [E] 105: read error (input) "//seqfile
      close(50)
      call exit(105)

 202  write(*,*) " [E] 202: scan number error" 
      close(50)
      call exit(202)

 600  format(A13,1X,A8,1X,I5,3X,A5,1X,I3,5X,F9.4)

      END
      
