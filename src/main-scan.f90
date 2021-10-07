!
!  mstm main program
!
!
!  original release: 15 January 2011
!  21 February 2011: modifications to fixed orientation efficiency factor
!  November 2012: multiple run input file options, and other housekeeping
!  February 2013: mpi modifications, and more housekeeping.
!
      program main
      use mpidefs
      use mpidata
      use intrinsics
      use spheredata
      use numconstants
      use specialfuncs
      use miecoefdata
      use translation
      use solver
      use scatprops
      use nearfield
      use interpolate
      implicit none
      integer :: nsphere,neqns,nodrmax,nodrt,i,k,niter,istat,numtheta, &
                 nblkt,nodrg,m,n,w,fixedorrandom,numargs,calctmatrix,maxiter,&
                 nodrta(1),calcnf,calcamn,ma,na,nsend,nfplane,nfoutunit,&
                 nfoutdata,maxmbperproc,trackiterations,nonactive,normalizesm, &
                 storetranmat,niterstep,j,smnumprocs,runcomm,runnumprocs, &
                 excitedsphere,runnum,appendafile,appendnffile,noff,nblk, &
                 writespheredata,incidentortargetframe,amnunit,rank, &
                 printinputdata,runprintunit,numprocs,fseekstatus,&
                 amnfileposition,azimuthaverage
      integer, allocatable :: nodr(:),ntran(:),sphereblk(:),sphereoff(:), &
                              hostsphere(:),numberfieldexp(:)
      logical :: moreruns,fftranpresent
      logical, allocatable :: tmonfile(:)
      real (8) :: alphadeg,betadeg,alpha,beta,epsmie,epstran,epssoln, &
                  qexttot,qabstot,xv,qscatot,asymparm, &
                  phideg,theta1d,theta2d,thetad,costheta,sintheta,phi, &
                  time1,time2,epstcon,scatrotvec(3),&
                  cbeam,gbfocus(3),maxerr,nfplanepos,nfplanevert(2,2), &
                  deltax,gammadeg,epspw,gamma,qexttotpar,qexttotper, &
                  qabstotpar,qabstotper,qscatotpar,qscatotper,cphi,sphi,s11, &
                  nfdistance,rposi(3),dummy(4),qeff(4),xgeff
      real(8), allocatable :: xsp(:), rpos(:,:),qext(:,:),qabs(:,:), &
               qsca(:,:),smc(:,:,:),smt(:,:,:),qabsvol(:,:),smccf(:,:,:), &
               smtcf(:,:,:),smticf(:,:,:),scatteringanglearray(:,:), &
               s00(:,:,:),s02(:,:,:),sp22(:,:,:),sm22(:,:,:)
      complex(8) :: sa(4),rimedium(2)
      complex(8), allocatable :: amnp(:),amnp0(:,:,:,:),ri(:,:), &
                  gmn(:),amnp1(:,:,:),amnp2(:,:,:),pmnp0(:,:,:,:)
      character*60 :: inputfile,outfile,tmatrixfile,&
                      amnfile,nfoutfile,oldamnfile
      character*60, allocatable :: tmfile(:)

      !CWH Added 4-7-2017
      !Added in chronological order...so don't make much sense and will
      !reorder when things seem complete
      !Variables for multiple wavelengths and interpolation
      integer :: nwav, nrefwav, idvout, ilam, iwav, dpcalc, mn, nn1
      integer :: nfnewcalc, nfcalcindex, neqnsnf, writescatdata
      real(8) :: xdp0(3), rdp0(3), apol(3)
      real(8), allocatable :: refmwav(:,:), wavlist(:),&
                             xdp(:,:), rdp(:,:), kmed(:)
      complex(8) :: dpmom(3)
      !For medium refractive index and wavenumber
      real(8) ::  xvn
      !wavenumbers for scattering calculation and near-field calculation
      complex(8) :: medk
      complex(8) :: knf
      complex(8), allocatable :: epfield(:,:)
      !refind will be the refractive index of each particle at each
      !wavelength with the 0-th element being the medium
      !riw will replace ri in things like the Mie calculation, might not
      !need both..we'll see
      !I define a bunch of extra variables to fill at each wavelength
      !for some reason
      complex(8), allocatable :: refind(:,:,:), riw(:,:)
      real(8),allocatable :: rposn(:,:), xspn(:)
      real(8), allocatable :: xspnf(:), nodrnf(:)
      !integer :: nfnewcalc, nfcalcindex
      complex(8), allocatable :: amnpnf(:)
      real(8) :: nfint
      real(8), allocatable :: xaccept0(:), xaccept(:)
      character(60) :: epfile, acceptfile
      complex(8) :: eincacc(3),hincacc(3),escatacc(3),hscatacc(3),etotacc(3),htotacc(3)
      complex(8) :: eincdon(3),hincdon(3),escatdon(3),hscatdon(3),etotdon(3),htotdon(3)
      real(8) :: incaccint, scataccint, totaccint, xgpmax(3)
      real(8) :: acceptxyz(3), accorient(3)
      complex(8) :: incaccint2, scataccint2, totaccint2
      real(8) :: dummyD,pi2,eps0,mu0
      !data pi2/3.141592653589793d0/, eps0/8.85418781762039d-12/
      !YJ added 4-1-2019
      complex(8), allocatable :: amnp3(:,:,:)
      complex(8), allocatable :: amnpt(:,:)
      complex(8), allocatable :: pmnp(:)
      complex(8), allocatable :: pmnp1(:,:,:),pmnp2(:,:,:),pmnp3(:,:,:)
      complex(8), allocatable :: pmnpt(:,:)
      integer :: b11, b12
      integer, allocatable :: nblk_nf(:),noff_nf(:)
      real(8) :: dpmommag, omega, lambda
      real(8) :: powerdecay, powerself, powerscat ! decay power components
      real(8) :: ratedecay, rateself, ratescat ! decay power components
      real(8) :: totqy, vacqy, intqy
      real(8) :: totlifetime, vaclifetime
      real(16) :: abspower, radpower ! absoprtion by sphere, radiated energy, and sum
      real(8) :: c, cmicro
      real(8) :: hbar, dpauconv
      data c/299792458/, cmicro/299792458d6/
      !YJ added 2-27-2020
      data hbar/1.054571817d-34/, dpauconv/8.4783536255d-30/
      !data mu0/1.2566370614d-6/
      !YJ added END
 
      character*60 jobname, unpolscat, parascat, perpscat, commorigin
      character*25 dummyCh

      pi2=4.d0*atan(1.d0)
      mu0=4.d0*pi2*1.d-7
      eps0=1.d0/mu0/c**2

      !Mackowski's default is 1
      printinputdata=0
      !printinputdata=0
      amnunit=20

! initialize the mpi environment
!
      call mstm_mpi(mpi_command='init')
      call mstm_mpi(mpi_command='barrier')
      call mstm_mpi(mpi_command='size',mpi_size=numprocs)
      call mstm_mpi(mpi_command='rank',mpi_rank=rank)
      call mstm_mpi(mpi_command='barrier')
!
!  the loop through input file runs
!
!  Removing Mackowski's loop through input filename.  I will add the
!  loop for wavelengths and other variations can be done with python
!  scripting or something...I guess
!  He uses the runnum variable for some useful things, so I may keep
!  that, we'll see.

!      moreruns=.true.
!      do while(moreruns)
!         runnum=runnum+1
!
!  command line argument retrieval for input file
!
      runnum=1
      do i=0,numprocs-1
         if(i.eq.rank) then
            numargs=mstm_nargs()
            if(numargs.eq.0) then
               inputfile='mstm.inp'
            else
               call mstm_getarg(inputfile)
            endif
            call inputdata(inputfile,printinputdata)
         endif
         call mstm_mpi(mpi_command='barrier')
      enddo
!
!  reading of run and sphere data, setting up of arrays
!
      call getspheredata(number_spheres=nsphere)
      if(allocated(xsp)) deallocate(xsp,rpos,nodr,ntran,ri,sphereblk, &
               sphereoff,hostsphere,numberfieldexp,tmfile,tmonfile)
      allocate(xsp(nsphere),rpos(3,nsphere),nodr(nsphere),ntran(nsphere), &
               ri(2,0:nsphere),sphereblk(nsphere),sphereoff(nsphere+1), &
               hostsphere(nsphere),numberfieldexp(nsphere),tmfile(nsphere), &
               tmonfile(nsphere))
      call getspheredata(sphere_size_parameters=xsp,sphere_positions=rpos, &
           sphere_refractive_indices=ri,volume_size_parameter=xv, &
           host_spheres=hostsphere,tmatrix_file=tmfile,tmatrix_on_file=tmonfile)
      call getrunparameters(mie_epsilon=epsmie,translation_epsilon=epstran, &
           solution_epsilon=epssoln,max_number_iterations=niter, &
           fixed_or_random_orientation=fixedorrandom,output_file=outfile, &
           min_scattering_angle_deg=theta1d,max_scattering_angle_deg=theta2d, &
           number_scattering_angles=numtheta,gaussian_beam_constant=cbeam, &
           gaussian_beam_focal_point=gbfocus,run_print_unit=runprintunit, &
           max_memory_per_processor=maxmbperproc, &
           normalize_scattering_matrix=normalizesm, &
           store_translation_matrix=storetranmat, &
           near_field_distance=nfdistance, &
           iterations_per_correction=niterstep, &
           write_sphere_data=writespheredata)


      !Following Mackowski's format, I store the input variables as
      !private variables in the module spheredata, then extract them
      !with a function getscandata, where everything is an optional
      !variable, so you can use it to get anything as needed.  I get the
      !variables that determine array sizes and then initialize them,
      !then call it again
      call getscandata(nwav=nwav, dpcalc=dpcalc, &
            write_scat_data = writescatdata,  accorient2=accorient)
      !call getscandata(nwav=nwav)
      allocate(wavlist(nwav),refind(2,0:nsphere,nwav),kmed(nwav), &
               rposn(3, nsphere), xspn(nsphere))
      call getscandata(wavlist=wavlist, kmed=kmed, refind=refind)
      if(dpcalc.ne.0) then
          allocate(xdp(3,nsphere), rdp(3,nsphere),epfield(nwav,3))
          cbeam=-1.0
          allocate(xaccept0(3))
          allocate(xaccept(3))
          call getscandata(xdp0=xdp0,rdp0=rdp0, xdp=xdp, rdp=rdp,&
                      epfield=epfield, xaccept0=xaccept0, acceptfile=acceptfile)
          open(8,file=acceptfile,status='replace',action='write')
          close(8)
!  Going to mirror print options from nfoutdata for the dipole field at
!  the acceptor
!  nfoutdata = 0 -> print intensity          
!  nfoutdata = 1 -> print E-field components
!  nfoutdata = 2 -> print E-field and H-field components
          call getrunparameters(near_field_output_data=nfoutdata)
      endif 
      !endif
      !calculation would start with fresh amn file...otherwise the
      !append would do exactly that to an old calculation
      call getrunparameters(scattering_coefficient_file=amnfile)

      appendafile = 0
      if(amnfile.ne.' ') then
          !appendafile = 1
         open(3,file=amnfile)
         close(3,status='delete')
      endif
      if(outfile.ne.' ') then
         open(3,file=outfile)
         close(3,status='delete')
         do i=1,len_trim(outfile)
            if(outfile(i:i+3).eq.'.out') then
                jobname = outfile(1:(i-1))
                exit
            else
                jobname=outfile
            endif
         enddo

         unpolscat = trim(jobname) // ".unpolscat"
         parascat = trim(jobname) // ".parascat"
         perpscat = trim(jobname) // ".perpscat"
         commorigin = trim(jobname) // ".commorigin"
         open(20,file=unpolscat)
         close(20,status='delete')
         open(21,file=parascat)
         close(21,status='delete')
         open(22,file=perpscat)
         close(22,status='delete')
         open(23,file=commorigin)
         close(23,status='delete')
      endif
      appendafile = 0
      !Mackowski now has an rimedium.  He did not in v2.2.  Therefore, I
      !should be able to get rid of medri
      rimedium(:)=ri(:,0)

      !Interpolating the refractive index using the same method from DDA
      !It fits two data sets, the real and imaginary parts.  As a
      !result, it needs x, y1, y2 where y1 = real and y2 = imaginary
      !The arguments are x_i, y1(x_i), y2(x_i), array of x, array of y1,
      !array of y2.  x, y1, y2 are the data from Johnson and Christy,
      !etc.
      !Notes: To add multiple materials you can...
      !Separate interpolation from filling refind()
      !Add loop over ncomponent, make refindre and refindim arrays of
      !length ncomponent.  Add dimension to refmwav so you you pull out
      !the material you want.  Then fill refind using some index
      !icomponent that tells you what material each sphere is.  Would be
      !easiest to add an integer to the end of the geometry
      !specification, probably      

      !Should broadcast appropriately.  Still feels strange doing it in
      !Fortran so the explicit loop is just commented ou
      !refind(1,:,:) = refind(1,:,:)/rimedium(1)
      !refind(2,:,:) = refind(2,:,:)/rimedium(2)
      !!do i=1,nwav      
      !!  do k=1,nsphere
      !!      refind(1,k,i) = refind(1,k,i)/rimedium(1)
      !!      refind(2,k,i) = refind(2,k,i)/rimedium(2)
      !!  enddo
      !!enddo
      !refind(1:2,0,:) is the wavelength-independent refractive index of the
      !medium
      !do i=1,nwav
      !  refind(1,0,i) = rimedium(1)
      !  refind(2,0,i) = rimedium(2)
      !enddo



!
!  determine if optical activity is present
!
      nonactive=1
      do i=1,nsphere
         !if(cdabs(ri(1,i)-ri(2,i)).gt.1.d-10) then
         if(cdabs(refind(1,i,1)-refind(2,i,1)).gt.1.d-10) then
            nonactive=0
            exit
         endif
      enddo

!
!  Check for nf calc option
!  CWH 4-7-2017  To accomodate the wavelength scanning with near-field
!  plots, I added a wavelength specification for the near-field plot in
!  the input.  The calculation should automatically select the
!  scattering calculation at the wavelength closest to the one you
!  specified and generate the near-field for that.
!
      call getrunparameters(calculate_near_field=calcnf)
      call getscandata(nfcalcindex=nfcalcindex)
      !if(calcnf.eq.1) then
      !    call getscandata(nfcalcindex=nfcalcindex)
      !endif

!  Raman polarizability
!  Dipole moment p(:) = apol(:)*efield(:)
      apol(1) = 1.0
      apol(2) = 1.0
      apol(3) = 1.0
!  In my dipole code I just deleted everything for fixedorrandom=1
!  because we do not need random orientation calculations, again, hope
!  they still work, too
!  Loop over wavelengths

!  There is a chance I broke the random orientation calculation.  I only
!  use the fixed orientation, so if running random orientation, may want
!  to double check everything


      accorient(:) = accorient(:) * dpauconv
      do ilam=1,nwav
         !Define ri for given wavelength
         do k=0,nsphere
           ri(1,k)=refind(1,k,ilam)
           ri(2,k)=refind(2,k,ilam)
         !write(*,*) "refractive indices", k, ri(1,k), ri(2,k)
         enddo
         !Dipole moment
         if(dpcalc.ne.0) then
            do k=1,3 
              !dpmom(k) = epfield(ilam,k)*apol(k)
              dpmom(k) = epfield(ilam,k)*apol(k)*dpauconv
            enddo
         else
            do k=1,3 
              dpmom(k) = cmplx(0.0,0.0)
            enddo
         endif

!  Define wavenumber for the given wavelength
!  Because it seemed like a good idea, kmed(:) doesn't use the
!  refractive index of the medium and is simply 2pi/lamda
!  I guess this is good if your medium is optically active and you ant
!  to separate wavenumbers for L- and R-, but this all still may need
!  some tweaking
!  For now just define medk using the rimedium(1)
!  In things like the size parameter, which Mackowski defines as being
!  real, we can't just use the complex refractive index for the medium.
         !medk = dble(rimedium(1))*kmed(ilam)
         !write(*,*) "rimedium", rimedium
         !medk = rimedium(1)*kmed(ilam)

         medk = kmed(ilam)

         !Define k*r radii and cartesian coordinates.  These are used
         !everywhere except in nearfieldgridcalc where both r and kr are
         !required to produce the correct map
         !If you find a place with it not changed, it's probably wrong
         !xvn        = dble(rimedium(1))*dble(medk)*xv
         rposn(:,:) = dble(medk)*rpos(:,:)
         xspn(:)  = dble(medk)*xsp(:)
         !Modified size parameter
         xvn = dble(medk)*xv

         !I never do anything with scattering matrices, but I do hope to
         !not mess up any of those things
         if(numtheta.gt.0) then
            if(allocated(smt)) deallocate(smt,smtcf,smticf)
            allocate(smt(4,4,numtheta),smtcf(4,4,numtheta),smticf(4,4,numtheta))
            if(allocated(scatteringanglearray)) deallocate(scatteringanglearray)
            allocate(scatteringanglearray(2,numtheta))
            call getrunparameters(scattering_angle_array=scatteringanglearray)
         endif
         if(cbeam.eq.0.d0) then
            xgeff=xvn
         else
               xgeff=1.d0/cbeam
         endif
!
!  host configurations
!
         call hostconfiguration(nsphere,hostsphere,numberfieldexp)
!
!  calculation of sphere mie coefficients, order limits
!
         call miecoefcalc(nsphere,xspn,ri, & 
                            hostsphere,numberfieldexp,epsmie, &
              tmatrix_file=tmfile)
         call getmiedata(sphere_order=nodr,max_order=nodrmax,number_equations=neqns)
!
!  determine orders required to expand scattered fields about target origin
!
         call tranorders(nsphere,nodr,rposn,ri,hostsphere,epstran,ntran,nodrt)
!
!  determine the size of the parallel run and set it up
!
!  CWH:  I do not see where runnumprocs is defiend, worst case, it looks
!  like you should be able to use sm_number_processors if you move this
!  about 20 lines down
         call mstm_mpi(mpi_command='barrier')
         call mpisetup(nsphere,runnumprocs,runcomm)
         if(fixedorrandom.le.1) then
            call rottranmtrxsetup(nsphere,nodr,rposn,ri,storetranmat,&
                 nfdistance,ntran,runprintunit,fftranpresent,runcomm)
         endif
!
!  report the size of the run
!
      !  if(rank.eq.0) then
      !     write(runprintunit,'('' maximum sphere order:'',i5)') nodrmax
      !     write(runprintunit,'('' estimated T matrix order:'',i5)') nodrt
      !     write(runprintunit,'('' number of equations:'',i9)') neqns
      !     call flush(runprintunit)
      !  endif!
!
         if(fixedorrandom.eq.1) then
!
!  random orientation option
!
            call getrunparameters(calculate_t_matrix=calctmatrix,&
                 t_matrix_file=tmatrixfile, &
                 t_matrix_convergence_epsilon=epstcon, &
                 excited_sphere=excitedsphere, &
                 sm_number_processors=smnumprocs)
            if(allocated(qext)) deallocate(qext,qabs,qsca,qabsvol)
            allocate(qext(nsphere,1), qabs(nsphere,1), qsca(nsphere,1),qabsvol(nsphere,1))
!
!  this option calculates the T matrix either from the beginning or where left off
!
            if(calctmatrix.ge.1) then
               if(rank.eq.0) time1=mytime()
               call tmatrixsoln(neqns,nsphere,nodr,ntran,nodrt,xspn,rposn,hostsphere,numberfieldexp, &
                    rimedium,epssoln,epstcon,niter,calctmatrix,tmatrixfile,fftranpresent, &
                    niterstep,excitedsphere,qext,qabs,qsca,istat,mpi_comm=runcomm)
               call mstm_mpi(mpi_command='barrier')
               nodrta(1)=nodrt
               call mstm_mpi(mpi_command='bcast',mpi_number=1,mpi_send_buf_i=nodrta, &
                    mpi_rank=0)
               nodrt=nodrta(1)
               if(rank.eq.0) then
                  time2=mytime()-time1
                  write(dummyCh,'(I5,A)') ilam,' Execution time:'
                  call timewrite(runprintunit,dummyCh,time2)
               endif
               call rottranmtrxclear()
            else
!
!  and this has the T matrix already calculated and stored in the file.
!
!  read the order of the T matrix and broadcast to the processors.
!
               if(rank.eq.0) then
                  open(3,file=tmatrixfile)
                  read(3,*) i,nodrt,i,i
                  close(3)
                  write(runprintunit,'('' t matrix order:'',i5)') nodrt
                  call flush(runprintunit)
               endif
               nodrta(1)=nodrt
               call mstm_mpi(mpi_command='bcast',mpi_send_buf_i=nodrta, &
                  mpi_number=1,mpi_rank=0)
               nodrt=nodrta(1)
               call mstm_mpi(mpi_command='barrier')
            endif
!
!  the T matrix is available; calculate the random orientation scattering matrix
!
            nblkt=nodrt*(nodrt+2)
            if(numtheta.eq.0) then
               nodrg=2
            else
               nodrg=nodrt*2
            endif
            if(allocated(smc)) deallocate(smc,smccf)
            allocate(smc(4,4,0:nodrg),smccf(4,4,0:nodrg))
            call ranorientscatmatrix(xgeff,nsphere,nodrt,nodrg,cbeam,tmatrixfile, &
                 smc,smccf,qext,qabs,qsca, &
                 number_processors=smnumprocs)
            if(rank.eq.0) then
               qexttot=0.
               qabstot=0.
               do i=1,nsphere
                  if(hostsphere(i).eq.0) then
                     qexttot=qexttot+qext(i,1)*xspn(i)*xspn(i)/xgeff/xgeff
                     qabstot=qabstot+qabs(i,1)*xspn(i)*xspn(i)/xgeff/xgeff
                  endif
               enddo
               qscatot=qexttot-qabstot
               asymparm=dble(smc(1,1,1)/smc(1,1,0))/3.d0
               do i=1,numtheta
                  costheta=cos(scatteringanglearray(1,i)*pi/180.d0)
                  call ranorienscatmatrixcalc(costheta,smc,nodrg,smt(:,:,i))
                  call ranorienscatmatrixcalc(costheta,smccf,nodrg,smtcf(:,:,i))
                  smticf(:,:,i)=smt(:,:,i)-smtcf(:,:,i)
               enddo
            endif
!
! volume absorption and correction
!
            do i=1,nsphere
               qabsvol(i,:)=qabs(i,:)
               do j=1,nsphere
                  if(hostsphere(j).eq.i) then
                     qabsvol(i,:)=qabsvol(i,:)-qabs(j,:)*xspn(j)*xspn(j)/xsp(i)/xsp(i)
                  endif
               enddo
            enddo
            do i=1,nsphere
               if(.not.tmonfile(i)) then
                  if(dimag(ri(1,i)).eq.0.d0.and.dimag(ri(2,i)).eq.0.d0) then
                     qabsvol(i,:)=0.d0
                  endif
               endif
            enddo
            qabs=0.d0
            do i=1,nsphere
               qabs(i,:)=qabsvol(i,:)
               do j=1,nsphere
                  if(hostsphere(j).eq.i) then
                     qabs(i,:)=qabs(i,:)+qabsvol(j,:)*xspn(j)*xspn(j)/xsp(i)/xsp(i)
                  endif
               enddo
            enddo
         elseif(fixedorrandom.eq.0) then
!
!  fixed orientation option
!
            call mstm_mpi(mpi_command='barrier')
            if(runnum.eq.1) oldamnfile=' '
            call getrunparameters(calculate_scattering_coefficients=calcamn, &
                    scattering_coefficient_file=amnfile, &
                    incident_azimuth_angle_deg=alphadeg, &
                    incident_polar_angle_deg=betadeg, &
                    track_iterations=trackiterations, &
                    !append_a_file=appendafile, &
                    incident_or_target_frame=incidentortargetframe, &
                    azimuth_average_scattering_matrix=azimuthaverage)
            alpha=alphadeg*pi/180.d0
            beta=betadeg*pi/180.d0
            phideg=0.d0
            if(calcamn.le.1) then
               if(allocated(amnp)) deallocate(amnp,pmnp) !deallocate(amnp)
               if(allocated(qext)) deallocate(qext,qabs,qsca,qabsvol)
               !allocate(amnp(neqns*2),qext(nsphere,3), qabs(nsphere,3), &
               !   qsca(nsphere,3),qabsvol(nsphere,3))
               allocate(amnp(neqns*2),pmnp(neqns*2),qext(nsphere,3),&
                  qabs(nsphere,3),qsca(nsphere,3),qabsvol(nsphere,3))
            endif
            if(calcamn.eq.1) then
!
!  this option calculates the scattering coefficients
!
               if(rank.eq.0) time1=mytime()
               call mstm_mpi(mpi_command='barrier')
               !TEST
               !write(*,*) "fftranpresent: ", fftranpresent
               !write(*,*) "niter: ", niter
               !TEST_END
               !call fixedorsoln(neqns,nsphere,nodr,alpha,beta,cbeam,xspn,rposn,hostsphere,&
               !    numberfieldexp,rimedium,epssoln,epstran,niter,amnp,qext,qabs,qsca, &
               !    maxerr,maxiter,trackiterations,fftranpresent,niterstep,istat, &
               !    medk, refind(:,:,ilam), rdp, dpmom, mpi_comm=runcomm)
               call fixedorsoln(neqns,nsphere,nodr,alpha,beta,cbeam,xspn,rposn,hostsphere,&
                   numberfieldexp,rimedium,epssoln,epstran,niter,amnp,pmnp,qext,qabs,qsca, &
                   maxerr,maxiter,trackiterations,fftranpresent,niterstep,istat, &
                   medk, refind(:,:,ilam), rdp, dpmom, mpi_comm=runcomm)
               !TEST
               !write(*,*) "amnp after fixedorsoln: ", amnp
               !write(*,*) "pmnp after fixedorsoln: ", pmnp
               !TEST_END
            !Save scattering coefficients for nea-field calculation
            !if(calcnf.eq.1.and.ilam.eq.nfcalcindex) then
            if(ilam.eq.nfcalcindex) then
                allocate(amnpnf(neqns*2))
                allocate(xspnf(nsphere))
                allocate(nodrnf(nsphere))
                neqnsnf = neqns*2
                nodrnf(:) = nodr(:)
                xspnf(:) = xspn(:)
                amnpnf(:) = amnp(:)
            endif

!
!  broadcast the scattering coefficients to the other processors
!
               nsend=neqns*2
               call mstm_mpi(mpi_command='barrier')
               call mstm_mpi(mpi_command='bcast',mpi_send_buf_dc=amnp,mpi_number=nsend,&
                    mpi_rank=0)
!
!  write the scattering coefficients to the file
!
               if(rank.eq.0) then
                  time2=mytime()-time1
                  write(runprintunit,'('' max iterations, soln &
                              &   error:'',i6,e13.5)')  maxiter,maxerr
                  write(dummyCh,'(I5,A)') ilam,' Execution time:'
                  call timewrite(runprintunit,dummyCh,time2)
                  if(amnfile.ne.' ') then
                     if(appendafile.eq.0) then
                        open(3,file=amnfile,status='replace', &
                            action='write')
                     else
                        open(3,file=amnfile,position='append')
                     endif
                     !Mackowski's statement
                     !write(3,'(2i9,6e13.5)') nsphere,neqns,alpha,beta,cbeam,rimedium(1),wavlist(ilam) 
                     !Mine
                     write(3,'(2i9,3e13.5)') nsphere,neqns, wavlist(ilam), rimedium(1)
                     noff=0
                     do i=1,nsphere
                        write(3,'(2i5,8e13.5)') nodr(i),numberfieldexp(i),xsp(i), &
                             rpos(:,i),ri(:,i)
                        write(3,'(9e13.5)') qext(i,:),qabs(i,:),qsca(i,:)
                        allocate(amnp1(0:nodr(i)+1,nodr(i),2),amnp2(0:nodr(i)+1,nodr(i),2))
                        nblk=2*nodr(i)*(nodr(i)+2)
                        do k=1,numberfieldexp(i)
                           amnp1=reshape(amnp(noff+1:noff+nblk), &
                                        (/nodr(i)+2,nodr(i),2/))
                           amnp2=reshape(amnp(noff+nblk+1:noff+2*nblk), &
                                        (/nodr(i)+2,nodr(i),2/))
                           do n=1,nodr(i)
                              do m=-n,n
                                 if(m.le.-1) then
                                    ma=n+1
                                    na=-m
                                 else
                                    ma=m
                                    na=n
                                 endif
                                 write(3,'(4e17.9)') amnp1(ma,na,1),amnp2(ma,na,1)
                                 write(3,'(4e17.9)') amnp1(ma,na,2),amnp2(ma,na,2)
                              enddo
                           enddo
                           noff=noff+2*nblk
                        enddo
                        deallocate(amnp1,amnp2)
                     enddo
                     close(3)
                  endif
               endif
            elseif(calcamn.le.0.and.amnfile.ne.' ') then
!
!  this option reads the scattering coefficients from the file
!
               do j=0,numprocs-1
                  if(rank.eq.j) then
                     open(amnunit,file=amnfile,action='read')
                     if(runnum.eq.1.or.oldamnfile.ne.amnfile) then
                        amnfileposition=0
                     endif
                     call mstm_fseek(amnunit,amnfileposition,0,fseekstatus)
                     read(amnunit,'(2i9,5e13.5)',iostat=fseekstatus) &
                           nsphere,neqns,alpha,beta,cbeam,dummy(1:2)
                     if(fseekstatus.ne.0) then
                        rewind(amnunit)
                        amnfileposition=0
                     else
                        backspace(amnunit)
                     endif
                     read(amnunit,'(2i9,5e13.5)') nsphere,neqns,alpha,beta,cbeam,dummy(1:2)
                     amnfileposition=amnfileposition+1
                     rimedium=dcmplx(dummy(1),dummy(2))
                     noff=0
                     do i=1,nsphere
                        read(amnunit,'(2i5,8e13.5)') nodr(i),numberfieldexp(i),xsp(i), &
                             rpos(:,i),dummy(1:4)
                        ri(1,i)=dcmplx(dummy(1),dummy(2))
                        ri(2,i)=dcmplx(dummy(3),dummy(4))
                        read(amnunit,'(9e13.5)') qext(i,:),qabs(i,:),qsca(i,:)
                        allocate(amnp1(0:nodr(i)+1,nodr(i),2),amnp2(0:nodr(i)+1,nodr(i),2))
                        nblk=2*nodr(i)*(nodr(i)+2)
                        do k=1,numberfieldexp(i)
                           do n=1,nodr(i)
                              do m=-n,n
                                 if(m.le.-1) then
                                    ma=n+1
                                    na=-m
                                 else
                                    ma=m
                                    na=n
                                 endif
                                 read(amnunit,'(4e17.9)') dummy(1:4)
                                 amnp1(ma,na,1)=dcmplx(dummy(1),dummy(2))
                                 amnp2(ma,na,1)=dcmplx(dummy(3),dummy(4))
                                 read(amnunit,'(4e17.9)') dummy(1:4)
                                 amnp1(ma,na,2)=dcmplx(dummy(1),dummy(2))
                                 amnp2(ma,na,2)=dcmplx(dummy(3),dummy(4))
                              enddo
                           enddo
                           amnp(noff+1:noff+nblk) &
                             =reshape(amnp1(0:nodr(i)+1,1:nodr(i),1:2),(/nblk/))
                           amnp(noff+nblk+1:noff+2*nblk) &
                             =reshape(amnp2(0:nodr(i)+1,1:nodr(i),1:2),(/nblk/))
                           noff=noff+2*nblk
                        enddo
                        deallocate(amnp1,amnp2)
                        amnfileposition=amnfileposition+2+nblk*numberfieldexp(i)
                     enddo
                     close(amnunit)
                  endif
                  call mstm_mpi(mpi_command='barrier')
               enddo
            endif
            oldamnfile=amnfile
!
!  calculate the efficiency factors
!
            call mstm_mpi(mpi_command='barrier')
            if(rank.eq.0) then
               cphi=cos(phi)
               sphi=sin(phi)
               qexttot=0.d0
               qabstot=0.d0
               qexttotpar=0.
               qexttotper=0.
               qabstotpar=0.
               qabstotper=0.
               do i=1,nsphere
                  qabsvol(i,:)=qabs(i,:)
                  do j=1,nsphere
                     if(hostsphere(j).eq.i) then
                        qabsvol(i,:)=qabsvol(i,:)-qabs(j,:)*xspn(j)*xspn(j)/xspn(i)/xspn(i)
                     endif
                  enddo
               enddo
               do i=1,nsphere
                  if(.not.tmonfile(i)) then
                     if(dimag(ri(1,i)).eq.0.d0.and.dimag(ri(2,i)).eq.0.d0) then
                        qabsvol(i,:)=0.d0
                     endif
                  endif
               enddo
               qabs=0.d0
               do i=1,nsphere
                  qabs(i,:)=qabsvol(i,:)
                  do j=1,nsphere
                     if(hostsphere(j).eq.i) then
                        qabs(i,:)=qabs(i,:)+qabsvol(j,:)*xspn(j)*xspn(j)/xspn(i)/xspn(i)
                     endif
                  enddo
               enddo
               do i=1,nsphere
                  if(hostsphere(i).eq.0) then
                     qexttotpar=qexttotpar+(qext(i,1)*cphi*cphi+2.d0*qext(i,3)*cphi*sphi &
                         +qext(i,2)*sphi*sphi)*xspn(i)*xspn(i)/xgeff/xgeff
                     qexttotper=qexttotper+(qext(i,1)*sphi*sphi-2.d0*qext(i,3)*cphi*sphi &
                         +qext(i,2)*cphi*cphi)*xspn(i)*xspn(i)/xgeff/xgeff
                     qabstotpar=qabstotpar+(qabs(i,1)*cphi*cphi+2.d0*qabs(i,3)*cphi*sphi &
                         +qabs(i,2)*sphi*sphi)*xspn(i)*xspn(i)/xgeff/xgeff
                     qabstotper=qabstotper+(qabs(i,1)*sphi*sphi-2.d0*qabs(i,3)*cphi*sphi &
                         +qabs(i,2)*cphi*cphi)*xspn(i)*xspn(i)/xgeff/xgeff
                     qexttot=qexttot+(qext(i,1)+qext(i,2))*0.5d0*xspn(i)*xspn(i)/xgeff/xgeff
                     qabstot=qabstot+(qabs(i,1)+qabs(i,2))*0.5d0*xspn(i)*xspn(i)/xgeff/xgeff
                  endif
               enddo
               qscatotpar=qexttotpar-qabstotpar
               qscatotper=qexttotper-qabstotper
               qscatot=qexttot-qabstot
            endif
            call rottranmtrxclear()
!
!  calculate the target-based expansion and rotate to the incident field frame
!
            if(numtheta.gt.0) then
               allocate(amnp0(0:nodrt+1,nodrt,2,2))
               call amncommonorigin(nsphere,nodr,ntran,nodrt,rposn, &
                       hostsphere,numberfieldexp,rimedium, &
                       amnp(1:neqns*2),amnp0,number_rhs=2)
               do k=1,2
                  call rotvec(alpha,beta,0.d0,nodrt,nodrt,amnp0(0:,1:,1:,k),1)
               enddo
!
!  calculate the asymmetry parameter and the scattering matrix
!
               if(rank.eq.0) then
                  allocate(gmn(0:2))
                  call s11expansion(amnp0,nodrt,0,1,gmn)
                  asymparm=dble(gmn(1)/gmn(0))/3.d0
                  deallocate(gmn)
               endif
               if(azimuthaverage.eq.1) then
                  if(allocated(s00)) deallocate(s00,s02,sp22,sm22)
                  nodrg=2*nodrt
                  allocate(s00(4,4,0:nodrg),s02(4,4,0:nodrg), &
                    sp22(4,4,0:nodrg),sm22(4,4,0:nodrg))
                  call fosmexpansion(nodrt,amnp0,s00,s02,sp22,sm22)
               endif
               smt=0.d0
               do i=1,numtheta
                  if(mod(i-1,numprocs).eq.rank) then
                     costheta=cos(scatteringanglearray(1,i)*pi/180.d0)
                     phi=scatteringanglearray(2,i)*pi/180.d0
!
! incidentortargetframe=1: rotate scattering direction to incident frame
!
                     if(incidentortargetframe.eq.1) then
                        sintheta=sin(scatteringanglearray(1,i)*pi/180.d0)
                        cphi=cos(phi)
                        sphi=sin(phi)
                        scatrotvec=(/sintheta*cphi,sintheta*sphi,costheta/)
                        call eulerrotation(scatrotvec,(/alpha,beta,0.d0/),1,scatrotvec)
                        costheta=scatrotvec(3)/sqrt(dot_product(scatrotvec,scatrotvec))
                        if(scatrotvec(1).eq.0.d0.and.scatrotvec(2).eq.0.d0) then
                           phi=0.d0
                        else
                           phi=datan2(scatrotvec(2),scatrotvec(1))
                        endif
                     endif
                     if(azimuthaverage.eq.0) then
                        call scatteringmatrix(amnp0,nodrt,costheta,phi,sa,smt(:,:,i))
                     else
                        call fosmcalc(nodrt,s00,s02,sp22,sm22,costheta,smt(:,:,i))
                     endif
                  endif
               enddo
               nsend=16*numtheta
               call mstm_mpi(mpi_command='allreduce',mpi_recv_buf_dp=smt(1:4,1:4,1:numtheta), &
                    mpi_number=nsend,mpi_operation=mstm_mpi_sum)
               call mstm_mpi(mpi_command='allreduce',mpi_recv_buf_dp=smtcf(1:4,1:4,1:numtheta), &
                    mpi_number=nsend,mpi_operation=mstm_mpi_sum)
               deallocate(amnp0)
            endif
         endif
!
!  output file operations
!
         if(rank.eq.0.and.fixedorrandom.le.1) then
     !   if(rank.eq.0.and.fixedorrandom.le.-100) then
            open(1,file=outfile,position='append')
            !write(1,'('' calculation results for wavelength '',/,f10.8)') wavlist(ilam)
            if(writespheredata.eq.1) then
               if(nonactive.eq.0) then
                  !write(1,'('' sphere host   ka     x-x(host) y-y(host) &
                  !      &z-z(host)'',&
                  !      & ''  Re(mL)    Im(mL)   '',&
                  !      & ''  Re(mR)    Im(mR)     '',&
                  !      & ''Qext         Qsca'',&
                  !      & ''        Qabs         Qabs(V)'')')
               else
                  !write(1,'(''sphere host   ka     x-x(host) y-y(host) &
                  !      &z-z(host)'',&
                  !      & ''  Re(m)     Im(m)      '',&
                  !      & ''Qext         Qsca'',&
                  !      & ''        Qabs         Qabs(V)'')')
                  !write(1,'(''xdp    ydp    zdp  wvln    sphere host   ka     x-x(host) y-y(host) &
                  !      &z-z(host)'',&
                  !      & ''  Re(m)     Im(m)      '',&
                  !      & ''Qext         Qsca'',&
                  !      & ''        Qabs         Qabs(V)'')')
               endif
               do i=1,nsphere
                  k=hostsphere(i)
                  if(fixedorrandom.eq.1) then
                     qeff(1)=qext(i,1)
                     qeff(2)=qsca(i,1)
                     qeff(3)=qabs(i,1)
                     qeff(4)=qabsvol(i,1)
                  else
                     qeff(1)=(qext(i,1)+qext(i,2))*0.5d0
                     qeff(2)=(qsca(i,1)+qsca(i,2))*0.5d0
                     qeff(3)=(qabs(i,1)+qabs(i,2))*0.5d0
                     qeff(4)=(qabsvol(i,1)+qabsvol(i,2))*0.5d0
                     !write(*,*) "qext(i,1), qext(i,2) = ", qext(i,1), qext(i,2)
                     !write(*,*) "qsca(i,1), qsca(i,2) = ", qsca(i,1), qsca(i,2)
                     !write(*,*) "qabs(i,1), qabs(i,2) = ", qabs(i,1), qabs(i,2)
                     !write(*,*) "qabsvol(i,1), qabsvol(i,2) = ", qabsvol(i,1), qabsvol(i,2)
                  endif
                  !if(k.eq.0) then
                  !   rposi=rpos(:,i)+gbfocus
                  !else
                  !   rposi=rpos(:,i)-rpos(:,k)
                  !endif
                  !I think these should change?
                  !The gbfocus needs to be scaled, but for printing you
                  !don't want them scaled...
                  if(k.eq.0) then
                     rposi=rpos(:,i)+gbfocus
                  else
                     rposi=rpos(:,i)-rpos(:,k)
                  endif
                  if(.not.tmonfile(i)) then
                     if(nonactive.eq.0) then
                        write(1,'(2i5,4f10.4,4f10.6,4e13.5)') i, k,xsp(i),rposi, ri(:,i), &
                           qeff(1),qeff(2),qeff(3),qeff(4)
                     else
                        !write(1,'(2i5,4f10.4,2f10.6,4e13.5)') i, k,xsp(i),rposi, ri(1,i), &
                        !   qeff(1),qeff(2),qeff(3),qeff(4)
                        write(1,'(3f10.6,A2,f5.3,2i5,4f10.4,2f10.6, &
                              4e13.5)') xdp(:,i), "  ", wavlist(ilam), &
                              i, k,xsp(i),rposi, ri(1,i), &
                              qeff(1),qeff(2),qeff(3),qeff(4)
                     endif
                  else
                     if(nonactive.eq.0) then
                        write(1,'(2i5,4f10.4,5x,a35,4e13.5)') i, k,xsp(i),rposi, tmfile(i), &
                           qeff(1),qeff(2),qeff(3),qeff(4)
                     else
                        write(1,'(2i5,4f10.4,3x,a17,4e13.5)') i, k,xsp(i),rposi, tmfile(i), &
                           qeff(1),qeff(2),qeff(3),qeff(4)
                     endif
                  endif
               enddo
            endif
            if(fixedorrandom.eq.1) then

               open(20,file=unpolscat,position='append')
               write(20,'(7e13.5)') wavlist(ilam), qexttot,qabstot,qscatot,asymparm
               !Original print statement
               !!write(1,'('' total ext, abs, scat efficiencies, &
               !!          & w.r.t. xv, and asym. parm'')')
               !!write(1,'(7e13.5)') wavlist(ilam), qexttot,qabstot,qscatot,asymparm
               close(20)
            else
               if(writescatdata.eq.1) then
                  open(20,file=unpolscat,position='append')
                  open(21,file=parascat, position='append')
                  open(22,file=perpscat, position='append')


                  write(20,'(7e13.5)')  wavlist(ilam), qexttot,qabstot,qscatot,asymparm
                  write(21,'(7e13.5)')  wavlist(ilam), qexttotpar,qabstotpar,qscatotpar
                  write(22,'(7e13.5)')  wavlist(ilam), qexttotper,qabstotper,qscatotper
                  close(20)
                  close(21)
                  close(22) 
                  !Original write statements
                  !write(1,'('' unpolarized total ext, abs, &
                  !               &  scat efficiencies, w.r.t. xv, &
                  !               &  and asym. parm'')')
                  !write(1,'(7e13.5)')  wavlist(ilam), qexttot,qabstot,qscatot,asymparm
                  !write(1,'('' parallel total ext, abs, &
                  !             & scat efficiencies'')')
                  !write(1,'(7e13.5)')  wavlist(ilam), qexttotpar,qabstotpar,qscatotpar
                  !write(1,'('' perpendicular total ext, abs, &
                  !             &  scat efficiencies'')')
                  !write(1,'(7e13.5)')  wavlist(ilam), qexttotper,qabstotper,qscatotper
               endif
            endif

            ! YJ 3-26-2019 disabling scattering matrix output
            !if(numtheta.gt.0) then
            !   if(normalizesm.eq.0) then
            !      write(1,'('' scattering matrix elements'')')
            !   else
            !      write(1,'('' scattering matrix elements &
            !                & (normalized w/ S11)'')')
            !   endif
            !   if(fixedorrandom.eq.0) then
            !      if(azimuthaverage.eq.0) then
            !         write(1,'(''   theta     phi'',$)')
            !      else
            !         write(1,'(''   theta'',$)')
            !      endif
            !      do i=1,4
            !         do j=1,4
            !            write(1,'(''      '',2i1,''     '',$)') i,j
            !         enddo
            !      enddo
            !   else
            !      write(1,'(''   theta'',$)')
            !      do i=1,4
            !         do j=i,4
            !            write(1,'(''      '',2i1,''     '',$)') i,j
            !         enddo
            !      enddo
            !   endif
            !   write(1,*)
            !   do i=1,numtheta
            !      thetad=scatteringanglearray(1,i)
            !      phi=scatteringanglearray(2,i)
            !      if(normalizesm.eq.1) then
            !         s11=smt(1,1,i)
            !         smt(:,:,i)=smt(:,:,i)/s11
            !         smt(1,1,i)=s11
            !      endif
            !      if(fixedorrandom.eq.0) then
            !         if(azimuthaverage.eq.0) then
            !            write(1,'(2f8.2,$)') thetad,phi
            !         else
            !            write(1,'(f8.2,$)') thetad
            !         endif
            !         do j=1,4
            !            do k=1,4
            !               write(1,'(e13.5,$)') smt(j,k,i)
            !            enddo
            !         enddo
            !      else
            !         write(1,'(f8.2,$)') thetad
            !         do j=1,4
            !            do k=j,4
            !               write(1,'(e13.5,$)') smt(j,k,i)
            !            enddo
            !         enddo
            !      endif
            !      write(1,*)
            !   enddo
            !   if(fixedorrandom.eq.1) then
            !      write(1,'('' scattering matrix expansion & 
            !             & coefficients'')')
            !      write(1,'(''    w  a11         a22         a33   '',&
            !       &''a23         a32         a44         a12      '',&
            !       &''a34         a13         a24         a14'')')
            !      do w=0,nodrg
            !         write(1,'(i5,11e12.4)') w,smc(1,1,w),smc(2,2,w),&
            !           smc(3,3,w),smc(2,3,w),smc(3,2,w),smc(4,4,w),&
            !           smc(1,2,w),smc(3,4,w),smc(1,3,w),smc(2,4,w),&
            !           smc(1,4,w)
            !      enddo
            !   endif
            !endif
        runnum = runnum+1
        endif
!
!  end of input file loop
!
!  CWH: Not sure what these mpi calls do but doing my best to not mess
!  with them.  This is the end of his input file loop, it'll be the end
!  of our wavelength loop
!  
!  I moved the near-field calculation out of the wavelength loop so you
!  don't run nearfieldgridcalc every time.  Probably never going to
!  bother modifying it to save the scattering coefficients from all
!  wavelengths and allow a user-selected near-field to be plotted.  Just
!  need to run a separate calculation
!
!
!  Write the E-field...intensity? to the acceptfile when dpcalc is 3
!  For use with Wendu's stuff
!
          ! Calculate RET rate and quantum yield
          if(dpcalc.eq.3) then
              open(8,file=acceptfile,position='append')
              open(80,position='append')
              open(100,position='append')
              call getrunparameters(plane_wave_epsilon=epspw)

              intqy = 0.80 ! donor intrinsic quantum yield
              lambda=wavlist(ilam)*1.d-6
              omega = 2*pi2*c/lambda
              totaccint = 0.0
              scataccint = 0.0
              incaccint = 0.0
              escatacc(1) = cmplx(0.0,0.0)
              escatacc(2) = cmplx(0.0,0.0)
              escatacc(3) = cmplx(0.0,0.0)
              eincacc(1) = cmplx(0.0,0.0)
              eincacc(2) = cmplx(0.0,0.0)
              eincacc(3) = cmplx(0.0,0.0)
              etotacc(1) = cmplx(0.0,0.0)
              etotacc(2) = cmplx(0.0,0.0)
              etotacc(3) = cmplx(0.0,0.0)
              
              xaccept=xaccept0
              xgpmax=10.0*xaccept

              nfnewcalc = 1
              call  nearfieldpointcalc(neqns,nsphere,nodr,alpha,beta, &
                     cbeam,xsp,rpos,ri,hostsphere, numberfieldexp,  &
                     amnp,phi,epspw,xaccept, nfnewcalc,escatacc,hscatacc, medk)
              call  nearfieldpointcalc(neqns,nsphere,nodr,alpha,beta, &
                     cbeam,xsp,rpos,ri,hostsphere, numberfieldexp,  &
                     amnp,phi,epspw,xaccept, nfnewcalc,escatacc,hscatacc, medk)
              !write(8,'(f8.5,6e16.8)') wavlist(ilam), eaccept(:)
              call nearfielddipolepart(xaccept, ri(1,0), &
                                       eincacc,hincacc,medk,dpmom)
              !einc(:) = einc(:)/sqrt(4*pi2*eps0)/sqrt(4*pi2*eps0)*2.d8
              !escat(:) = -1.d0*escat(:)/sqrt(4*pi2*eps0)/sqrt(4*pi2*eps0)*2.d8
              eincacc(:) = 2.d0*eincacc(:)
              escatacc(:) = 2.d0*escatacc(:)
              eincacc(:) = eincacc(:) / (4*pi2*eps0) * 1.0d18 * (-1.d0) !SI
              escatacc(:) = escatacc(:) / (4*pi2*eps0) * 1.0d18 * (-1.d0) !SI
              etotacc(:) = eincacc(:) + escatacc(:)
              totaccint = sqrt(abs(etotacc(1))**2 + abs(etotacc(2))**2 + abs(etotacc(3))**2)
              scataccint = sqrt(abs(escatacc(1))**2 + abs(escatacc(2))**2 + abs(escatacc(3))**2)
              incaccint = sqrt(abs(eincacc(1))**2 + abs(eincacc(2))**2 + abs(eincacc(3))**2)
              if(nfoutdata.eq.0) then
!                 Do k = 1,3
!                    write(*,*)  epfield(1,k)
!                 End Do
                  write(8,'(f10.5,4es17.9)') 1000*wavlist(ilam), totaccint,&
                                            incaccint, scataccint, totaccint/incaccint
              elseif(nfoutdata.eq.1) then
                  write(90,'(A6,23es17.9)') "dddd",wavlist(ilam)*1.d-6,totaccint/incaccint, &
                    totaccint, (dble(etotacc(k)), aimag(etotacc(k)),k=1,3), &
                    incaccint, (dble(eincacc(k)), aimag(eincacc(k)),k=1,3), &
                    scataccint, (dble(escatacc(k)), aimag(escatacc(k)),k=1,3)
!                 write(90,'(A6,9es17.9)')"dddd",wavlist(ilam)*1.d-6, &
!                   incint, (dble(einc(k)), aimag(einc(k)),k=1,3), totint/incint

!                   write(8,'(f10.5,7e12.4)') wavlist(ilam)*1.d-6, escat(:)

              elseif(nfoutdata.eq.2) then
                  write(8,'(f10.5,10e12.4)') 1000*wavlist(ilam), &
                                               etotacc(:), htotacc(:)
              elseif(nfoutdata.eq.3) then
!                 Do k = 1,3
!                    write(*,*)  epfield(1,k), accorient(k)
!                 End Do
                  incaccint2=cmplx(0.d0,0.d0)
                  totaccint2=cmplx(0.d0,0.d0)
                  scataccint2=cmplx(0.d0,0.d0)
                  do k = 1, 3
                    incaccint2=incaccint2+eincacc(k)*accorient(k)
                    totaccint2=totaccint2+etotacc(k)*accorient(k)
!                   escat(k)=cmplx(0.d0,0.d0)
!                   do i = 1, npart
!                     escat(k)=escat(k)+es(i,k)
!                   enddo
                    scataccint2=scataccint2+escatacc(k)*accorient(k)
                  enddo
                  write(*,*) totaccint2, incaccint2, scataccint2
                  totaccint=2*pi/hbar*abs(totaccint2)**2.d0
                  incaccint=2*pi/hbar*abs(incaccint2)**2.d0
                  scataccint=2*pi/hbar*abs(scataccint2)**2.d0
                  write(*,*) totaccint, incaccint, scataccint
!                 write(80,'(2A6,4es17.9)') "bbbb", "cccc",totint,
!           &        incint, scatint, totint/incint
                  totaccint = totaccint / (hbar * omega)
                  incaccint = incaccint / (hbar * omega)
                  scataccint = scataccint / (hbar * omega)
                  write(*,*) totaccint, incaccint, scataccint
                  write(80,'(A6,F9.3,4es17.9)') "dddd",wavlist(ilam)*1.d3, &
                     totaccint,incaccint, scataccint, totaccint/incaccint
              endif

              escatdon(1) = cmplx(0.0,0.0)
              escatdon(2) = cmplx(0.0,0.0)
              escatdon(3) = cmplx(0.0,0.0)
              eincdon(1) = cmplx(0.0,0.0)
              eincdon(2) = cmplx(0.0,0.0)
              eincdon(3) = cmplx(0.0,0.0)
              etotdon(1) = cmplx(0.0,0.0)
              etotdon(2) = cmplx(0.0,0.0)
              etotdon(3) = cmplx(0.0,0.0)
              
              nfnewcalc = 1
              xaccept=xdp0 ! the evaluation point is now at the donor
              xgpmax=10.0*xaccept
              call  nearfieldpointcalc(neqns,nsphere,nodr,alpha,beta, &
                     cbeam,xsp,rpos,ri,hostsphere, numberfieldexp,  &
                     amnp,phi,epspw,xaccept, nfnewcalc,escatdon,hscatdon, medk)
              call  nearfieldpointcalc(neqns,nsphere,nodr,alpha,beta, &
                     cbeam,xsp,rpos,ri,hostsphere, numberfieldexp,  &
                     amnp,phi,epspw,xaccept, nfnewcalc,escatdon,hscatdon, medk)
              !write(8,'(f8.5,6e16.8)') wavlist(ilam), eaccept(:)
              !call nearfielddipolepart(xaccept, ri(1,0), &
              !                         einc,hinc,medk,dpmom)
              !eincdon(:) = 2.d0*eincdon(:)
              escatdon(:) = 2.d0*escatdon(:)
              !eincdon(:) = eincdon(:) / (4*pi2*eps0) * 1.0d18 * (-1.d0) !SI
              escatdon(:) = escatdon(:) / (4*pi2*eps0) * 1.0d18 * (-1.d0) !SI
              !etotdon(:) = eincdon(:) + escatdon(:)

              ! Caution: only real part
              dpmommag = sqrt(RealPart(dpmom(1))**2+ &
                              ImagPart(dpmom(1))**2+ &
                              RealPart(dpmom(2))**2+ &
                              ImagPart(dpmom(2))**2+ &
                              RealPart(dpmom(3))**2+ &
                              ImagPart(dpmom(3))**2)
              !powerself = 4*pi2*dpmommag**2*omega*medk**3/3/ri(1,0)**2/cmicro
              powerself = dpmommag**2*omega**4/3/ri(1,0)**2/c**3/(4*pi2*eps0)
              !escat(:) = escat(:) * (-1.d0)
              !write(*,*) "escat = ", escat
              !powerscat = (2*pi2*omega)/cmicro*ImagPart(conjg(dpmom(1))*escat(1) + &
              !            conjg(dpmom(2))*escat(2) + &
              !            conjg(dpmom(3))*escat(3))
              powerscat = omega/2*ImagPart(conjg(dpmom(1))*escatdon(1) + &
                          conjg(dpmom(2))*escatdon(2) + &
                          conjg(dpmom(3))*escatdon(3))
              powerdecay = powerself + powerscat
         
              ratedecay = powerdecay / (hbar * omega)
              rateself = powerself / (hbar * omega)
              ratescat = powerscat / (hbar * omega)
              totqy = totaccint / &
                      (ratedecay + totaccint + &
                       (1-intqy) / intqy * rateself)
              vacqy = incaccint / &
                      (rateself + incaccint + &
                       (1-intqy) / intqy * rateself)
 
              totlifetime = 1 / totaccint ! caculated in second, print in picosecond
              vaclifetime = 1 / incaccint ! caculated in second, print in picosecond
              write(100,'(f10.2,6es12.4,2e12.2,2es12.4,2es12.2,es14.2)') wavlist(ilam)*1000, totaccint, &
                                           incaccint, scataccint, ratedecay, &
                                           rateself, ratescat, totqy, vacqy, &
                                           totlifetime*1.d12, vaclifetime*1.d12, &
                                           totaccint/incaccint
              close(8) 
          endif

         ! Calculate the total decay rate and radiative decay rate of the donor dipole
         ! YJ 4-2-2019
         ! There are many prefactor issues here. pmnpt and amnpt are (-1/2) times of
         ! what's in the Heaps and Schatz (2017) paper. Also, nmn's are 1/2 times, too.
         ! However, (-2)*amnpt*nmn gives the correct scattered electric field (compared to
         ! FDTD), instead of (-4)*amnpt*nmn. So I presume that plm and qlm in the paper
         ! is missing a 1/2 factor. However, using the 1/2 smaller amnpt and pmnpt results
         ! in 1/2 times smaller absoprtion power by sphere than the total decay rate, for
         ! small sphere cases. Without justification I multiply the absoprtion power by 2.
         ! However, this must be solved at some point.

         ! YJ 4-7-2019 The reference VSHs I was using were sqrt(2) times of what I was thinking.
         ! So here nmn, amnpt, pmnpt's are (-sqrt(2)) times of the standard definition (excluding
         ! sqrt(4*pi) and (-1)^m factors.
          if(dpcalc.eq.4) then
              !open(8,file=acceptfile,position='append')
              call getrunparameters(plane_wave_epsilon=epspw)

              nfnewcalc = 1
              totaccint = 0.0
              scataccint = 0.0
              incaccint = 0.0
              escatacc(1) = cmplx(0.0,0.0)
              escatacc(2) = cmplx(0.0,0.0)
              escatacc(3) = cmplx(0.0,0.0)
              eincacc(1) = cmplx(0.0,0.0)
              eincacc(2) = cmplx(0.0,0.0)
              eincacc(3) = cmplx(0.0,0.0)
              etotacc(1) = cmplx(0.0,0.0)
              etotacc(2) = cmplx(0.0,0.0)
              etotacc(3) = cmplx(0.0,0.0)
              
              nfnewcalc = 1
              xgpmax=10.0*xaccept

              call  nearfieldpointcalc(neqns,nsphere,nodr,alpha,beta, &
                     cbeam,xsp,rpos,ri,hostsphere, numberfieldexp,  &
                     amnp,phi,epspw,xaccept, nfnewcalc,escatacc,hscatacc, medk)
              call  nearfieldpointcalc(neqns,nsphere,nodr,alpha,beta, &
                     cbeam,xsp,rpos,ri,hostsphere, numberfieldexp,  &
                     amnp,phi,epspw,xaccept, nfnewcalc,escatacc,hscatacc, medk)
              call  nearfielddipolepart(xaccept, ri(1,0), &
                                       eincacc,hincacc,medk,dpmom)
              eincacc(:) = 2.d0*eincacc(:)
              escatacc(:) = 2.d0*escatacc(:)
              etotacc(:) = eincacc(:) + escatacc(:)

              lambda=wavlist(ilam)*1.d-6

              dpmommag = sqrt(RealPart(dpmom(1))**2+ &
                              RealPart(dpmom(2))**2+ &
                              RealPart(dpmom(3))**2)
              omega = 2*pi2*c/lambda
              !powerself = 4*pi2*dpmommag**2*omega*medk**3/3/ri(1,0)**2/cmicro
              powerself = dpmommag**2*omega**4/3/ri(1,0)**2/c**3/(4*pi2*eps0)
              !escat(:) = escat(:) * (-1.d0)
              escatacc(:) = escatacc(:) / (4*pi2*eps0) * 1.0d18 * (-1.d0)
              !write(*,*) "escat = ", escat
              !powerscat = (2*pi2*omega)/cmicro*ImagPart(conjg(dpmom(1))*escat(1) + &
              !            conjg(dpmom(2))*escat(2) + &
              !            conjg(dpmom(3))*escat(3))
              powerscat = omega/2*ImagPart(conjg(dpmom(1))*escatacc(1) + &
                          conjg(dpmom(2))*escatacc(2) + &
                          conjg(dpmom(3))*escatacc(3))
              powerdecay = powerself + powerscat

              if(allocated(noff_nf)) deallocate(noff_nf,nblk_nf)
              allocate(noff_nf(nsphere),nblk_nf(nsphere))
              noff_nf(1)=0
              do i=1,nsphere
                nblk_nf(i)=2*nodr(i)*(nodr(i)+2)
                if(i.lt.nsphere) noff_nf(i+1)=noff_nf(i)+nblk_nf(i)
              enddo
              do i=1,nsphere
                write(70,*) "examination"
                if(allocated(amnp1)) deallocate(amnp1,amnp2,amnp3)
                if(allocated(pmnp1)) deallocate(pmnp1,pmnp2,pmnp3)
                if(allocated(amnpt)) deallocate(amnpt)
                if(allocated(pmnpt)) deallocate(pmnpt)
                !  write(*,*) "nblk_nf(",i,") = ", nblk_nf(i)
                ! write(*,*) "nodr(",i,") = ", nodr(i)
                nblk = nodr(i) * (nodr(i) + 2)
                b11=(2*noff_nf(i))+1
                b12=(2*noff_nf(i))+nblk_nf(i)
                !allocate(nmn(2, nblk,3))
                allocate(amnpt(nblk, 2))
                allocate(amnp1(0:nodr(i)+1,nodr(i),2))
                allocate(amnp2(0:nodr(i)+1,nodr(i),2))
                allocate(amnp3(0:nodr(i)+1,nodr(i),2))
                allocate(pmnpt(nblk, 2))
                allocate(pmnp1(0:nodr(i)+1,nodr(i),2))
                allocate(pmnp2(0:nodr(i)+1,nodr(i),2))
                allocate(pmnp3(0:nodr(i)+1,nodr(i),2))
                amnp1=reshape(amnp(b11:b12), &
                             (/nodr(i)+2,nodr(i),2/))
                amnp2=reshape(amnp(b11+nblk_nf(i):b12+nblk_nf(i)), &
                             (/nodr(i)+2,nodr(i),2/))
                pmnp1=reshape(pmnp(b11:b12), &
                             (/nodr(i)+2,nodr(i),2/))
                pmnp2=reshape(pmnp(b11+nblk_nf(i):b12+nblk_nf(i)), &
                             (/nodr(i)+2,nodr(i),2/))
                !noff = noff + nblk
                amnp3(:,:,1) = amnp1(:,:,1) + amnp1(:,:,2)
                amnp3(:,:,2) = amnp1(:,:,1) - amnp1(:,:,2)
                pmnp3(:,:,1) = pmnp1(:,:,1) + pmnp1(:,:,2)
                pmnp3(:,:,2) = pmnp1(:,:,1) - pmnp1(:,:,2)
                !write(*,*) "amnp3 = ", amnp3
                call reorder_sph_vec(amnpt, amnp3, 2, nodr(i))
                call reorder_sph_vec(pmnpt, pmnp3, 2, nodr(i))

                !write(*,*) "amnpt before making it SI = ", amnpt
                !write(*,*) "pmnpt before making it SI = ", pmnpt
                
                !amnpt(:,:) = amnpt(:,:)*dsqrt(4*pi2)*dsqrt(2.d0) ! cgs
                !pmnpt(:,:) = pmnpt(:,:)*dsqrt(4*pi2)*dsqrt(2.d0)
                amnpt(:,:) = amnpt(:,:)/dsqrt(4*pi2)/eps0*1.d18*dsqrt(2.d0) ! SI
                pmnpt(:,:) = pmnpt(:,:)/dsqrt(4*pi2)/eps0*1.d18*dsqrt(2.d0)
                !write(*,*) "amnpt = ", amnpt
                !write(*,*) "pmnpt = ", pmnpt
                !write(50,*) "amnpt examination"
                do n=1,nodr(i)
                   nn1=n*(n+1)
                   write(50,*) "n = ", n, "coefficients"
                   do m=-n,n
                      mn=nn1+m
                      write(50, &
                      '(A,ES20.10,A,ES20.10,A,A,ES20.10,A,ES20.10,A)') &
                      "(",RealPart(amnpt(mn,1)),",", &
                      ImagPart(amnpt(mn,1)),")", "(", &
                      RealPart(amnpt(mn,2)),",", &
                      ImagPart(amnpt(mn,2)),")"
                   enddo
                enddo
                write(60,*) "pmnpt examination"
                do n=1,nodr(i)
                   nn1=n*(n+1)
                   write(60,*) "n = ", n, "coefficients"
                   do m=-n,n
                      mn=nn1+m
                      write(60, &
                      '(A,ES20.10,A,ES20.10,A,A,ES20.10,A,ES20.10,A)') &
                      "(",RealPart(pmnpt(mn,1)),",", &
                      ImagPart(pmnpt(mn,1)),")", "(", &
                      RealPart(pmnpt(mn,2)),",", &
                      ImagPart(pmnpt(mn,2)),")"
                   enddo
                enddo
                call calcpower(amnpt,pmnpt,nodr(i),abspower)
                !abspower=abspower * cmicro/(2.d0*omega* &
                !medk*RealPart(ri(1,0))) ! hostsphere(i)
                abspower=abspower/(2.d0*omega*mu0* &
                (2.d0*pi2/lambda)*RealPart(ri(1,0))) ! hostsphere(i)
                !abspower=abspower*(-1.d0)
                radpower=powerdecay-abspower
                !write(*,*) "abspower = ", abspower

                write(100,'(f10.2,5es17.9)') wavlist(ilam)*1000, &
                powerdecay,powerself,powerscat,abspower,radpower
              enddo
          endif

         call mstm_mpi(mpi_command='barrier')
      enddo !End loop over wavelengths










!  Goal is to save coefficients for nf wavelength


!  Definition in module
!         subroutine amncommonorigin(nsphere,nodr,ntran,nodrt,rpos,hostsphere, &
!                    numberfieldexp,rimed,amnp,amnp0,number_rhs,mpi_comm)

! Calling it in main
               !allocate(amnp0(0:nodrt+1,nodrt,2,2))
               !call amncommonorigin(nsphere,nodr,ntran,nodrt,rposn, &
               !        hostsphere,numberfieldexp,rimedium, &
               !        amnp(1:neqns*2),amnp0,number_rhs=2)
!Me calling it
!
!      if(allocated(amnp0)) deallocate(amnp0)
!      allocate(amnp0(0:nodrt+1,nodrt,2,2))
!      !allocate(amnp0(0:nodrt+1,nodrt,2,1))
!      knf = 2.*rimedium(1)*pi/wavlist(nfcalcindex)
!      rposn(:,:) = dble(knf)*rpos(:,:)
!      call amncommonorigin(nsphere,nodr,ntran,nodrt,rposn, &
!                         hostsphere,numberfieldexp,rimedium, &
!                          amnpnf(1:neqns*2),amnp0,number_rhs=2)
!
!      open(23,file=commorigin,status='replace', &
!             action='write')
!
!      write(23,'(i9,3e13.5)') nodrt, wavlist(nfcalcindex), rimedium(1)
!
!
!!  Mackowski says amnp0 is in the TE/TM basis.  Can I just write them
!!  normally?
!      do n=1,nodrt
!         do m=-n,n
!            if(m.le.-1) then
!               ma=n+1
!               na=-m
!            else
!               ma=m
!               na=n
!            endif
!           !write(23,'(2i5, 4e17.9)') n, m,amnp0(ma,na,1,1),amnp0(ma,na,2,1)
!           !write(23,'(2i5, 4e17.9)') n, m,sqrt(2.0)*amnp0(ma,na,:,:)
!           write(23,'(2i5, 4e17.9)') n, m,sqrt(2.0)*amnp0(ma,na,1,1), &
!                                          sqrt(2.0)*amnp0(ma,na,2,1)
!            !write(3,'(2i, 4e17.9)') n, m, amnp1(ma,na,2),amnp2(ma,na,2)
!         enddo
!      enddo
!      close(23)
!
!
!
!


!      noff=0
!      do i=1,nsphere
!         allocate(amnp1(0:nodr(i)+1,nodr(i),2),amnp2(0:nodr(i)+1,nodr(i),2))
!         nblk=2*nodr(i)*(nodr(i)+2)
!         do k=1,numberfieldexp(i)
!            amnp1=reshape(amnp(noff+1:noff+nblk), &
!                         (/nodr(i)+2,nodr(i),2/))
!            amnp2=reshape(amnp(noff+nblk+1:noff+2*nblk), &
!                         (/nodr(i)+2,nodr(i),2/))
!            do n=1,nodr(i)
!               do m=-n,n
!                  if(m.le.-1) then
!                     ma=n+1
!                     na=-m
!                  else
!                     ma=m
!                     na=n
!                  endif
!                  write(3,'(4e17.9)') amnp1(ma,na,1),amnp2(ma,na,1)
!                  write(3,'(4e17.9)') amnp1(ma,na,2),amnp2(ma,na,2)
!               enddo
!            enddo
!            noff=noff+2*nblk
!         enddo
!         deallocate(amnp1,amnp2)
!      enddo









!
!  near field calculation options
!
!  Moved outside loop

      if(fixedorrandom.eq.0.and.calcnf.eq.1) then

         open(1,file=outfile,position='append')
            write(1,*) "Calculating the near-field for &
                                & lambda = ", wavlist(nfcalcindex)
         close(1)

         call getrunparameters(near_field_plane_coord=nfplane, &
              near_field_plane_position=nfplanepos,near_field_plane_vertices=nfplanevert, &
              spacial_step_size=deltax,polarization_angle_deg=gammadeg, &
              near_field_output_file=nfoutfile,near_field_output_data=nfoutdata, &
              plane_wave_epsilon=epspw,append_nf_file=appendnffile)
         nfoutunit=2
         if(rank.eq.0) then
            if(appendnffile.eq.0) then
               open(2,file=nfoutfile,status='replace',action='write')
            else
               open(2,file=nfoutfile,position='append')
            endif
         endif
         gamma=gammadeg*pi/180.d0
         !medium wavenumber at the correct wavelength
         knf = 2.*rimedium(1)*pi/wavlist(nfcalcindex)
         call nearfieldgridcalc(neqns,nsphere,nodr,alpha,beta,cbeam,xsp,rpos,ri, &
                    hostsphere,numberfieldexp,amnpnf,nfplane,nfplanepos,nfplanevert, &
                    gbfocus,deltax,gamma,nfoutunit,epspw,nfoutdata,runprintunit, knf)
         if(rank.eq.0) then
            close(nfoutunit)
         endif
      endif
!
!  all done!
!

      if(allocated(rposn)) deallocate(rposn,xspn)
      if(allocated(refmwav)) deallocate(refmwav)
      if(allocated(wavlist)) deallocate(wavlist)
      if(allocated(amnpnf)) deallocate(amnpnf)
      if(rank.eq.0) close(1)
      call mstm_mpi(mpi_command='barrier')
      call mstm_mpi(mpi_command='finalize')
      end
