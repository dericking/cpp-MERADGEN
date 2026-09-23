      subroutine ptrc_unit(u)
      implicit none
      integer u
      integer uu, ios
      character*1024 path
      logical opened
      save uu, opened
      data uu/0/, opened/.false./
      if (.not. opened) then
         opened=.true.
         path=' '
         call get_environment_variable('MERADGEN_PARITY_TRACE', path)
         if (len_trim(path) .eq. 0) then
            uu=0
         else
            uu=77
            open(uu, file=trim(path), status='replace', iostat=ios)
            if (ios .ne. 0) uu=0
         endif
      endif
      u=uu
      end

      subroutine ptrc_event()
      implicit none
      integer n, u
      save n
      data n/0/
      n=n+1
      call ptrc_unit(u)
      if (u .eq. 0) return
      write(u,'(A,1X,I0)') 'EVENT', n
      end

      subroutine ptrc_d(tag, x)
      implicit none
      character*(*) tag
      real*8 x
      integer*8 ih
      integer u
      call ptrc_unit(u)
      if (u .eq. 0) return
      ih=transfer(x, ih)
      write(u,'(A,1X,A,1X,A,Z16.16)') 'D', trim(tag), '0x', ih
      end

      subroutine ptrc_i(tag, iv)
      implicit none
      character*(*) tag
      integer iv, u
      call ptrc_unit(u)
      if (u .eq. 0) return
      write(u,'(A,1X,A,1X,I0)') 'I', trim(tag), iv
      end

      subroutine ptrc_v(iv, vvn, sinv, ds, da)
      implicit none
      integer iv, u
      real*8 vvn, sinv, ds, da
      integer*8 h1, h2, h3, h4
      call ptrc_unit(u)
      if (u .eq. 0) return
      h1=transfer(vvn, h1)
      h2=transfer(sinv, h2)
      h3=transfer(ds, h3)
      h4=transfer(da, h4)
      write(u,'(A,1X,I0,4(1X,A,Z16.16))') 'V', iv, '0x', h1, '0x', h2,
     &  '0x', h3, '0x', h4
      end

      subroutine ptrc_t(idx, tt1n, sinv, ds, da)
      implicit none
      integer idx, u
      real*8 tt1n, sinv, ds, da
      integer*8 h1, h2, h3, h4
      call ptrc_unit(u)
      if (u .eq. 0) return
      h1=transfer(tt1n, h1)
      h2=transfer(sinv, h2)
      h3=transfer(ds, h3)
      h4=transfer(da, h4)
      write(u,'(A,1X,I0,4(1X,A,Z16.16))') 'T', idx, '0x', h1, '0x', h2,
     &  '0x', h3, '0x', h4
      end

      subroutine ptrc_z(iz, zzn, sinv, ds, da)
      implicit none
      integer iz, u
      real*8 zzn, sinv, ds, da
      integer*8 h1, h2, h3, h4
      call ptrc_unit(u)
      if (u .eq. 0) return
      h1=transfer(zzn, h1)
      h2=transfer(sinv, h2)
      h3=transfer(ds, h3)
      h4=transfer(da, h4)
      write(u,'(A,1X,I0,4(1X,A,Z16.16))') 'Z', iz, '0x', h1, '0x', h2,
     &  '0x', h3, '0x', h4
      end
