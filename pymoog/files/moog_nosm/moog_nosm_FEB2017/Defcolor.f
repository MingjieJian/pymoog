
      subroutine defcolor (icolor)
c******************************************************************************
c     This routine decides on whether to call for a black & white hardcopy
c     plot or a colorful screen plot
c******************************************************************************

      include 'Pstuff.com'
      character*7 colors(8)


c*****assign colors to character arrays
      colors(1) = 'white  '
      colors(2) = 'red    '
      colors(3) = 'cyan   '
      colors(4) = 'yellow '
      colors(5) = 'green  '
      colors(6) = 'magenta'
      colors(7) = 'blue   '
      colors(8) = 'black  '


c      if (choice.eq.'h' .or. choice.eq.'f' .or.
c     .    choice.eq.'g') then
c         call sm_ctype (colors(8))
c      else
c         call sm_ctype (colors(icolor))
c      endif

      
      return
      end




