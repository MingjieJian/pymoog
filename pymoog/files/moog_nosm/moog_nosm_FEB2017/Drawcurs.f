
      subroutine drawcurs
c******************************************************************************
c     This subroutine draws an arrow and writes the user x- and y-positions
c     of the cursor upon a click
c******************************************************************************

      include 'Pstuff.com'
      integer ichr


c      call sm_graphics
c      if     (whichwin .eq. '1of1') then
c         call sm_window (1,1,1,1,1,1)
c      elseif (whichwin .eq. '2of2') then
c         call sm_defvar ('y_gutter','0.0')
c         call sm_window (1,2,1,1,1,1)
c      endif
c      call sm_curs (xplotpos,yplotpos,ichr)


c      call sm_relocate (xplotpos,yplotpos-0.10*(yhi-ylo))
c      call sm_draw (xplotpos,yplotpos)
c      call sm_draw (xplotpos-0.01*(xhi-xlo),yplotpos-0.03*(yhi-ylo))
c      call sm_relocate (xplotpos,yplotpos)
c      call sm_draw (xplotpos+0.01*(xhi-xlo),yplotpos-0.03*(yhi-ylo))
      call writenumber (xplotpos)
c      call sm_expand (0.6)
c      call sm_relocate (xplotpos,yplotpos-0.11*(yhi-ylo))
c      call sm_putlabel (5,array)
      call writenumber (yplotpos)
c      if     (whichwin(4:4) .eq. '1') then
c         call sm_relocate (xplotpos,yplotpos-0.15*(yhi-ylo))
c      elseif (whichwin(4:4) .eq. '2') then
c         call sm_relocate (xplotpos,yplotpos-0.18*(yhi-ylo))
c      elseif (whichwin(4:4) .eq. '3') then
c         call sm_relocate (xplotpos,yplotpos-0.21*(yhi-ylo))
c      endif
c      call sm_putlabel (5,array)


c      call sm_gflush
c      call sm_alpha


      
      return
      end



