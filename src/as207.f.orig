      subroutine gllm(istop, ni, nid, nj, nk, nkp, ji, y, c, conv, 
     *  w, v, e, f, cspr, cslr, ifault)
C
C       Algorithm AS 207 Appl. Statist. (1984) vol.33, no.3
C       modified by David Duffy to stop after specified number (istop) of
C       EM iterations
C
C       Fitting a generalized log-linear model
C       to fully or partially classified frequencies.
C
      integer ji(ni), istop, it
      double precision y(nj), c(nid, nkp), w(4, ni), v(2, nkp), e(ni),
     *  f(nj), eps, zero, one, cmin, conv, cslr, cspr, csum, ctmax, ff,
     *  yy
      logical lconv
      data eps, zero, one/0.00001d0, 0.0d0, 1.0d0/
C
      if (istop.le.0) istop=1
      ifault = 1
C
C       check the ji array
C
      do 100 i = 1, ni
      if (ji(i) .lt. 1 .or. ji(i) .gt. nj) return
  100 continue
      ifault = 0

C
C       initialize
C
      do 110 i = 1, ni
  110 e(i) = one
      do 120 j = 1, nj
  120 w(3, j) = zero
C
C       standardize the c matrix
C
      cmin = zero
      do 130 i = 1, ni
      do 130 k = 1, nk
      cmin = min(cmin, c(i, k))
  130 continue
      if (cmin .eq. zero) goto 150
      do 140 i = 1, ni
      do 140 k = 1, nk
        c(i, k) = c(i, k) - cmin
  140 continue
  150 ctmax = zero
      do 170 i = 1, ni
        csum = zero
        do 160 k = 1, nk
  160   csum = csum + c(i, k)
        ctmax = max(ctmax, csum)
        w(4, i) = csum
  170 continue
      if (ctmax .gt. eps) goto 180
      ifault = 2
      return
  180 if (abs(ctmax - one) .le. eps) goto 200
      do 190 i = 1, ni
        do 190 k = 1, nk
          c(i, k) = c(i, k) / ctmax
  190 continue
      goto 150
  200 do 210 i = 1, ni
        if (abs(w(4, i) - one) .gt. eps) goto 220
  210 continue
      nkk = nk
      goto 300
  220 nkk = nkp
      do 230 i = 1, ni
  230 c(i, nkk) = one - w(4, i)
C
C       enter the EM algorithm
C
       it=0

  300  continue
       it=it+1
       do 310 j = 1,nj
  310  f(j) = zero
       do 320 i = 1,ni
         j = ji(i)
         f(j) = f(j) + e(i)
  320  continue
C
C       check for convergence
C
      lconv = .true.
      do 330 j = 1, nj
        if (abs(f(j) - w(3, j)) .gt. conv) lconv = .false.
        w(3, j) = f(j)
  330 continue
      if (lconv .or. it.gt.istop) goto 500
C
      do 340 i = 1, ni
        j = ji(i)
        w(1, i) = y(j)
        if (f(j) .gt. eps) w(1, i) = e(i) * y(j) / f(j)
  340 continue
      do 360 k = 1, nkk
        v(1, k) = zero
        do 350 i = 1, ni
  350   v(1, k) = v(1, k) + c(i, k) * w(1, i)
  360 continue
C
C       enter the IPF algorithm
C
  400 do 410 i = 1, ni
  410 w(2, i) = e(i)
      do 440 k = 1, nkk
        v(2, k) = zero
        do 420 i = 1, ni
  420 v(2, k) = v(2, k) + c(i,k) * e(i)
      do 430 i = 1, ni
  430 if (c(i, k) .gt. eps .and. v(2, k) .gt. eps)
     *  e(i) = e(i) * (v(1, k) / v(2, k)) ** c(i, k)
  440 continue
      do 450 i = 1, ni
        if (abs(e(i) - w(2, i)) .gt. conv) goto 400
  450 continue
      goto 300
C
C       calculate the goodness-of-fit statistics
C
  500 cspr = zero
      cslr = zero
      do 530 j = 1, nj
        yy = y(j)
        ff = f(j)
        if (ff .le. eps) goto 530
        cspr = cspr + (yy - ff) ** 2 / ff
        if (yy .le. eps) goto 530
        cslr = cslr + yy * log(yy / ff)
  530 continue
      cslr = 2.0d0 * cslr
      return
      end
