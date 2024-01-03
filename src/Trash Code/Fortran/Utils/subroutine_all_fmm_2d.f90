!!!!============================================================
!!!!==== subroutine for biotsavart-fmm 2D ======================
!!!!============================================================

!!!!============================================================
!!!!============================================================
subroutine par_loc(n0,n1,npi,xi,yi, &
xb_mini,xb_maxi,yb_mini,yb_maxi,nip,ipp)

implicit real(8)(a-h,o-z), integer(i-n)
dimension xi(npi), yi(npi)
dimension ipp(npi)

nip = 0
do i = n0,n1
    if( (xi(i)>=xb_mini).and.(xi(i)<xb_maxi).and. &
        (yi(i)>=yb_mini).and.(yi(i)<yb_maxi) )then
        nip = nip + 1
        ipp(nip) = i
    end if
end do

if(nip==0)then
    ipp(1) = 0
end if

return
end

!!!!==============================================================
!!!!==============================================================
subroutine amount_inbox(n0,n1,npi,xi,yi,n2,n3,npj,xj,yj, &
xb_mini,xb_maxi,yb_mini,yb_maxi,npr_in,n_inter)

implicit real(8)(a-h,o-z), integer(i-n)

dimension xi(npi), yi(npi)
dimension xj(npj), yj(npj)

if (n_inter==1)then
    npr_in = 0
    do i = n0,n1
        if( (xi(i)>=xb_mini).and.(xi(i)<xb_maxi).and. &
            (yi(i)>=yb_mini).and.(yi(i)<yb_maxi) )then        
            npr_in = npr_in + 1  
        end if
    end do

else if (n_inter==2)then
    npr_in = 0
    do j = n2,n3
        if( (xj(j)>=xb_mini).and.(xj(j)<xb_maxi).and. &
            (yj(j)>=yb_mini).and.(yj(j)<yb_maxi) )then        
            npr_in = npr_in + 1 
        end if
    end do
    do i = n0,n1
        if( (xi(i)>=xb_mini).and.(xi(i)<xb_maxi).and. &
            (yi(i)>=yb_mini).and.(yi(i)<yb_maxi) )then           
            npr_in = npr_in + 1
        end if
    end do

else if (n_inter==3)then
    npr_in = 0
    do i = n0,n1
        if( (xi(i)>=xb_mini).and.(xi(i)<xb_maxi).and. &
            (yi(i)>=yb_mini).and.(yi(i)<yb_maxi) )then        
            npr_in = npr_in + 1 
        end if
    end do
    do j = n2,n3
        do i = n0,n1
            isame = 0
            if( (xj(j)==xi(i)).and.(yj(j)==yi(i)) )then
                isame = 1
            end if
        end do
        if( (xj(j)>=xb_mini).and.(xj(j)<xb_maxi).and. &
            (yj(j)>=yb_mini).and.(yj(j)<yb_maxi).and.(isame==0) )then 
            npr_in = npr_in + 1 
        end if
    end do

end if

return
end

!!!=========================================================
!!!=========================================================
subroutine hierarchy_mesh &
(n0,n1,npi,xi,yi,n2,n3,npj,xj,yj,n_s,n_inter,xmin,ymin,xmax,ymax, &
 nb,xb_min,yb_min,xb_max,yb_max,xb_cen,yb_cen,lev, &
 nchild,ichild,iparent)

!use parameters_2d, only: npmax
use memory_fmm_2d                          
implicit real(8)(a-h,o-z), integer(i-n)
real(8),parameter:: tol2 = 1.0d-15
                
dimension xi(npi), yi(npi), xj(npj), yj(npj)
dimension xb_min(nbmrl,lmax), xb_max(nbmrl,lmax)
dimension yb_min(nbmrl,lmax), yb_max(nbmrl,lmax)
dimension xb_cen(nbmrl,lmax), yb_cen(nbmrl,lmax)
dimension nb(lmax), nprin(nbmrl,lmax)
dimension nchild(nbmrl,lmax), ichild(nbmrl,lmax,4),iparent(nbmrl,lmax)

xmin =  1.0d0/tol2
xmax = -1.0d0/tol2
ymin =  1.0d0/tol2
ymax = -1.0d0/tol2

do i = n0,n1
    xmin = min( xmin,xi(i) )
    xmax = max( xmax,xi(i) )
    ymin = min( ymin,yi(i) )
    ymax = max( ymax,yi(i) )
end do

do j = n2,n3
    xmin = min( xmin,xj(j) )
    xmax = max( xmax,xj(j) )
    ymin = min( ymin,yj(j) )
    ymax = max( ymax,yj(j) )
end do

tlout = 1.0d-5
dr = min( xmax-xmin, ymax-ymin )/dble(2**lmax)
xmin = xmin - dr*tlout
xmax = xmax + dr*tlout
ymin = ymin - dr*tlout
ymax = ymax + dr*tlout

rdx = (xmax-xmin)
rdy = (ymax-ymin)

if (rdy>=rdx) then
    ratio = rdy/rdx
    if (ratio<2.75d0)then
        cent = (xmin+xmax)/2.0d0
        xmin = cent - rdy/2.0d0
        xmax = cent + rdy/2.0d0
        db = rdy/2.0d0
        nb1 = 2
        nb2 = 2
    else if (ratio>=2.75d0)then
        cent = (ymin+ymax)/2.0d0
        nbt = ceiling(rdy/rdx)
        ymin = cent - rdx*nbt/2.0d0
        ymax = cent + rdx*nbt/2.0d0
        db = rdx
        nb1 = 1
        nb2 = nbt
    end if
else if (rdx>rdy) then
    ratio = rdx/rdy
    if (ratio<2.75d0)then
        cent = (ymin+ymax)/2.0d0
        ymin = cent - rdx/2.0d0
        ymax = cent + rdx/2.0d0
        db = rdx/2.0d0
        nb1 = 2
        nb2 = 2
    else if (ratio>=2.75d0)then
        cent = (xmin+xmax)/2.0d0
        nbt = ceiling(rdx/rdy)
        xmin = cent - rdy*nbt/2.0d0
        xmax = cent + rdy*nbt/2.0d0
        db = rdy
        nb1 = nbt
        nb2 = 1
    end if
end if

ib = 0
do ib1 = 1,nb1
    do ib2 = 1,nb2
        xb_mint = xmin + (ib1-1)*db
        xb_maxt = xmin + ib1*db
        yb_mint = ymin + (ib2-1)*db
        yb_maxt = ymin + ib2*db
        call amount_inbox(n0,n1,npi,xi,yi,n2,n3,npj,xj,yj,xb_mint,xb_maxt, &
                          yb_mint,yb_maxt,npr_in,n_inter)
        if(npr_in>0) then
            ib = ib + 1
            nb(1) = ib
            nprin(ib,1) = npr_in
            iparent(ib,1) = 0
            xb_min (ib,1) = xb_mint
            xb_max (ib,1) = xb_maxt
            yb_min (ib,1) = yb_mint
            yb_max (ib,1) = yb_maxt
            xb_cen (ib,1) = ( xb_maxt+xb_mint )/2.0d0
            yb_cen (ib,1) = ( yb_maxt+yb_mint )/2.0d0
            nchild(ib,1) = 0
        end if
    end do
end do
nb = ib

i_stop = 0
do k = 1,lmax-1
    k1 = k + 1
    nb(k1) = 0
    n_stop = 0
    do ib = 1,nb(k)
        nchild(ib,k) = 0       
        if (nprin(ib,k)>n_s) then
            n_stop = n_stop + 1
            dbx = (xb_max(ib,k)-xb_min(ib,k))/2.0d0
            dby = (yb_max(ib,k)-yb_min(ib,k))/2.0d0
            do ib1 = 1,2
                do ib2 = 1,2
                    xb_mint = xb_min(ib,k) + (ib1-1)*dbx
                    xb_maxt = xb_min(ib,k) + ib1*dbx
                    yb_mint = yb_min(ib,k) + (ib2-1)*dby
                    yb_maxt = yb_min(ib,k) + ib2*dby
                    call amount_inbox(n0,n1,npi,xi,yi,n2,n3,npj,xj,yj,xb_mint,xb_maxt, &
                                      yb_mint,yb_maxt,npr_in,n_inter)
                    if(npr_in>0) then
                        nb(k1) = nb(k1) + 1
                        nchild(ib,k) = nchild(ib,k) + 1
                        ichild(ib,k,nchild(ib,k)) = nb(k1)
                        nprin(nb(k1),k1) = npr_in
                        xb_min (nb(k1),k1) = xb_mint
                        xb_max (nb(k1),k1) = xb_maxt
                        yb_min (nb(k1),k1) = yb_mint
                        yb_max (nb(k1),k1) = yb_maxt
                        xb_cen (nb(k1),k1) = ( xb_maxt+xb_mint )/2.0d0
                        yb_cen (nb(k1),k1) = ( yb_maxt+yb_mint )/2.0d0
                        iparent(nb(k1),k1) = ib
                    end if
                end do
            end do
        end if
    end do
    
    lev = k1
    if(n_stop==0)then
        lev = k
        exit
    end if
    
end do
nchild(1:nb(lev),lev) = 0

! open(55,file='grid_fmm.dat',status='replace')
! do k = 1,lev
!     do ib = 1,nb(k)
!         write(55,*) xb_min(ib,k),xb_max(ib,k),yb_min(ib,k),yb_max(ib,k),xb_cen(ib,k),yb_cen(ib,k)
!     end do
! end do
! close(55)
! open(57,file='grid_fmm2.dat',status='replace')
! do ip = n0,n1
!     write(57,*) xi(ip), yi(ip)
! end do
! if( (n_inter==2).or.(n_inter==3) )then
!     do ip = n2,n3
!         write(57,*) xj(ip), yj(ip)
!     end do
! end if
! close(57)

return
end
 
!!!=================================================================
!!! list 1    
subroutine list_one(ib,k,nb,lev,xb_min,yb_min,xb_max,yb_max, &
                    nchild,nls1,ils1,kls1)

use memory_fmm_2d
implicit real(8)(a-h,o-z), integer(i-n)
real(8),parameter:: tol2 = 1.0d-15

dimension xb_min(nbmrl,lmax), xb_max(nbmrl,lmax)
dimension yb_min(nbmrl,lmax), yb_max(nbmrl,lmax)
dimension nb(lmax), nchild(nbmrl,lmax)
dimension ils1(nbl1), kls1(nbl1)


if(nchild(ib,k)==0)then
    nls1 = 1
    ils1(nls1) = ib
    kls1(nls1) = k

    do k2 = 1,lev
        do ib2 = 1,nb(k2)
            istop = 0
            if((k2==k).and.(ib2==ib))then
                istop = 1
            end if
            if( (nchild(ib2,k2)==0).and.(istop==0) )then
            
                dx_xn =dabs( xb_max(ib2,k2)-xb_min(ib,k) )
                dx_nx =dabs( xb_min(ib2,k2)-xb_max(ib,k) )
                dy_xn =dabs( yb_max(ib2,k2)-yb_min(ib,k) )
                dy_nx =dabs( yb_min(ib2,k2)-yb_max(ib,k) )
                                        
                if( ( ((dy_xn<=tol2).or.(dy_nx<=tol2)).and. &
                      (xb_min(ib2,k2)<=xb_max(ib,k)+tol2).and. &
                      (xb_max(ib2,k2)>=xb_min(ib,k)-tol2) ).or. &
                    ( ((dx_xn<=tol2).or.(dx_nx<=tol2)).and. &
                      (yb_min(ib2,k2)<=yb_max(ib,k)+tol2).and. &
                      (yb_max(ib2,k2)>=yb_min(ib,k)-tol2) ) )then  
                    nls1 = nls1 + 1
                    ils1(nls1) = ib2
                    kls1(nls1) = k2
                end if
                                   
            end if
        end do
    end do
                  
else if(nchild(ib,k)>0)then
    nls1 = 0
end if

return
end

!!!=================================================================
!!! list 2  
subroutine list_two(ib,k,nb,xb_min,yb_min,xb_max,yb_max, &
                    iparent,nchild,ichild,nls2,ils2,kls2)

use memory_fmm_2d
implicit real(8)(a-h,o-z), integer(i-n)
real(8),parameter:: tol2 = 1.0d-15

dimension xb_min(nbmrl,lmax), xb_max(nbmrl,lmax)
dimension yb_min(nbmrl,lmax), yb_max(nbmrl,lmax), nb(lmax)
dimension iparent(nbmrl,lmax), nchild(nbmrl,lmax), ichild(nbmrl,lmax,4)
dimension ils2(nbl2), kls2(nbl2)


if(k==1)then
    nls2 = 0
    do ib1 =1,nb(1)
        if(ib1/=ib)then
           if( (xb_min(ib1,1)>xb_max(ib,1)+tol2).or. &
               (xb_max(ib1,1)<xb_min(ib,1)-tol2).or. &
               (yb_min(ib1,1)>yb_max(ib,1)+tol2).or. &
               (yb_max(ib1,1)<yb_min(ib,1)-tol2) )then
                nls2 = nls2 + 1
                ils2(nls2) = ib1
                kls2(nls2) = 1
            end if
        end if
    end do


else if (k>1)then
    nls2 = 0
    ipt = iparent(ib,k)
    kpt = k-1
    do ib2 = 1,nb(kpt)
        if( (nchild(ib2,kpt)>0).and.(ib2/=ipt) )then
            dx_xn =dabs( xb_max(ib2,kpt)-xb_min(ipt,kpt) )
            dx_nx =dabs( xb_min(ib2,kpt)-xb_max(ipt,kpt) )
            dy_xn =dabs( yb_max(ib2,kpt)-yb_min(ipt,kpt) )
            dy_nx =dabs( yb_min(ib2,kpt)-yb_max(ipt,kpt) )
                                        
            if( ( ((dy_xn<=tol2).or.(dy_nx<=tol2)).and. &
                  (xb_min(ib2,kpt)<=xb_max(ipt,kpt)+tol2).and. &
                  (xb_max(ib2,kpt)>=xb_min(ipt,kpt)-tol2) ).or. &
                ( ((dx_xn<=tol2).or.(dx_nx<=tol2)).and. &
                  (yb_min(ib2,kpt)<=yb_max(ipt,kpt)+tol2).and. &
                  (yb_max(ib2,kpt)>=yb_min(ipt,kpt)-tol2) ) )then
                do ic2 = 1,nchild(ib2,kpt)
                    ic = ichild(ib2,kpt,ic2)  
                    if (ic/=ib)then
                        if( (xb_min(ic,k)>xb_max(ib,k)+tol2).or. &
                            (xb_max(ic,k)<xb_min(ib,k)-tol2).or. &
                            (yb_min(ic,k)>yb_max(ib,k)+tol2).or. &
                            (yb_max(ic,k)<yb_min(ib,k)-tol2) )then
                            nls2 = nls2 + 1
                            ils2(nls2) = ic
                            kls2(nls2) = k
                        end if
                    end if
                end do
            end if
            
        end if
    end do
end if

return
end

!!!=================================================================
!!! list 3  
subroutine list_three(ib,k,nb,lev,xb_min,yb_min,xb_max,yb_max, &
                      nchild,ichild,nls3,ils3,kls3)

use memory_fmm_2d
implicit real(8)(a-h,o-z), integer(i-n)
real(8),parameter:: tol2 = 1.0d-15
                    
dimension xb_min(nbmrl,lmax), xb_max(nbmrl,lmax)
dimension yb_min(nbmrl,lmax), yb_max(nbmrl,lmax)
dimension nchild(nbmrl,lmax), ichild(nbmrl,lmax,4)
dimension ils3(nbl3), kls3(nbl3), nb(lmax)
                      

if (nchild(ib,k)==0) then
    nls3 = 0
    
    do k2 = k,lev
        do ib2 = 1,nb(k2)
        
            istop = 0
            if((k2==k).and.(ib2==ib))then
                istop = 1
            end if
            
            if( (nchild(ib2,k2)>0).and.(istop==0) )then
                dx_xn =dabs( xb_max(ib2,k2)-xb_min(ib,k) )
                dx_nx =dabs( xb_min(ib2,k2)-xb_max(ib,k) )
                dy_xn =dabs( yb_max(ib2,k2)-yb_min(ib,k) )
                dy_nx =dabs( yb_min(ib2,k2)-yb_max(ib,k) )
                                    
                if( ( ((dy_xn<=tol2).or.(dy_nx<=tol2)).and. &
                      (xb_min(ib2,k2)<=xb_max(ib,k)+tol2).and. &
                      (xb_max(ib2,k2)>=xb_min(ib,k)-tol2) ).or. &
                    ( ((dx_xn<=tol2).or.(dx_nx<=tol2)).and. &
                      (yb_min(ib2,k2)<=yb_max(ib,k)+tol2).and. &
                      (yb_max(ib2,k2)>=yb_min(ib,k)-tol2) ) )then  
                    do ic2 = 1,nchild(ib2,k2)
                        ic = ichild(ib2,k2,ic2)
                        kc = k2+1
                        if( (xb_min(ic,kc)>xb_max(ib,k)+tol2).or. &
                            (xb_max(ic,kc)<xb_min(ib,k)-tol2).or. &
                            (yb_min(ic,kc)>yb_max(ib,k)+tol2).or. &
                            (yb_max(ic,kc)<yb_min(ib,k)-tol2) )then
                            nls3 = nls3+1
                            ils3(nls3) = ic
                            kls3(nls3) = kc
                        end if
                    end do
                end if
                
            end if
        end do
    end do
    
else if(nchild(ib,k)>0)then        
    nls3 = 0
    
end if       

return
end

!!!=================================================================
!! list 4  
subroutine list_four(ib,k,nb,xb_min,yb_min,xb_max,yb_max, &
                     iparent,nchild,nls4,ils4,kls4)

use memory_fmm_2d
implicit real(8)(a-h,o-z), integer(i-n)
real(8),parameter:: tol2 = 1.0d-15

dimension xb_min(nbmrl,lmax), xb_max(nbmrl,lmax)
dimension yb_min(nbmrl,lmax), yb_max(nbmrl,lmax)
dimension nchild(nbmrl,lmax)
dimension iparent(nbmrl,lmax), nb(nbmrl)               
dimension ils4(nbl4), kls4(nbl4)


ipt = iparent(ib,k)
kpt = k-1

if(k==1)then
    nls4 = 0

else if(k>1)then
    nls4 = 0
    
    do k2 = 1,k-1
        do ib2 = 1,nb(k2)
            if(nchild(ib2,k2)==0)then

                dx_mxmn =dabs( xb_max(ib2,k2)-xb_min(ipt,kpt) )
                dx_mnmx =dabs( xb_min(ib2,k2)-xb_max(ipt,kpt) )
                dy_mxmn =dabs( yb_max(ib2,k2)-yb_min(ipt,kpt) )
                dy_mnmx =dabs( yb_min(ib2,k2)-yb_max(ipt,kpt) )
                
                if( ( ((dy_mxmn<=tol2).or.(dy_mnmx<=tol2)).and. &
                      (xb_min(ib2,k2)<=xb_max(ipt,kpt)+tol2).and. &
                      (xb_max(ib2,k2)>=xb_min(ipt,kpt)-tol2) ).or. &
                    ( ((dx_mxmn<=tol2).or.(dx_mnmx<=tol2) ).and. &
                      (yb_min(ib2,k2)<=yb_max(ipt,kpt)+tol2).and. &
                      (yb_max(ib2,k2)>=yb_min(ipt,kpt)-tol2) ) ) then
                
                    if( (xb_min(ib2,k2)>xb_max(ib,k)+tol2).or. &
                        (xb_max(ib2,k2)<xb_min(ib,k)-tol2).or. &
                        (yb_min(ib2,k2)>yb_max(ib,k)+tol2).or. &
                        (yb_max(ib2,k2)<yb_min(ib,k)-tol2) )then
                        nls4 = nls4 + 1
                        ils4(nls4) = ib2
                        kls4(nls4) = k2
                    end if
                end if
                     
            end if
        end do
    end do
    
end if
                            
return
end

!!!!==============================================================
!!!!==============================================================
subroutine bico(n,k,c)

implicit real(8)(a-h,o-z), integer(i-n)

if( (k<0).or.(k>n).or.(n<0) )then
    c = 0.0d0
    write(*,*) 'n = ', n, 'k = ', k
    write(*,*) 'error: binomial input are out of range'
else if( (k==0).or.(k==n).or.(n==0) )then
    c = 1.0d0   
else if( (k==1).or.(k==n-1) ) then
    c = dble(n)
else
    if (k>n/2) then 
        k = n-k
    end if
    n_k =  n-k
    c = 1.0d0
    do i = 1,k
        c = c * ( dble(n_k+i)/dble(i) )
    end do
end if

return
end

!!!!==============================================================
!!!!==============================================================
! 2D Regularization Function
subroutine regul_func_2d(rij,sij,q)

! use parameters_2d
! use memory_fmm_2d

implicit none

intent(in)  rij,sij
intent(out) q

real(8)     :: rij, sij, q,qt,qte
real(8)     :: rho,rho2,rij2,sij2,sij3,pi


            rij2 = rij**2
            sij2 = sij**2
            sij3 = dsqrt( sij2**3)
            rho2 = rij2/sij2
            rho  = dsqrt( rho2 )


pi          = 2.0d0*dacos(0.0d0)
! icutoff =  0.singular ; 1.super (high-oder) algebraic  ; 2. Gaussian  ; 3 super Gaussian 

! if     (icutoff == 0) then  ! Singular
!         q  = 1.0d0 /(2.0d0*pi)
! else if(icutoff == 1)then   ! High order algebraic
!         q  = ( (rho2*(rho2+2.0d0))/(rho2+1.0d0)**2 )/(2.0d0*pi)
! elseif (icutoff == 2) then  ! Gaussian 2nd order
        qt = -rho2/(2.0d0)
        q  = ( 1.0d0 - dexp(qt) )/(2.0d0*pi)
! elseif (icutoff == 3) then  ! Gaussian 3rd order
!         qt = -rho2/(2.0d0)
!         qte  = (1-qt)*dexp(qt)
!         q  = ( 1.0d0 - qte )/( 2.0d0*pi )
! endif

return
end subroutine regul_func_2d

!!!!==============================================================
!!!!==============================================================
subroutine direct_sum(nip,npi,ipp,xi,yi,si,ui,vi, &
                      njp,npj,jpp,xj,yj,sj,gj,icutoff)

implicit real(8)(a-h,o-z), integer(i-n)
                   
dimension xi(npi),yi(npi),si(npi),ui(npi),vi(npi)
dimension xj(npj),yj(npj),sj(npj),gj(npj)
dimension ipp(npi), jpp(npj)

pi = 2.0d0*dacos(0.0d0)

!$OMP PARALLEL DO PRIVATE(i2,j2,i,j,dxij,dyij,  &
!$OMP                     rij2,sij2,q,qt)
do i2 = 1,nip
    i = ipp(i2)
    do j2 = 1,njp
        j = jpp(j2)
    
        dxij = xi(i) - xj(j)
        dyij = yi(i) - yj(j)
        rij2 = dxij**2 + dyij**2
        sij2 = 0.5d0*( si(i)**2 + sj(j)**2 )

        rij     = dsqrt( rij2 )
        sij  = dsqrt(sij2)
        
        if(rij2>0.0d0)then
            call regul_func_2d(rij,sij,q)
            ui(i) = ui(i) - q*dyij*gj(j)/rij2
            vi(i) = vi(i) + q*dxij*gj(j)/rij2
        end if

    end do
end do
!$OMP END PARALLEL DO

return
end

!!!! =============end of subroutines ============================
!!!! =============end of subroutines ============================
!!!! =============end of subroutines ============================
