subroutine biotsavart_fmm_2d &
(n0,n1,npi,xi,yi,si,ui,vi, &
 n2,n3,npj,xj,yj,sj,gj,icutoff,n_s,n_inter,ndp)

use memory_fmm_2d
implicit real(8)(a-h,o-z), integer(i-n)

dimension xi(npi),yi(npi),si(npi),ui(npi),vi(npi)
dimension xj(npj),yj(npj),sj(npj),gj(npj)

complex(8) cmw, cmw2, ci, z, zi, zj, zp, zc, zo, zg, zm

allocatable :: xb_min(:,:), xb_max(:,:)
allocatable :: yb_min(:,:), yb_max(:,:)
allocatable :: xb_cen(:,:), yb_cen(:,:)
allocatable :: nb(:), iparent(:,:)
allocatable :: nchild(:,:), ichild(:,:,:)
allocatable :: ipp(:), jpp(:)
allocatable :: ils1(:), kls1(:)
allocatable :: ils2(:), kls2(:)
allocatable :: ils3(:), kls3(:)
allocatable :: ils4(:), kls4(:)
allocatable :: ibt(:,:)
complex(8), allocatable ::  ak(:,:), bk(:,:)

allocate( xb_min(nbmrl,lmax), xb_max(nbmrl,lmax) )
allocate( yb_min(nbmrl,lmax), yb_max(nbmrl,lmax) )
allocate( xb_cen(nbmrl,lmax), yb_cen(nbmrl,lmax) )
allocate( nb(lmax), iparent(nbmrl,lmax) )
allocate( nchild(nbmrl,lmax), ichild(nbmrl,lmax,4) )
allocate( ipp(npi), jpp(npj)  )
allocate( ils1(nbl1), kls1(nbl1) )
allocate( ils2(nbl2), kls2(nbl2) )
allocate( ils3(nbl3), kls3(nbl3) )
allocate( ils4(nbl4), kls4(nbl4) )
allocate( ibt(nbmrl,lmax) )
allocate( ak(nbmax,npcm), bk(nbmax,npcm) )

xb_min(:,:) = 0.0d0; xb_max(:,:) = 0.0d0
yb_min(:,:) = 0.0d0; yb_max(:,:) = 0.0d0
xb_cen(:,:) = 0.0d0; yb_cen(:,:) = 0.0d0
nb(:)       = 0;  iparent(:,:) = 0 
nchild(:,:) = 0; ichild(:,:,:) = 0
ipp(:)   = 0; jpp(:)   = 0
ils1(:)  = 0; kls1(:) = 0
ils2(:)  = 0; kls2(:) = 0
ils3(:)  = 0; kls3(:) = 0
ils4(:)  = 0; kls4(:) = 0
ibt(:,:) = 0;
ak(:,:)= dcmplx(0.0d0,0.0d0) 
bk(:,:)= dcmplx(0.0d0,0.0d0) 


pi = 2.0d0*dacos(0.0d0)
ci = dcmplx(0.0d0,1.0d0)

do i = n0,n1
    ui(i) = 0.0d0
    vi(i) = 0.0d0
end do


!!!!=== step 1 === tree(?)

call hierarchy_mesh &
(n0,n1,npi,xi,yi,n2,n3,npj,xj,yj,n_s,n_inter,xmin,ymin,xmax,ymax, &
 nb,xb_min,yb_min,xb_max,yb_max,xb_cen,yb_cen,lev, &
 nchild,ichild,iparent)

!!!!=============================================== 
ibtt = 0
do k = 1,lev
    do ib = 1,nb(k)
        ibtt = ibtt + 1
        ibt(ib,k) = ibtt
        do kp = 1,ndp+1
            if(kp<=ndp)then
                ak(ibt(ib,k),kp) = dcmplx(0.0d0,0.0d0)
            end if
            bk(ibt(ib,k),kp) = dcmplx(0.0d0,0.0d0)  
        end do
    end do
end do

!!!!=== step 2.1 ======================
do k = 1,lev
    do ib = 1,nb(k)
        if(nchild(ib,k)==0)then
            call par_loc(n2,n3,npj,xj,yj,xb_min(ib,k),xb_max(ib,k), &
                         yb_min(ib,k),yb_max(ib,k),njp,jpp)
            
            zm = dcmplx( xb_cen(ib,k) , yb_cen(ib,k) )            
            if(njp>0)then
                do kp = 1,ndp
                    do jp2 = 1,njp
                        jp = jpp(jp2)
                        zj = dcmplx( xj(jp) , yj(jp) )
                        ak(ibt(ib,k),kp) = ak(ibt(ib,k),kp) + gj(jp)*((zj-zm)**(kp-1))
                    end do
                end do
            end if
        end if
    end do
end do

!!!!=== step 2.2 ===  
do k = lev-1,1,-1
    do ib = 1,nb(k)
        if( nchild(ib,k)>0 )then
        
            zp = dcmplx( xb_cen(ib,k) , yb_cen(ib,k) ) 
            do ic2 = 1,nchild(ib,k)
                ic = ichild(ib,k,ic2)      
                zc = dcmplx( xb_cen(ic,k+1) , yb_cen(ic,k+1) )
                do kl = 1,ndp
                    do kp = 1,kl
                        call bico(kl-1,kp-1,cbi)
                        ak(ibt(ib,k),kl) = ak(ibt(ib,k),kl) + cbi* &
                                           ((zc-zp)**(kl-kp))*ak(ibt(ic,k+1),kp)                         
                    end do
                end do
            end do
            
        end if
    end do
end do

!!!!=== step 4 & 6 ===         
do k = 1,lev
    do ib = 1,nb(k)
        zg = dcmplx( xb_cen(ib,k),yb_cen(ib,k) ) 
        
        !!!! ==== step 4 ===
        call list_two(ib,k,nb,xb_min,yb_min,xb_max,yb_max, &
                      iparent,nchild,ichild,nls2,ils2,kls2)
        if(nls2>0)then 
            do ib1 = 1,nls2
                ib2 = ils2(ib1)
                k2 = kls2(ib1)
                zo  = dcmplx( xb_cen(ib2,k2),yb_cen(ib2,k2) )
                do kl =1,ndp+1
                    ll = kl - 1
                    do kp = 1,ndp
                        call bico(ll+kp-1,kp-1,cbi)
                        bk(ibt(ib,k),kl) = bk(ibt(ib,k),kl) + ((-1.0d0)**kp)*cbi* &
                                           ak(ibt(ib2,k2),kp)/( ((zo-zg)**(kp+ll)) )
                    end do
                end do      
            end do
        end if
        !!!! ==== step 6 ===
        call list_four(ib,k,nb,xb_min,yb_min,xb_max,yb_max, &
                       iparent,nchild,nls4,ils4,kls4)
        if(nls4>0)then
            do ib1 = 1,nls4
                ib2 = ils4(ib1)
                k2 = kls4(ib1)
                call par_loc(n2,n3,npj,xj,yj,xb_min(ib2,k2),xb_max(ib2,k2), &
                             yb_min(ib2,k2),yb_max(ib2,k2),njp,jpp)
                do jp2 = 1,njp
                    jp = jpp(jp2)
                    zo  = dcmplx( xj(jp),yj(jp) )
                    do kl =1,ndp+1
                        ll = kl - 1
                        bk(ibt(ib,k),kl) = bk(ibt(ib,k),kl) - gj(jp)/((zo-zg)**(1+ll))
                    end do
                end do      
            end do
        end if
        
    end do
end do            
       
!!!!=== step 7 ===
do k = 1,lev-1
    do ib = 1,nb(k)
        if(nchild(ib,k)>0)then
            zp  = dcmplx( xb_cen(ib,k),yb_cen(ib,k) ) 
            
            do ic = 1,nchild(ib,k)
                ib2 = ichild(ib,k,ic)
                zc  = dcmplx( xb_cen(ib2,k+1),yb_cen(ib2,k+1) )     
                do kl = 1,ndp+1
                    ll = kl - 1 
                    do kp2 = kl,ndp+1   
                        kp = kp2 - 1 
                        call bico(kp,ll,cbi)
                        bk(ibt(ib2,k+1),kl) = bk(ibt(ib2,k+1),kl) + cbi* &
                                         ((zc-zp)**(kp-ll))*bk(ibt(ib,k),kp2)

                    end do
                end do
            end do
            
        end if
    end do
end do

!!!!=== step 3 & 5 & 8 ===
do k = 1,lev
    do ib = 1,nb(k)
        if(nchild(ib,k)==0)then
            call par_loc(n0,n1,npi,xi,yi,xb_min(ib,k),xb_max(ib,k), &
                         yb_min(ib,k),yb_max(ib,k),nip,ipp)
            zg = dcmplx( xb_cen(ib,k),yb_cen(ib,k) )
            
            !!!! ==== step 3 ===
            call list_one(ib,k,nb,lev,xb_min,yb_min,xb_max,yb_max, &
                          nchild,nls1,ils1,kls1)
            if( (nip>0).and.(nls1>0) )then
                do il = 1,nls1
                    ib2 = ils1(il)
                    k2 = kls1(il)
                    call par_loc(n2,n3,npj,xj,yj,xb_min(ib2,k2),xb_max(ib2,k2), &
                                 yb_min(ib2,k2),yb_max(ib2,k2),njp,jpp)
                    if(njp>0)then
                        call direct_sum(nip,npi,ipp,xi,yi,si,ui,vi,njp,npj,jpp,xj,yj,sj,gj,icutoff)
                    end if  
                end do
            end if
            !!!! ==== step 5 ===
            call list_three(ib,k,nb,lev,xb_min,yb_min,xb_max,yb_max, &
                            nchild,ichild,nls3,ils3,kls3)
            if((nip>0).and.(nls3>0))then
                do ip2 = 1,nip
                    ip = ipp(ip2)
                    cmw = dcmplx(0.0d0,0.0d0)
                    z = dcmplx( xi(ip),yi(ip) )
                    
                    do il = 1,nls3
                        ic = ils3(il)
                        ikc = kls3(il)
                        zm  = dcmplx( xb_cen(ic,ikc),yb_cen(ic,ikc) )
                        do kp = 1,ndp  
                            cmw = cmw + ak(ibt(ic,ikc),kp)/((z-zm)**kp)
                        end do
                    end do
                    cmw2 = -ci*cmw/(2.0d0*pi)
                    ui(ip) = ui(ip) + dble(cmw2)
                    vi(ip) = vi(ip) - dimag(cmw2)
                end do
            end if
!            !!!! ==== step 8 ===
            if(nip>0)then            
                do ip2 = 1,nip
                    ip = ipp(ip2)
                    zi  = dcmplx( xi(ip),yi(ip) ) 
                    cmw = dcmplx(0.0d0,0.0d0)
                    do kl = 1,ndp+1
                        ll = kl - 1
                        cmw = cmw + bk(ibt(ib,k),kl)*((zi-zg)**ll)
                    end do
                    cmw2 = -ci*cmw/(2.0d0*pi)
                    ui(ip) = ui(ip) + dble(cmw2)
                    vi(ip) = vi(ip) - dimag(cmw2)
                end do
            end if

        end if
    end do
end do

deallocate( xb_min, xb_max )
deallocate( yb_min, yb_max )
deallocate( xb_cen, yb_cen )
deallocate( nb, iparent )
deallocate( nchild, ichild )
deallocate( jpp, ipp )
deallocate( ils1, kls1 )
deallocate( ils2, kls2 )
deallocate( ils3, kls3 )
deallocate( ils4, kls4 )
deallocate( ibt, ak, bk )

return
end
!!!!==============================================================
