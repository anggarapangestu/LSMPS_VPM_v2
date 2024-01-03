module memory_fmm_2d

integer,parameter:: &
    lmax        = 6,                        & ! max of FMM levels
    npcm        = 100,                      &
    nfwrd       = 8,                        & ! 8, ori
    nbmax       = nfwrd*(4**(lmax+1)-4)/3,  &
    nbl1        = 8*lmax,                   &
    nbl2        = 256,                      &
    nbl3        = 256,                      &
    nbl4        = 16,                       &
    nbmrl       = nfwrd*1024
    ! icutoff     = 2,                        &
    ! n_s         = 10,                       & !
    ! ndp         = 10,                       & !
    ! n_inter     = 1
 end module memory_fmm_2d         
 

!!!    nwmax       = 1600,                     & ! max of wall elements
!!!    npmax_init  = 5d4,                      & ! max of particles
!!!    ngsmax      = 20,                       &  
