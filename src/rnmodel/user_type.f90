module user_type
    use constants, only: KINDR
    implicit none
    type user_para
     ! interval (x-range) to sample x values from
     real(KINDR), dimension(2) :: interval
     ! interval (x-range) to scale estimated x values to
     ! these values will be used only if doing a reaction norm analysis
     ! but two real values must be input regardless
     real(KINDR), dimension(2) :: chalvl
     integer :: nFarm
     real(KINDR) :: farmRange
     integer :: allocation ! allocation scenario 1,random;2,clust
     integer :: analysisType ! 1=covariate analysis, 2 = single trait
     ! do a reaction norm for covariate analysis 
     ! or take farms as fixed effect for single trait analysis
     logical :: RN
    end type user_para
    ! global variable userParam only accessible to the user. i.e., 
    !  all the subroutines written by the user can inherit this object
    !  by "use user_type" and modify/read it.
    ! The main program does not and need not know about this type
    type(user_para) :: userParam
end module user_type