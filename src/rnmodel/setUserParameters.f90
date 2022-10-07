subroutine setUserParameters(xmin, xmax, nfarm, farmRange, allocation)
    ! Almost everything in the subroutine is upto user. There may not even be a need
    ! to have a parameter file or the userParameter derived type can be empty
    ! hence, this subroutine can also be empty. But for compatibility reasons
    ! it needs to stay and be called
    use constants
    use global_module
    use user_type
    implicit none

    real(KINDR), intent(in) :: xmin, xmax, farmrange
    integer, intent(in) :: nfarm, allocation

    33 format(a34,": ", l1)
    34 format(a34,": ", i0)
    35 format(a34,": ", g0.15)
    301 format(18x,f3.1)

    !!!!!!!!!! reading parameter file !!!!!!!!!!!!!!

     userParam%interval(1) = xmin
     userParam%interval(2) = xmax
     userParam%chalvl(1) = xmin
     userParam%chalvl(2) = xmax
     
    ! number of farms
     userParam%nFarm = nfarm
    
    ! farm range (input as proportion)
     userParam%farmRange = farmRange 

    ! allocation scenario
     userParam%allocation = allocation

    ! re-input analysis type (covariate, single) and reaction-norm
     userParam%analysisType = 1
    ! reaction norm analysis
     userParam%RN = .false.
end subroutine setUserParameters