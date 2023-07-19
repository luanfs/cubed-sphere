module deallocation
!========================================================================
!
! Module for memory deallocation routines
!
! Routines based on iModel (https://github.com/pedrospeixoto/iModel)
!
! Luan Santos 2022
!========================================================================

use datastruct, only: &
    cubedsphere

implicit none

contains 


subroutine meshdeallocation(mesh)
    !---------------------------------------------------
    ! MESHDEALLOCATION
    ! deallocate all the needed mesh atributtes
    ! except the tangent vectors
    !--------------------------------------------------
    type(cubedsphere):: mesh

    deallocate(mesh%po) ! Vertices
    deallocate(mesh%pc) ! Centers
    deallocate(mesh%pu) ! Midpoints at u
    deallocate(mesh%pv) ! Midpoints at v
    deallocate(mesh%ll2contra_pu) ! latlon 2 contravariant conversion at u
    deallocate(mesh%ll2contra_pv) ! latlon 2 contravariant conversion at v
    deallocate(mesh%contra2ll_pu) ! contravariant 2 latlon conversion at u
    deallocate(mesh%contra2ll_pv) ! contravariant 2 latlon conversion at v 
    deallocate(mesh%mt_pc)! Metric tensor at centers
    deallocate(mesh%mt_pu) ! Metric tensor at u
    deallocate(mesh%mt_pv) ! Metric tensor at v
    deallocate(mesh%ix_ll, mesh%jy_ll, mesh%panels_ll) ! Latlon grid indexes on cubedsphere

end subroutine meshdeallocation 

subroutine tgvectors_deallocation(mesh)
    !---------------------------------------------------
    ! TGVECTORS_DEALLOCATION
    ! deallocate all the tangent vectors and local coordinates
    !--------------------------------------------------
    type(cubedsphere):: mesh

    deallocate(mesh%tgx_pu) ! Tangent vector at pu in x direction
    deallocate(mesh%tgy_pu) ! Tangent vector at pu in y direction
    deallocate(mesh%tgx_pv) ! Tangent vector at pv in x direction
    deallocate(mesh%tgy_pv) ! Tangent vector at pv in y direction
    deallocate(mesh%tgx_pc) ! Tangent vector at pc in x direction
    deallocate(mesh%tgy_pc) ! Tangent vector at pc in y direction
  end subroutine tgvectors_deallocation 

subroutine adv_deallocation()
    use advection_vars
    !---------------------------------------------------
    ! ADV_DEALLOCATION
    ! deallocate all the variables of the advection model
    !--------------------------------------------------

    deallocate(Q%f)
    deallocate(Q_exact%f)
    deallocate(div_ugq%f)
    deallocate(div_ugq_exact%f)
    deallocate(div_ugq_error%f)
    deallocate(Qx%f, Qy%f)

    deallocate(px%q_L, px%q_R, px%dq, px%q6, px%f_upw, px%df, px%Q%f)
    deallocate(py%q_L, py%q_R, py%dq, py%q6, py%f_upw, py%df, py%Q%f)

    deallocate(wind_pu%u%f)
    deallocate(wind_pu%v%f)
    deallocate(wind_pu%ucontra%f)
    deallocate(wind_pu%vcontra%f)
    deallocate(wind_pu%ucontra_old%f)
    deallocate(wind_pu%vcontra_old%f)
    deallocate(wind_pu%ucontra_time_av%f)
    deallocate(wind_pu%vcontra_time_av%f)
    deallocate(wind_pu%ucontra_time_centered%f)
    deallocate(wind_pu%vcontra_time_centered%f) 
    deallocate(wind_pu%ucovari%f)
    deallocate(wind_pu%vcovari%f)

    deallocate(wind_pv%u%f)
    deallocate(wind_pv%v%f)
    deallocate(wind_pv%ucontra%f)
    deallocate(wind_pv%vcontra%f)
    deallocate(wind_pv%ucontra_old%f)
    deallocate(wind_pv%vcontra_old%f)
    deallocate(wind_pv%ucontra_time_av%f)
    deallocate(wind_pv%vcontra_time_av%f)
    deallocate(wind_pv%ucontra_time_centered%f)
    deallocate(wind_pv%vcontra_time_centered%f) 
    deallocate(wind_pv%ucovari%f)
    deallocate(wind_pv%vcovari%f)

end subroutine adv_deallocation
end module deallocation 

