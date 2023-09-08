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
    deallocate(mesh%ll2contra_pc) ! latlon 2 contravariant conversion at c
    deallocate(mesh%ll2contra_po) ! latlon 2 contravariant conversion at o
    deallocate(mesh%contra2ll_pu) ! contravariant 2 latlon conversion at u
    deallocate(mesh%contra2ll_pv) ! contravariant 2 latlon conversion at v 
    deallocate(mesh%contra2ll_pc) ! contravariant 2 latlon conversion at c
    deallocate(mesh%contra2ll_po) ! latlon 2 contravariant conversion at o

    deallocate(mesh%ll2covari_pu) ! latlon 2 covariant conversion at u
    deallocate(mesh%ll2covari_pv) ! latlon 2 covariant conversion at v
    deallocate(mesh%ll2covari_pc) ! latlon 2 covariant conversion at c
    deallocate(mesh%ll2covari_po) ! latlon 2 covariant conversion at o
    deallocate(mesh%covari2ll_pu) ! covariant 2 latlon conversion at u
    deallocate(mesh%covari2ll_pv) ! covariant 2 latlon conversion at v 
    deallocate(mesh%covari2ll_pc) ! covariant 2 latlon conversion at c
    deallocate(mesh%covari2ll_po) ! covariant 2 latlon conversion at o

    deallocate(mesh%covari2contra_pu) ! covari 2 contravariant conversion at u
    deallocate(mesh%covari2contra_pv) ! covari 2 contravariant conversion at v
    deallocate(mesh%covari2contra_pc) ! covari 2 contravariant conversion at c
    deallocate(mesh%covari2contra_po) ! covari 2 contravariant conversion at o
    deallocate(mesh%contra2covari_pu) ! contravariant 2 covari conversion at u
    deallocate(mesh%contra2covari_pv) ! contravariant 2 covari conversion at v 
    deallocate(mesh%contra2covari_pc) ! contravariant 2 covari conversion at c
    deallocate(mesh%contra2covari_po) ! contravariant 2 covari conversion at o

    deallocate(mesh%mt_pc)! Metric tensor at centers
    deallocate(mesh%mt_pu) ! Metric tensor at u
    deallocate(mesh%mt_pv) ! Metric tensor at v
    deallocate(mesh%mt_po) ! Metric tensor at o
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

    deallocate(L_pc%y_support, L_pc%f_support, L_pc%x_nodes, L_pc%y_nodes)
    deallocate(L_pc%p_nodes, L_pc%f_nodes, L_pc%k0, L_pc%kend) 
    deallocate(L_pc%halodata_east, L_pc%halodata_west, L_pc%halodata_north, L_pc%halodata_south) 


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

    deallocate(wind_pc%u%f)
    deallocate(wind_pc%v%f)
    deallocate(wind_pc%ucontra%f)
    deallocate(wind_pc%vcontra%f)
    deallocate(wind_pc%ucontra_old%f)
    deallocate(wind_pc%vcontra_old%f)
    deallocate(wind_pc%ucontra_time_av%f)
    deallocate(wind_pc%vcontra_time_av%f)
    deallocate(wind_pc%ucontra_time_centered%f)
    deallocate(wind_pc%vcontra_time_centered%f) 
    deallocate(wind_pc%ucovari%f)
    deallocate(wind_pc%vcovari%f)
    deallocate(wind_pc%ucovari_old%f)
    deallocate(wind_pc%vcovari_old%f)

    !deallocate(wind_po%u%f)
    !deallocate(wind_po%v%f)
    !deallocate(wind_po%ucontra%f)
    !deallocate(wind_po%vcontra%f)
    !deallocate(wind_po%ucontra_old%f)
    !deallocate(wind_po%vcontra_old%f)
    !deallocate(wind_po%ucontra_time_av%f)
    !deallocate(wind_po%vcontra_time_av%f)
    !deallocate(wind_po%ucontra_time_centered%f)
    !deallocate(wind_po%vcontra_time_centered%f) 
    !deallocate(wind_po%ucovari%f)
    !deallocate(wind_po%vcovari%f)
    !deallocate(wind_po%ucovari_old%f)
    !deallocate(wind_po%vcovari_old%f)



end subroutine adv_deallocation


subroutine swm_deallocation()
    use swm_vars
    !---------------------------------------------------
    ! SWM_DEALLOCATION
    ! deallocate all the variables of the shallow water model
    !--------------------------------------------------

    deallocate(H%f)
    deallocate(H_exact%f)
    deallocate(div_ugH%f, rel_vort%f, abs_vort%f, fcoriolis_pc%f, div_abs_vort%f)
    deallocate(dx_H_pv%f, dy_H_pu%f)
    deallocate(dx_div_ugH_pv%f, dy_div_ugH_pu%f)
    deallocate(dy_K_pu%f, dx_K_pv%f)
    deallocate(Qx%f, Qy%f)

    deallocate(px%q_L, px%q_R, px%dq, px%q6, px%f_upw, px%df, px%Q%f)
    deallocate(py%q_L, py%q_R, py%dq, py%q6, py%f_upw, py%df, py%Q%f)

    deallocate(Ku_px%q_L, Ku_px%q_R, Ku_px%dq, Ku_px%q6, Ku_px%f_upw, Ku_px%df, Ku_px%Q%f)
    deallocate(Kv_py%q_L, Kv_py%q_R, Kv_py%dq, Kv_py%q6, Kv_py%f_upw, Kv_py%df, Kv_py%Q%f)

    deallocate(L_pc%y_support, L_pc%f_support, L_pc%x_nodes, L_pc%y_nodes)
    deallocate(L_pc%p_nodes, L_pc%f_nodes, L_pc%k0, L_pc%kend) 
    deallocate(L_pc%halodata_east, L_pc%halodata_west, L_pc%halodata_north, L_pc%halodata_south) 

    if(swm_simul%ic==0)then
        deallocate(div_ugH_exact%f)
        deallocate(div_ugH_error%f)
        deallocate(div_ugH_pu%f)
        deallocate(div_ugH_pv%f)
        deallocate(div_ugH_po%f)

        deallocate(rel_vort_exact%f)
        deallocate(rel_vort_error%f)

        deallocate(abs_vort_exact%f)
        deallocate(abs_vort_error%f)
        deallocate(abs_vort_flux_exact_pu%f)
        deallocate(abs_vort_flux_error_pu%f)
        deallocate(abs_vort_flux_exact_pv%f)
        deallocate(abs_vort_flux_error_pv%f)

        deallocate(H_po_exact%f)
        deallocate(H_pu_exact%f)
        deallocate(H_pv_exact%f)

        deallocate(Ku_po_exact%f)
        deallocate(Kv_po_exact%f)
        deallocate(K_po_exact%f)

    end if

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

    deallocate(wind_pc%u%f)
    deallocate(wind_pc%v%f)
    deallocate(wind_pc%ucontra%f)
    deallocate(wind_pc%vcontra%f)
    deallocate(wind_pc%ucontra_old%f)
    deallocate(wind_pc%vcontra_old%f)
    deallocate(wind_pc%ucontra_time_av%f)
    deallocate(wind_pc%vcontra_time_av%f)
    deallocate(wind_pc%ucontra_time_centered%f)
    deallocate(wind_pc%vcontra_time_centered%f) 
    deallocate(wind_pc%ucovari%f)
    deallocate(wind_pc%vcovari%f)
    deallocate(wind_pc%ucovari_old%f)
    deallocate(wind_pc%vcovari_old%f)

    deallocate(wind_po%u%f)
    deallocate(wind_po%v%f)
    deallocate(wind_po%ucontra%f)
    deallocate(wind_po%vcontra%f)
    deallocate(wind_po%ucontra_old%f)
    deallocate(wind_po%vcontra_old%f)
    deallocate(wind_po%ucontra_time_av%f)
    deallocate(wind_po%vcontra_time_av%f)
    deallocate(wind_po%ucontra_time_centered%f)
    deallocate(wind_po%vcontra_time_centered%f) 
    deallocate(wind_po%ucovari%f)
    deallocate(wind_po%vcovari%f)
    deallocate(wind_po%ucovari_old%f)
    deallocate(wind_po%vcovari_old%f)

    deallocate(ucovari_pv_exact%f)
    deallocate(vcovari_pu_exact%f)
    deallocate(ucovari_pv_error%f)
    deallocate(vcovari_pu_error%f)

end subroutine swm_deallocation


end module deallocation 

