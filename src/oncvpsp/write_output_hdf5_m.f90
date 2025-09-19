module write_output_hdf5_m
   use hdf5_utils, only: hdf_open_file, hdf_close_file, &
      hdf_create_group, hdf_open_group, hdf_close_group, &
      hdf_create_dataset, hdf_write_dataset, hdf_write_attribute, &
      hdf_write_vector_to_dataset, &
      HID_T
   use precision_m, only: dp
   use lsch_m, only: lschfb, lschfs, lschvkbb, lschvkbs
   implicit none
   private
   public :: write_output_hdf5
contains

subroutine write_output_hdf5(filename, lmax, npa, epa, lloc, irc, &
                             vkb, evkb, nproj, rr, vfull, vp, vpuns, zz, mmax, mxprj, drl, nrl, &
                             rho, rhoc, rhomod, srel, cvgplt)
   ! Input variables
   !> HDF5 filename
   character(*) :: filename
   !> Maximum angular momentum
   integer :: lmax
   !> Local potential angular momentum
   integer :: lloc
   !> Number of points in logarithmic radial mesh
   integer :: mmax
   !> Maximum number of projectors
   integer :: mxprj
   !> Number of points in linear radial mesh
   integer :: nrl
   !> Principal quantum number for corresponding all-electron state
   integer :: npa(mxprj, 6)
   !> Indices of core radii
   integer :: irc(6)
   !> Number of VKB projectors for each l
   integer :: nproj(6)
   !> Atomic number
   real(dp) :: zz
   !> Spacing of linear radial mesh
   real(dp) :: drl
   !> Logarithmic radial grid
   real(dp) :: rr(mmax)
   !> Semi-local pseudopotentials
   real(dp) :: vp(mmax, 5)
   !> Unscreened pseudopotentials
   real(dp) :: vpuns(mmax, 5)
   !> All-electron potential
   real(dp) :: vfull(mmax)
   !> VKB projectors
   real(dp) :: vkb(mmax, mxprj, 4)
   !> Valence pseudocharge
   real(dp) :: rho(mmax)
   !> Core charge
   real(dp) :: rhoc(mmax)
   !> Model core charge
   real(dp) :: rhomod(mmax, 5)
   !> Bound-state or scattering state reference energies for VKB potentials
   real(dp) :: epa(mxprj, 6)
   !> Coefficients of VKB projectors
   real(dp) :: evkb(mxprj, 4)
   !> Energy per electron error vs. cutoff
   real(dp) :: cvgplt(2, 7, mxprj, 4)
   !> Scalar-relativistic flag
   logical :: srel

   ! Local variables
   !> Angular momentum
   integer :: ll
   !> Angular momentum index (l + 1)
   integer :: l1
   !> Loop index
   integer :: ii
   !> Loop index
   integer :: jj
   !> Error code
   integer :: ierr
   !> Matching index (dummy)
   integer :: mch
   !> ???
   integer :: n2
   !> ???
   integer :: nn
   !> Projector index
   integer :: iprj
   !> Number of projectors for given l
   integer :: nproj_l
   !> Logarithmic radial mesh spacing (log(amesh))
   real(dp) :: al
   !> Maximum energy for finding pseudo wavefunctions
   real(dp) :: emax
   !> Minimum energy for finding pseudo wavefunctions
   real(dp) :: emin
   !> Test energy for bound or scattering state
   real(dp) :: etest
   !> Sign of all-electron wavefunction at some radial point
   real(dp) :: sgnae
   !> Sign of pseudo wavefunction at some radial point
   real(dp) :: sgnps
   !> All-electron wavefunction
   real(dp), allocatable :: uae(:)
   !> Pseudo wavefunction
   real(dp), allocatable :: ups(:)
   !> Wavefunction derivative (dummy)
   real(dp), allocatable :: up(:)
   !> Erro message
   character(len=1024) :: error_message

   ! HDF5 variables
   !> HDF5 file identifier
   integer(HID_T) :: file_id
   !> HDF5 group and subgroup identifier
   integer(HID_T) :: loc_id
   !> HDF5 Bound wavefunction group identifier
   integer(HID_T) :: bound_wf_group_id
   !> HDF5 Scattering wavefunction group identifier
   integer(HID_T) :: scattering_wf_group_id
   !> HDF5 angular momentum subgroup identifier
   integer(HID_T) :: l_subgroup_id
   !> HDF5 principal quantum number or projector index subgroup identifier
   integer(HID_T) :: n_iproj_subgroup_id
   !> HDF5 VKB projector group identifier
   integer(HID_T) :: vkb_group_id

   ! For creating dataset / group names
   character(len=1024) :: l_group_name
   character(len=1024) :: n_iproj_group_name

   allocate (uae(mmax), ups(mmax), up(mmax))

   al = 0.01_dp * log(rr(101) / rr(1))
   n2 = int(log(dble(nrl) * drl / rr(1)) / al + 1.0_dp)

   ! Create HDF5 file
   call hdf_open_file(file_id, filename, 'REPLACE', 'WRITE')
   ! Logarithmic radial grid
   call hdf_write_dataset(file_id, 'logarithmic_mesh', rr)
   ! Valence pseudo charge density
   call hdf_write_dataset(file_id, 'ps_valence_charge_density', rho)
   ! Unscreened pseudopotentials
   call hdf_write_dataset(file_id, 'unscreened_pseudopotentials', vpuns)
   ! Local potential
   call hdf_write_dataset(file_id, 'local_potential', vpuns(:, lloc + 1))
   ! Core charge density
   call hdf_write_dataset(file_id, 'ae_core_charge_density', rhoc)
   ! Model core charge density
   call hdf_write_dataset(file_id, 'model_core_charge_density', rhomod)
   ! Bound and scattering all-electron and pseudo wavefunctions for each n / iprj:
   ! /
   !  bound_wavefunctions/
   !    l_0/
   !      n_1/
   !        ae_wavefunction(mmax)
   !        ps_wavefunction(mmax)
   !      ...
   !    ...
   !  scattering_wavefunctions/
   !    l_0/
   !      iprj_1/
   !        ae_wavefunction(mmax)
   !        ps_wavefunction(mmax)
   !      ...
   !    ...
   call hdf_create_group(file_id, 'bound_wavefunctions')
   call hdf_open_group(file_id, 'bound_wavefunctions', bound_wf_group_id)
   call hdf_create_group(file_id, 'scattering_wavefunctions')
   call hdf_open_group(file_id, 'scattering_wavefunctions', scattering_wf_group_id)
   do l1 = 1, lmax + 1
      ll = l1 - 1
      nproj_l = nproj(l1)
      if (ll == lloc) nproj_l = 0
      l_group_name = ''
      write(l_group_name, '(a,i0)') 'l_', ll
      call hdf_create_group(bound_wf_group_id, l_group_name)
      call hdf_create_group(scattering_wf_group_id, l_group_name)
      do iprj = 1, nproj(l1)
         if (epa(iprj, l1) < 0.0_dp) then
            ! n= xx, l= xx, all-electron wavefunction, pseudo wavefunction
            etest = epa(iprj, l1)
            call lschfb(npa(iprj, l1), ll, ierr, etest, rr, vfull, uae, up, zz, mmax, mch, srel)
            if (ierr /= 0) then
               write(error_message, '(a,i0,a,i0)') 'write_output_hdf5: ERROR in lschfb for l=', ll, ' n=', npa(iprj, l1)
               call hdf_close_group(bound_wf_group_id)
               call hdf_close_group(scattering_wf_group_id)
               call hdf_close_file(file_id)
               stop error_message
            end if
            emax = 0.9_dp * etest
            emin = 1.1_dp * etest
            call lschvkbb(ll + iprj, ll, nproj_l, ierr, etest, emin, emax, rr, vp(:, lloc + 1), vkb(:, :, l1), &
                          evkb(1, l1), ups, up, mmax, mch)
            if (ierr /= 0) then
               write(error_message, '(a,i0,a,i0)') 'write_output_hdf5: ERROR in lschvkbb for l=', ll, ' iprj=', iprj
               call hdf_close_group(bound_wf_group_id)
               call hdf_close_group(scattering_wf_group_id)
               call hdf_close_file(file_id)
               stop error_message
            end if
            call hdf_open_group(bound_wf_group_id, l_group_name, l_subgroup_id)
            n_iproj_group_name = ''
            write(n_iproj_group_name, '(a,i0)') 'n_', npa(iprj, l1)
            call hdf_create_group(l_subgroup_id, n_iproj_group_name)
            call hdf_open_group(l_subgroup_id, n_iproj_group_name, n_iproj_subgroup_id)
            call hdf_write_dataset(n_iproj_subgroup_id, 'ae_wavefunction', uae)
            call hdf_write_dataset(n_iproj_subgroup_id, 'ps_wavefunction', ups)
            call hdf_close_group(n_iproj_subgroup_id)
            call hdf_close_group(l_subgroup_id)
         else
            ! scattering, iprj= xx, l= xx, all-electron wavefunction, pseudo wavefunction
            etest = epa(iprj, l1)
            call lschfs(nn, ll, ierr, etest, rr, vfull, uae, up, zz, mmax, n2, srel)
            if (ierr /= 0) then
               write(error_message, '(a,i0,a,i0)') 'write_output_hdf5: ERROR in lschfs for l=', ll, ' iprj=', iprj
               call hdf_close_group(bound_wf_group_id)
               call hdf_close_group(scattering_wf_group_id)
               call hdf_close_file(file_id)
               stop error_message
            end if
            call lschvkbs(ll, nproj_l, etest, rr, vp(:, lloc + 1), vkb(:, :, l1), evkb(:, l1), ups, up, mmax, n2)
            if (ierr /= 0) then
               write(error_message, '(a,i0,a,i0)') 'write_output_hdf5: ERROR in lschvkbs for l=', ll, ' iprj=', iprj
               call hdf_close_group(bound_wf_group_id)
               call hdf_close_group(scattering_wf_group_id)
               call hdf_close_file(file_id)
               stop error_message
            end if
            sgnae = sign(1.0_dp, uae(n2))
            sgnps = sign(1.0_dp, ups(n2))
            do ii = 1, mmax
               uae(ii) = sgnae * uae(ii)
               ups(ii) = sgnps * ups(ii)
            end do
            n_iproj_group_name = ''
            write(n_iproj_group_name, '(a,i0)') 'iprj_', iprj
            call hdf_open_group(scattering_wf_group_id, l_group_name, l_subgroup_id)
            call hdf_create_group(l_subgroup_id, n_iproj_group_name)
            call hdf_open_group(l_subgroup_id, n_iproj_group_name, n_iproj_subgroup_id)
            call hdf_write_dataset(n_iproj_subgroup_id, 'ae_wavefunction', uae)
            call hdf_write_dataset(n_iproj_subgroup_id, 'ps_wavefunction', ups)
         end if
      end do  ! iprj
   end do  ! l1
   call hdf_close_group(bound_wf_group_id)
   call hdf_close_group(scattering_wf_group_id)
   ! VKB projectors:
   ! /vkb_projectors/
   !    l_0/
   !      vkb(mmax, nproj(l))
   !    ...
   call hdf_create_group(file_id, 'vkb_projectors')
   call hdf_open_group(file_id, 'vkb_projectors', vkb_group_id)
   do l1 = 1, lmax + 1
      ll = l1 - 1
      l_group_name = ''
      write(l_group_name, '(a,i0)') 'l_', ll
      call hdf_create_group(vkb_group_id, l_group_name)
      call hdf_open_group(vkb_group_id, l_group_name, l_subgroup_id)
      call hdf_write_dataset(l_subgroup_id, 'vkb', vkb(:, 1:nproj(l1), l1))
      call hdf_close_group(l_subgroup_id)
   end do
   call hdf_close_group(vkb_group_id)
   ! Write convergence profiles:
   call hdf_write_dataset(file_id, 'convergence_profiles', cvgplt)
   ! Close HDF5 file
   call hdf_close_file(file_id)

end subroutine write_output_hdf5

end module write_output_hdf5_m
