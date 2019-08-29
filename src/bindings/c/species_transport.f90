module species_transport_c
  implicit none
  private

contains

!!!###################################################################################

  subroutine initialise_transport_c() bind(C, name="initialise_transport_c")

    use species_transport, only: initialise_transport
    implicit none

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_initialise_transport
#else
    call initialise_transport
#endif

  end subroutine initialise_transport_c


end module species_transport_c
