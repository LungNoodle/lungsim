module species_transport_c
  implicit none
  private

contains

!!!###################################################################################

  subroutine initialise_transport_c() bind(C, name="initialise_transport_c")

    use species_transport, only: initialise_transport
    implicit none

    call initialise_transport

  end subroutine initialise_transport_c


end module species_transport_c
