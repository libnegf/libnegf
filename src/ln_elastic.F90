module ln_elastic
  use interactions
  implicit none
  private

  public :: TElastic

  type, abstract, extends(TInteraction) :: TElastic

  end type TElastic


end module ln_elastic

