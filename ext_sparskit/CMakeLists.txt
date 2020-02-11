set(sources
  BLASSM/blassm.f
  BLASSM/zblassm.f
  BLASSM/matvec.f
  BLASSM/zmatvec.f
  FORMATS/formats.f
  FORMATS/zformats.f
  FORMATS/unary.f
  FORMATS/zunary.f
  FORMATS/sort.f
  INFO/infofun.f
  INFO/zinfofun.f
  INOUT/inout.f
  INOUT/zinout.f
  ORDERINGS/ccncommon.f
  ORDERINGS/ccn.f
  ORDERINGS/zccn.f
  ORDERINGS/zcolor.f
  ORDERINGS/zdsepart.f)

add_library(sparskit_objlib OBJECT ${sources})
set_target_properties(sparskit_objlib PROPERTIES POSITION_INDEPENDENT_CODE ${BUILD_SHARED_LIBS})
