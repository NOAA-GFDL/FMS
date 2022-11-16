#undef  FMS_MO_KIND_
#define FMS_MO_KIND_ r8_kind

#undef  GET_MOSAIC_XGRID_
#define GET_MOSAIC_XGRID_ get_mosaic_xgrid_r8

#undef  CALC_MOSAIC_GRID_AREA_
#define CALC_MOSAIC_GRID_AREA_ calc_mosaic_grid_area_r8

#undef  CALC_MOSAIC_GRID_GREAT_CIRCLE_AREA_
#define CALC_MOSAIC_GRID_GREAT_CIRCLE_AREA_ calc_mosaic_grid_great_circle_area_r8

#undef  IS_INSIDE_POLYGON_
#define IS_INSIDE_POLYGON_ is_inside_polygon_r8

#include <mosaic2.fh>
