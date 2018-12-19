define(`dim_colons',`:ifelse(`$1',`$2',,`,dim_colons(incr(`$1'),`$2')')')dnl
define(`dim_sizes',`$3($1)ifelse(`$1',`$2',,`,dim_sizes(incr(`$1'),`$2',`$3')')')dnl
define(`dim_slices',`$3($1):$3($1)+$4($1)-1 ifelse(`$1',`$2',,`,dim_slices(incr(`$1'),`$2',`$3',`$4')')')dnl
define(`dim_declare',`ifelse(`$1',`0',,`dimension(dim_colons(`1',`$1')),')')dnl
define(`var_size',`ifelse(`$1',`0',,`size($2), &')')dnl
