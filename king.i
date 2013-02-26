%module king
/* %include "typemaps.i" */
%include "cpointer.i"
%{
#include "include/king/cte_phys.h"
#include "include/king/king.h"
#include "include/king/mod.h"
/*
#include "include/king/rand.h"
#include "include/king/rk4.h"
#include "include/king/utils.h"
*/
%}

%include "include/king/cte_phys.h"
%include "include/king/king.h"
%include "include/king/mod.h"
/*
%include "include/king/rand.h"
%include "include/king/rk4.h"
%include "include/king/utils.h"
*/

/* %apply double *OUTPUT { double *G }; */
/* %pointer_functions(double, doublep); */
/*%pointer_class(double, doublec);

%extend doublec {
        char *__str__(void)
        {
                static char tmp[1024];
                sprintf(tmp, "%.15g", doublec_value($self));
                return tmp;
        }
}
*/
%extend King {
        void __del__(void)
        {
                King_free($self);
        }
}

