/*59:*/
#line 15 "glue.w"

#include "roundoff.h"
/*10:*/
#line 62 "roundoff.w"

#ifndef __sparc__
inline double infinity(){return 1.0/0.0;}
#endif

#ifdef sgi
#include <sys/fpu.h>

void initialize_roundoff(){
union fpc_csr csr;
csr.fc_word= get_fpc_csr();
csr.fc_struct.en_underflow= 1;
set_fpc_csr(csr.fc_word);
}
#else
#ifdef __sparc__
#include <floatingpoint.h>
void initialize_roundoff(){
ieee_handler("set","underflow",SIGFPE_ABORT);
}
#else 
void initialize_roundoff(){
}
#endif 
#endif 


#if defined(sgi) || defined(__sparc__)
int roundoff_ok(){
return 1;
}
#else 
#if defined(AIX)
int roundoff_ok(){
return fp_underflow()==0;
}
#else 
int roundoff_ok(){
return 0;
}
#endif 
#endif 

#line 1 "Codes.w"
/*:10*/
#line 17 "glue.w"


/*:59*/
