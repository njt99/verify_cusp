#ifndef __Params_h
#define __Params_h
#include <math.h>
#include "SL2ACJ.h"

template<class N> struct Params {
	N lattice;
	N loxodromic_sqrt;
	N parabolic;
};

SL2ACJ construct_G(const Params<ACJ>& params);
SL2ACJ construct_T(const Params<ACJ>& params, int x, int y);
SL2ACJ construct_word(const Params<ACJ>& params, char* word);

#endif // __Params_h
