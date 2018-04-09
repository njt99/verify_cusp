#ifndef __Box_h
#define __Box_h
#include "Params.h"

struct Box {
	Box(char* where);
	Params<ACJ> cover() const;
	Params<XComplex> nearest() const; // returns closest to 0 or 0 if box overlaps
	Params<XComplex> furthest() const; // furthest from 0
	Params<XComplex> maximum() const; // maximizes all values
private:
	double center_digits[6];
	double size_digits[6];
    double box_center[6];
    double box_size[6];
    void compute_center_and_size();
};

#endif // __Box_h
