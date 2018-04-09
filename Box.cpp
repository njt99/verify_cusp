#include "Box.h"
#include <stdlib.h>

double scale[6];
static bool scale_initialized = false; 
Box::Box(char * where) {
	if (!scale_initialized) {
		scale_initialized = true;
		for (int i = 0; i < 6; ++i) {
			scale[i] = pow(2, -i / 6.0);
		}
	}
	for (int i = 0; i < 6; ++i) {
		center_digits[i] = 0;
		size_digits[i] = 8;
	}
    int pos = 0;
    int idx = 0;
    int dir;
    while (where[idx] != '\0') {
        if (where[idx] == '0') {
            dir = 0;
        } else if (where[idx] == '1') {
            dir = 1;
        } else {
            fprintf(stderr, "verify: fatal error at %s\n", where);
            exit(1);
        }
	    size_digits[pos] *= 0.5;
	    center_digits[pos] += (2*dir-1)*size_digits[pos];
	    ++pos;
        if (pos == 6) {
            pos = 0;
        }
	    ++idx;
    }
    compute_center_and_size();
}

void Box::compute_center_and_size()
{
	for (int i = 0; i < 6; ++i) {
        // GMT paper page 419 of Annals
        // box_size guarantees that :
        // box_center - box_size <= true_center - true_size
        // box_center + box_size >= true_center + true_size
        // where box operations are floating point. 
        box_center[i] = scale[i]*center_digits[i];
        box_size[i]= (1+2*EPS)*(size_digits[i]*scale[i]+HALFEPS*fabs(center_digits[i]));
    }
}

Params<ACJ> Box::cover() const
{
	Params<ACJ> result;
	result.lattice = ACJ(
		XComplex(box_center[3], box_center[0]),
		XComplex(box_size[3], box_size[0]),
		0.,
		0.
	);
	result.loxodromic_sqrt = ACJ(
		XComplex(box_center[4], box_center[1]),
		0.,
		XComplex(box_size[4], box_size[1]),
		0.
	);
	result.parabolic = ACJ(
		XComplex(box_center[5], box_center[2]),
		0.,
		0.,
		XComplex(box_size[5], box_size[2])
	);
	return result;
}

Params<XComplex> Box::nearest() const
{
	double m[6];
    double temp;
	for (int i = 0; i < 6; ++i) {
		if (center_digits[i] < 0) {
            // GMT paper page 419 of Annals
            // temp is guaranteed to be >= than true_center + true_size
            temp = box_center[i]+box_size[i];
            if (temp > 0 ) { // Check if we have overlapped 0
                m[i] = 0;
            } else { // We know scale is positive
		    	m[i] = temp;
            }  
        } else {
            // GMT paper page 419 of Annals
            // temp is guaranteed to be <= than true_center - true_size
            temp = box_center[i]-box_size[i];
            if (temp < 0 ) { // Check if we have overlapped 0
                m[i] = 0;
            } else { // We know scale is positive
		    	m[i] = temp;
            }  
        }
	}
	
	Params<XComplex> result;
	result.lattice = XComplex(m[3], m[0]);
	result.loxodromic_sqrt = XComplex(m[4], m[1]);
	result.parabolic = XComplex(m[5], m[2]);
	return result;
}

Params<XComplex> Box::furthest() const
{
	double m[6];
	for (int i = 0; i < 6; ++i) {
		if (center_digits[i] < 0) {
            // GMT paper page 419 of Annals
            // guaranteed to be <= than true_center - true_size
		    m[i] = box_center[i]-box_size[i];
        } else {
            // guaranteed to be >= than true_center + true_size
		    m[i] = box_center[i]+box_size[i];
        }
	}
	
	Params<XComplex> result;
	result.lattice = XComplex(m[3], m[0]);
	result.loxodromic_sqrt = XComplex(m[4], m[1]);
	result.parabolic = XComplex(m[5], m[2]);
	return result;
}

// Note: size_digits is always positive
Params<XComplex> Box::maximum() const
{
	double m[6];
	for (int i = 0; i < 6; ++i) {
        // GMT paper page 419 of Annals
        // guaranteed to be >= than true_center + true_size
        m[i] = box_center[i]+box_size[i];
	}
	
	Params<XComplex> result;
	result.lattice = XComplex(m[3], m[0]);
	result.loxodromic_sqrt = XComplex(m[4], m[1]);
	result.parabolic = XComplex(m[5], m[2]);
	return result;
}
