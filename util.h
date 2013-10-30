#ifndef _UTIL_H_
#define _UTIL_H_

#include "constant.h"

namespace edubpt {

inline bool is_invalid_value(const Color &col) {
	if (_isnan(col.x) || _isnan(col.y) || _isnan(col.z))
		return true;
	if (col.x < 0.0 || kINF < col.x)
		return true;
	if (col.y < 0.0 || kINF < col.y)
		return true;
	if (col.z < 0.0 || kINF < col.z)
		return true;

	return false;
}
inline bool is_valid_value(const Color &col) {
	return !is_invalid_value(col);
}

inline void clamp(int &x, const int begin, const int end) {
	if (x < begin)
		x = begin;
	if (x >= end)
		x = end - 1;
}

}

#endif