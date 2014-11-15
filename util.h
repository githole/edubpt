#ifndef _UTIL_H_
#define _UTIL_H_

#include <algorithm>
#include "constant.h"

namespace edubpt {

#if __STDC_VERSION__ >= 199901L
	// isnanが定義されていると予想
#else
	// isnanが定義されていないと予想
#if defined(_MSC_VER)
	inline int isnan(double arg) {
		return _isnan(arg);
	}
#else
	// C99コンパイラでもなく、MSVCでもないなら諦める
#endif // defined(_MSC_VER)

#endif // __STDC_VERSION__ >= 199901L

inline bool is_invalid_value(const Color &col) {
	if (isnan(col.x) || isnan(col.y) || isnan(col.z))
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