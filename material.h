#ifndef _MATERIAL_H_
#define _MATERIAL_H_

#include "vec.h"

namespace edubpt {

typedef Vec Color;

enum ReflectionType {
	REFLECTION_TYPE_DIFFUSE,	// 完全拡散面。いわゆるLambertian面。
	REFLECTION_TYPE_MIRROR,		// 理想的な鏡面。
	REFLECTION_TYPE_GLASS,		// 理想的なガラス的物質。
};

const double refractive_index_of_vaccum = 1.0; // 真空の屈折率
const double refractive_index_of_object = 1.5; // オブジェクトの屈折率
const Color BackgroundColor(0, 0, 0);

};

#endif
