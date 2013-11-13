#ifndef _SPECULAR_H_
#define _SPECULAR_H_

#include "material.h"
#include "vec.h"

namespace edubpt {

// 界面に向かう入射ベクトルがin
// 入射点における法線がnormal
// 返り値は入射点から出射する、反射ベクトル
Vec reflection_vector(const Vec &in, const Vec &normal) {
	return normalize(in - normal * 2.0 * dot(normal, in));
}

// falseなら全反射
// trueなら反射＋屈折
bool check_refraction(
	const bool into,
	const Vec& position,
	const Vec &in,
	const Vec &normal,
	const Vec &orienting_normal,
	Vec *reflection_dir,
	Vec *refraction_dir,
	double *fresnel_reflectance,
	double *fresnel_transmittance) {
	*reflection_dir = reflection_vector(in, normal);

	// Snellの法則
	const double nnt = into ? refractive_index_of_vaccum / refractive_index_of_object : refractive_index_of_object / refractive_index_of_vaccum;
	const double ddn = dot(in, orienting_normal);
	const double cos2t = 1.0 - nnt * nnt * (1.0 - ddn * ddn);
		
	if (cos2t < 0.0) { // 全反射
		*refraction_dir = Vec();
		*fresnel_reflectance = 1.0;
		*fresnel_transmittance = 0.0;
		return false;
	}
	// 屈折の方向
	*refraction_dir = 
		normalize(in * nnt - normal * (into ? 1.0 : -1.0) * (ddn * nnt + sqrt(cos2t)));
		
	// SchlickによるFresnelの反射係数の近似を使う
	const double a = refractive_index_of_object - refractive_index_of_vaccum;
	const double b = refractive_index_of_object + refractive_index_of_vaccum;
	const double R0 = (a * a) / (b * b);

	const double c = 1.0 - (into ? -ddn : dot(*refraction_dir, -1.0 * orienting_normal));
	*fresnel_reflectance = R0 + (1.0 - R0) * pow(c, 5.0); // 反射方向の光が反射してray.dirの方向に運ぶ割合。同時に屈折方向の光が反射する方向に運ぶ割合。
	*fresnel_transmittance = 1.0 - *fresnel_reflectance; // 屈折方向の光が屈折してray.dirの方向に運ぶ割合

	return true;
}

};

#endif

