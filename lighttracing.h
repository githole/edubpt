#ifndef _LIGHTTRACING_H_
#define _LIGHTTRACING_H_

#include "random.h"
#include "scene.h"
#include "util.h"
#include "camera.h"
#include "specular.h"
#include "sampling.h"
#include "camera.h"
#include "vertex.h"

namespace edubpt {

struct LighttracingResult {
	Color value;
	int imagebuffer_x, imagebuffer_y;
	bool is_lens_hit;

	LighttracingResult(const Color &value, const int imagebuffer_x, const int imagebuffer_y, const bool is_lens_hit) :
	value(value), imagebuffer_x(imagebuffer_x), imagebuffer_y(imagebuffer_y), is_lens_hit(is_lens_hit) {}
};
LighttracingResult generate_vertices_by_lighttracing(const Camera &camera, Random *rnd, std::vector<Vertex> &vertices) {
	// 光源上にサンプル点生成（y0）
	double pdf_A_on_light;
	const Vec position_on_light = spheres[LightID].position + sample_sphere(spheres[LightID].radius, rnd, &pdf_A_on_light);
	const Vec normal_on_light = normalize(position_on_light - spheres[LightID].position);
	double total_pdf_A = pdf_A_on_light; // 確率密度の積を保持（面積測度に関する確率密度）

	 // 光源上に生成された頂点を頂点リストに追加
	vertices.push_back(Vertex(position_on_light, normal_on_light, normal_on_light, LightID, Vertex::OBJECT_TYPE_LIGHT, total_pdf_A, Color(0, 0, 0)));

	// 現在の放射輝度（モンテカルロ積分のスループット）
	Color MC_throughput = spheres[LightID].emission; 

	// 完全拡散光源を仮定しているので、Diffuse面におけるサンプリング方法と同じものをつかって次の方向を決める
	double now_sampled_pdf_omega;
	const Vec next_dir = sample_hemisphere_cos_term(normal_on_light, rnd, &now_sampled_pdf_omega);
	Ray now_ray(position_on_light, next_dir);
	Vec previous_normal = normal_on_light;
	
	double russian_roulette_probability = 1.0;
	for (;;) {
		Intersection intersection;
		const bool scene_hit = intersect_scene(now_ray, &intersection);
		
		// レンズと交差判定
		Vec position_on_lens, position_on_objectplane, position_on_imagesensor, uv_on_imagesensor;
		const double lens_t = camera.intersect_lens(now_ray, &position_on_lens, &position_on_objectplane, &position_on_imagesensor, &uv_on_imagesensor);
		if (kEPS < lens_t && lens_t < intersection.hitpoint.distance) {
			// レイがレンズにヒット＆イメージセンサにヒット
			const Vec x0_xI = position_on_imagesensor - position_on_lens;
			const Vec x0_xV = position_on_objectplane - position_on_lens;
			const Vec x0_x1 = now_ray.org - position_on_lens;

			// イメージセンサ（ピクセル上）の位置、入射方向によって寄与率が変わる
			// （が、今回はWを定数にしているのでは寄与率は変わらない）
			int x = (int)uv_on_imagesensor.x;
			int y = (int)uv_on_imagesensor.y;
			clamp(x, 0, camera.image_width_px);
			clamp(y, 0, camera.image_height_px);
					
			// レンズの上の点のサンプリング確率を計算
			const double now_sampled_pdf_A = now_sampled_pdf_omega * (dot(normalize(x0_x1), normalize(camera.imagesensor_dir)) / x0_x1.length_squared());
			total_pdf_A *= now_sampled_pdf_A;

			// ジオメトリターム
			const double G =  dot(normalize(x0_x1), normalize(camera.imagesensor_dir)) * dot(normalize(-1.0 * x0_x1), previous_normal) / x0_x1.length_squared();
			MC_throughput = G * MC_throughput;
					
			// レンズ上に生成された点を頂点リストに追加（基本的に使わない）
			vertices.push_back(Vertex(position_on_lens, normalize(camera.imagesensor_dir), normalize(camera.imagesensor_dir), -1, Vertex::OBJECT_TYPE_LENS, total_pdf_A, MC_throughput));
					
			// 幾何的な係数計算 + センサーセンシティビティを計算してモンテカルロ積分して最終的な結果を得る
			const Color result = (camera.W_dash(x0_xV, x0_xI, x0_x1) * MC_throughput) / total_pdf_A;
			return LighttracingResult(result, x, y, true);
		}
		if (!scene_hit)
			break;

		const Sphere &now_object = spheres[intersection.object_id];
		const Hitpoint &hitpoint = intersection.hitpoint;
		const Vec orienting_normal = dot(hitpoint.normal , now_ray.dir) < 0.0 ? hitpoint.normal: (-1.0 * hitpoint.normal); // 交差位置の法線（物体からのレイの入出を考慮）
		russian_roulette_probability = russian_roulette(now_object);
		// ロシアンルーレット
		if (rnd->next01() >= russian_roulette_probability) {
			break;
		}
		
		// 新しい頂点がサンプリングされたので、トータルの確率密度に乗算する
		total_pdf_A *= russian_roulette_probability;
		
		const Vec to_next_vertex = now_ray.org - hitpoint.position;
		// 新しい頂点をサンプリングするための確率密度関数は立体角測度に関するものであったため、これを面積測度に関する確率密度関数に変換する
		const double now_sampled_pdf_A = now_sampled_pdf_omega * (dot(normalize(to_next_vertex), orienting_normal) / to_next_vertex.length_squared());
		// 全ての頂点をサンプリングする確率密度の総計を出す
		total_pdf_A *= now_sampled_pdf_A;

		// ジオメトリターム（G項）
		const double G =  dot(normalize(to_next_vertex), orienting_normal) * dot(normalize(-1.0 * to_next_vertex), previous_normal) / to_next_vertex.length_squared();
		MC_throughput = G * MC_throughput;

		// 新しい頂点を頂点リストに追加する
		vertices.push_back(Vertex(hitpoint.position, orienting_normal, hitpoint.normal, intersection.object_id, 
				now_object.reflection_type == REFLECTION_TYPE_DIFFUSE ? Vertex::OBJECT_TYPE_DIFFUSE :Vertex::OBJECT_TYPE_SPECULAR, 
				total_pdf_A, MC_throughput));


		
		// 材質ごとの処理
		// 次の頂点をサンプリングするために、新しいレイの方向を決める
		switch (now_object.reflection_type) {
		// 完全拡散面
		case REFLECTION_TYPE_DIFFUSE: {
			// 次の方向をサンプリング
			const Vec dir = sample_hemisphere_cos_term(orienting_normal, rnd, &now_sampled_pdf_omega);
			now_ray = Ray(hitpoint.position, dir);
			// スループット調整（dir <- hitpoint.position <- 一個前　の情報がそろっているのでBRDF評価できる)
			MC_throughput = multiply(now_object.color / kPI, MC_throughput);
		} break;
				
		// 完全鏡面
		case REFLECTION_TYPE_SPECULAR: {
			// 完全鏡面なのでレイの反射方向は決定的。
			// スペキュラ面におけるサンプリング確率密度はDiracのδ関数の形をとるが、分母とキャンセルされるため、係数のみ与える
			now_sampled_pdf_omega = 1.0;
			now_ray = Ray(hitpoint.position, reflection_vector(now_ray.dir, hitpoint.normal));

			// スループット調整
			// スペキュラBRDFはδ/cosθになる（立体角測度に関する関数）
			// また、このδはpdf_omegaのδとキャンセルされる（よって、δの扱いは無視してよい。now_sampled_pdf_omega = 1.0にしたのもそのため）
			// なお、pdf_omegaはδになる（立体角測度に関する確率密度）
			// BRDFはcosθをかけて半球積分すると1、pdf_omegaはそのまま半球積分すると1になることからこのようになる。
			MC_throughput = multiply(now_object.color / dot(normalize(to_next_vertex), orienting_normal), MC_throughput);
		} break;
				
		// ガラス
		case REFLECTION_TYPE_REFRACTION: {
			const bool into = dot(hitpoint.normal, orienting_normal) > 0.0; // レイがオブジェクトから出るのか、入るのか

			Vec reflection_dir, refraction_dir;
			double fresnel_reflectance, fresnel_transmittance;
			if (check_refraction(into, hitpoint.position, now_ray.dir, hitpoint.normal, orienting_normal, &reflection_dir, &refraction_dir, &fresnel_reflectance, &fresnel_transmittance)) {
				// 屈折 + 反射
				if (rnd->next01() < reflection_probability) { // 反射
					now_sampled_pdf_omega = 1.0;
					now_ray = Ray(hitpoint.position, reflection_dir);
					MC_throughput = fresnel_reflectance * multiply(now_object.color / dot(normalize(to_next_vertex), orienting_normal), MC_throughput);
					total_pdf_A *= reflection_probability;
				} else { // 屈折
					const double nnt2 = 
						pow(into ?
						refractive_index_of_object / refractive_index_of_vaccum :
						refractive_index_of_vaccum / refractive_index_of_object, 2.0); 
					now_sampled_pdf_omega = 1.0;
					now_ray = Ray(hitpoint.position, refraction_dir);
					MC_throughput = nnt2 * fresnel_transmittance * multiply(now_object.color / dot(normalize(to_next_vertex), orienting_normal), MC_throughput);
					total_pdf_A *= 1.0 - reflection_probability;
				}
			} else {
				// 全反射
				now_sampled_pdf_omega = 1.0;
				now_ray = Ray(hitpoint.position, reflection_dir);
				MC_throughput = multiply(now_object.color / dot(normalize(to_next_vertex), orienting_normal), MC_throughput);
				break;
			}
		} break;

		}

		previous_normal = orienting_normal;
	}
	
	return LighttracingResult(Color(), 0, 0, false);
}


};

#endif