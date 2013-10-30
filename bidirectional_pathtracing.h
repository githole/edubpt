#ifndef _BIDIRECTIONAL_PATHTRACING_H_
#define _BIDIRECTIONAL_PATHTRACING_H_

#include "hdr.h"
#include "random.h"
#include "scene.h"
#include "util.h"
#include "camera.h"
#include "specular.h"
#include "sampling.h"
#include "camera.h"
#include "vertex.h"

#include "pathtracing.h"
#include "lighttracing.h"

#include <omp.h>
#include <assert.h>

namespace edubpt {

//
// 双方向パストレーシング
//
double calc_pdf_A_by_pathtracing(const Camera &camera, const Vertex &xip1, const Vertex &xi) {
	const Vec to = xi.position - xip1.position;
	const Vec normalized_to = normalize(to);
	double coeff = 0.0;

	// 始点のオブジェクトの種類によって、その次の頂点のサンプリング確率の計算方法が変わる
	switch (xip1.type) {
	case Vertex::OBJECT_TYPE_LIGHT:
	case Vertex::OBJECT_TYPE_DIFFUSE:
		coeff = sample_hemisphere_pdf_omega(xip1.orienting_normal, normalized_to);
		break;
	case Vertex::OBJECT_TYPE_LENS:
		{
			const Ray test_ray(xi.position, -1.0 * normalized_to);
			Vec position_on_lens, position_on_objectplane, position_on_imagesensor, uv_on_imagesensor;
			const double lens_t = camera.intersect_lens(test_ray, &position_on_lens, &position_on_objectplane, &position_on_imagesensor, &uv_on_imagesensor);
			if (kEPS < lens_t) {
				const Vec x0_xI = position_on_imagesensor - position_on_lens;
				const Vec x0_xV = position_on_objectplane - position_on_lens;
				const Vec x0_x1 = test_ray.org - position_on_lens;
					
				const double P_Image = 1.0 / (camera.imagesensor_width * camera.imagesensor_height);
				const double PA_x1 = camera.P_Image_to_PA_x1(P_Image, x0_xV, x0_x1, xi.orienting_normal);
				return PA_x1;
			} 
			return 0.0;
		}
		break;
	case Vertex::OBJECT_TYPE_SPECULAR:
		if (spheres[xip1.objectID].reflection_type == REFLECTION_TYPE_REFRACTION) {
			// 屈折か反射かを判定
			const double cost = dot(xip1.object_normal, normalize(to));
			if (cost > 0.0) // 反射
				coeff = reflection_probability;
			else
				coeff = 1.0 - reflection_probability;
		} else
			coeff = 1.0;
		break;
	}

	double cost = dot(-1.0 * normalized_to, xi.orienting_normal);
	if (cost < 0.0)
		cost = dot(-1.0 * normalized_to, -1.0 * xi.orienting_normal);
	return coeff * (cost / to.length_squared());
}
double calc_pdf_A_by_lighttracing(const Vertex &xim1, const Vertex &xi) {
	const Vec to = xi.position - xim1.position;
	const Vec normalized_to = normalize(to);
	double coeff = 0.0;
	
	switch (xim1.type) {
	case Vertex::OBJECT_TYPE_LIGHT:
	case Vertex::OBJECT_TYPE_DIFFUSE:
		coeff = sample_hemisphere_pdf_omega(xim1.orienting_normal, normalized_to);
		break;
	case Vertex::OBJECT_TYPE_LENS:
		assert(false);
		break;
	case Vertex::OBJECT_TYPE_SPECULAR:
		if (spheres[xim1.objectID].reflection_type == REFLECTION_TYPE_REFRACTION) {
			// 屈折か反射かを判定
			const double cost = dot(xim1.object_normal, normalize(to));
			if (cost > 0.0) // 反射
				coeff = reflection_probability;
			else
				coeff = 1.0 - reflection_probability;
		} else
			coeff = 1.0;
		break;
	}

	double cost = dot(-1.0 * normalized_to, xi.orienting_normal);
	if (cost < 0.0)
		cost = dot(-1.0 * normalized_to, -1.0 * xi.orienting_normal);
	return  coeff * (cost / to.length_squared());
}
double calc_mis_weight(const Camera &camera, const double total_pdf_A, const std::vector<Vertex> &eye_vs, const std::vector<Vertex> &light_vs, const int num_eye_vertex, const int num_light_vertex) {
	std::vector<const Vertex*> vs(num_light_vertex + num_eye_vertex);
	std::vector<double> pi1_pi(num_eye_vertex + num_light_vertex);
	const double PA_y0 = 1.0 / (kPI * spheres[LightID].radius * spheres[LightID].radius); // 光源上のサンプリング確率
	const double PA_x0 = camera.sampling_pdf_on_lens(); // レンズ上のサンプリング確率

	// 頂点を一列に並べる
	const int k = num_eye_vertex + num_light_vertex - 1;
	for (int i = 0; i < num_light_vertex; ++i)
		vs[i] = &light_vs[i];
	for (int i = num_eye_vertex - 1; i >= 0; --i)
		vs[num_light_vertex + num_eye_vertex - 1 - i] = &eye_vs[i];

	// ロシアンルーレットの確率忘れずに！
	pi1_pi[0] = PA_y0 / (calc_pdf_A_by_pathtracing(camera, *vs[1], *vs[0]) * russian_roulette(spheres[vs[0]->objectID]));
	for (int i = 1; i < k; ++i) {
		const double a = calc_pdf_A_by_lighttracing(*vs[i - 1], *vs[i]);
		const double b = calc_pdf_A_by_pathtracing(camera, *vs[i + 1], *vs[i]);
		pi1_pi[i] = a / b;
	}
	pi1_pi[k] = (calc_pdf_A_by_lighttracing(*vs[k - 1], *vs[k]) * russian_roulette(spheres[vs[k]->objectID])) / PA_x0;

	// pを求める
	std::vector<double> p(num_eye_vertex + num_light_vertex + 1);
	p[num_light_vertex] = total_pdf_A;
	for (int i = num_light_vertex; i <= k; ++i) {
		p[i + 1] = p[i] * pi1_pi[i];
	}
	for (int i = num_light_vertex - 1; i >= 0; --i) {
		p[i] = p[i + 1] / pi1_pi[i];
	}

	// Specular処理
	for (int i = 0; i < vs.size(); ++i) {
		if (vs[i]->type == Vertex::OBJECT_TYPE_SPECULAR) {
			p[i] = 0.0;
			p[i + 1] = 0.0;
		}
	}

	// Power-heuristic
	const double beta = 2.0;
	double mis_weight = 0.0;
	for (int i = 0; i < p.size(); ++i) {
		mis_weight += pow(p[i] / p[num_light_vertex], beta);
//		mis_weight += (p[i] / p[num_light_vertex]) * (p[i] / p[num_light_vertex]);
	}
	return 1.0 / mis_weight;
}

struct BidirectionalPathtracingResult {
	struct Data {
		int imagebuffer_x, imagebuffer_y;
		Color value;
		
		Data(const int imagebuffer_x, const int imagebuffer_y, const Color &value) :
		imagebuffer_x(imagebuffer_x), imagebuffer_y(imagebuffer_y), value(value) {}
	};
	std::vector<Data> data;
};

BidirectionalPathtracingResult bidirectional_pathtracing(const Camera &camera, const int imagebuffer_x, const int imagebuffer_y, Random *rnd) {
	BidirectionalPathtracingResult bpt_result;

	std::vector<Vertex> eye_vs, light_vs;
	PathtracingResult  pt_result = generate_vertices_by_pathtracing(camera, imagebuffer_x, imagebuffer_y, rnd, eye_vs);
	LighttracingResult lt_result = generate_vertices_by_lighttracing(camera, rnd, light_vs);

	// num_light_vertex == 0のとき、カメラ側からのパスが光源にヒットしていれば、それをサンプルとして使用する
	if (pt_result.is_light_hit) {
		const double mis_weight = calc_mis_weight(camera, eye_vs[eye_vs.size()-1].total_pdf_A, eye_vs, light_vs, (const int)eye_vs.size(), 0);
		const Color result = mis_weight * pt_result.value;
		bpt_result.data.push_back(BidirectionalPathtracingResult::Data(imagebuffer_x, imagebuffer_y, result));
	}
	// num_eye_vertex == 0のとき、光源側からのパスがレンズにヒットしていれば、それをサンプルとして使用する
	if (lt_result.is_lens_hit) {
		const double mis_weight = calc_mis_weight(camera, light_vs[light_vs.size()-1].total_pdf_A, eye_vs, light_vs, 0, (const int)light_vs.size());
		const int lx = lt_result.imagebuffer_x, ly = lt_result.imagebuffer_y;
		const Color result =  mis_weight * lt_result.value;
		bpt_result.data.push_back(BidirectionalPathtracingResult::Data(lx, ly, result));
	}
	for (int num_eye_vertex = 1; num_eye_vertex <= eye_vs.size(); ++num_eye_vertex) {
		for (int num_light_vertex = 1; num_light_vertex <= light_vs.size(); ++num_light_vertex) {
			int target_x = imagebuffer_x, target_y = imagebuffer_y;
			const Vertex &eye_end = eye_vs[num_eye_vertex - 1];
			const Vertex &light_end = light_vs[num_light_vertex - 1];

			// トータルの確率密度計算
			const double total_pdf_A = eye_end.total_pdf_A * light_end.total_pdf_A;
			if (total_pdf_A == 0) // 端点がSupecularとかだとここに入る
				continue;
					
			// MCスループットと重み計算
			Color eye_throughput = eye_end.throughput;
			Color light_throughput = light_end.throughput;
			Color eye_weight(1, 1, 1), light_weight(1, 1, 1);

			// 端点間の交差判定（必要に応じてレンズとの交差判定も行う）
			Intersection intersection;
			const Vec light_end_to_eye_end = eye_end.position - light_end.position;
			Ray test_ray = Ray(light_end.position, normalize(light_end_to_eye_end));
			intersect_scene(test_ray, &intersection);
					
			// nun_light_vertex == 1のとき、非完全拡散光源の場合は相手の頂点の位置次第でMCスループットが変化するため別処理を挟む
			// 今回は完全拡散光源なので単純にemissionの値を入れる
			if (num_light_vertex == 1) {
				light_throughput = spheres[light_vs[0].objectID].emission;
			}
			
			// eyeパスの重み調整
			switch (eye_end.type) {
			case Vertex::OBJECT_TYPE_DIFFUSE:
				eye_weight = multiply(eye_weight, spheres[eye_end.objectID].color / kPI);

				// 端点同士が別の物体で遮蔽されるかどうかを判定する
				if ((intersection.hitpoint.position - eye_end.position).length() >= kEPS) {
					continue;
				}
				break;
			case Vertex::OBJECT_TYPE_LENS:
				// ここに入るのはnum_eye_vertex == 1のときだけ
				{
					// Lightパスを直接レンズにつなげるパターン		
					Vec position_on_lens, position_on_objectplane, position_on_imagesensor, uv_on_imagesensor;
					const double lens_t = camera.intersect_lens(test_ray, &position_on_lens, &position_on_objectplane, &position_on_imagesensor, &uv_on_imagesensor);		
					if (kEPS < lens_t && lens_t < intersection.hitpoint.distance) {
						// レイがレンズにヒット＆イメージセンサにヒット
						const Vec x0_xI = position_on_imagesensor - position_on_lens;
						const Vec x0_xV = position_on_objectplane - position_on_lens;
						const Vec x0_x1 = test_ray.org - position_on_lens;
							
						target_x = (int)uv_on_imagesensor.x;
						target_y = (int)uv_on_imagesensor.y;
						clamp(target_x, 0, camera.image_width_px);
						clamp(target_y, 0, camera.image_height_px);

						eye_weight = camera.W_dash(x0_xV, x0_xI, x0_x1) * eye_weight;
					} else {
						// Lightパスを直接レンズにつなげようとしたが、遮蔽されたりイメージセンサにヒットしなかった場合、終わり
						continue;
					}
				}
				break;
				// eyeパス側の端点がスペキュラや光源（反射率ゼロ）だった場合は重みがゼロになりパス全体の寄与もゼロになる
			case Vertex::OBJECT_TYPE_LIGHT:
			case Vertex::OBJECT_TYPE_SPECULAR:
				eye_weight = Color(0, 0, 0);
				break;
			}

			// lightパスの重み調整
			switch (light_end.type) {
			case Vertex::OBJECT_TYPE_DIFFUSE:
				light_weight = multiply(light_weight, spheres[light_end.objectID].color / kPI);
				break;
			case Vertex::OBJECT_TYPE_LIGHT:
				break;
				// lightパス側の端点がスペキュラやレンズだった場合は重みがゼロになりパス全体の寄与もゼロになる
			case Vertex::OBJECT_TYPE_LENS:
			case Vertex::OBJECT_TYPE_SPECULAR:
				light_weight = Color(0, 0, 0);
				break;
			}

			// 端点間のジオメトリターム
			const double G = (
				std::max(dot(normalize(-1.0*light_end_to_eye_end), eye_end.orienting_normal), 0.0) *
				std::max(dot(normalize(light_end_to_eye_end), light_end.orienting_normal), 0.0)) /
				light_end_to_eye_end.length_squared();

			// MIS重み
			const double mis_weight = calc_mis_weight(camera, total_pdf_A, eye_vs, light_vs, num_eye_vertex, num_light_vertex);

			if (_isnan(mis_weight)) {
				continue;
			}

			const Color result = mis_weight * multiply(multiply(G * multiply(eye_weight, light_weight), eye_throughput), light_throughput) / total_pdf_A;
			bpt_result.data.push_back(BidirectionalPathtracingResult::Data(target_x, target_y, result));
		}
	}

	return bpt_result;
}

};

#endif