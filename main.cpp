
#include "hdr.h"
#include "random.h"
#include "scene.h"
#include "util.h"
#include "camera.h"
#include "vertex.h"
#include "specular.h"
#include "sampling.h"
#include "camera.h"

#include "reference_pathtracing.h"
#include "pathtracing.h"
#include "lighttracing.h"
#include "bidirectional_pathtracing.h"

#include <omp.h>
#include <assert.h>

using namespace edubpt;


int render_by_refernce_pathtracing(const Camera &camera) {
	const int width = camera.image_width_px;
	const int height = camera.image_height_px;
	const int spp = camera.samples_per_pixel;

	HDRImage hdr(width, height);
	std::vector<Color> image_buffer(width * height);
	
	for (int sample = 0; sample < spp; ++sample) {
		std::cout << sample << " " << std::endl;
#pragma omp parallel for schedule(dynamic, 1) num_threads(10)
		for (int y = 0; y < height; ++y) {
			Random random(y * spp + sample);
			for (int x = 0; x < width; ++x) {
				Vec position_on_imagesensor, position_on_objectplane, position_on_lens;
				double P_Image, P_lens;
				camera.sample_points(x, y, &random, &position_on_imagesensor, &position_on_objectplane, &position_on_lens, &P_Image, &P_lens);

				// 方向
				const Vec dir = normalize(position_on_objectplane - position_on_lens);
				// L(x0 <- x1) (= L(xI <- x0) = L(x0 <- xV))
				const Color L = radiance_no_recursion(Ray(position_on_lens, dir), &random);
				if (!is_valid_value(L)) // 異常値は除外する
					continue;

				const Vec x0_xI = position_on_imagesensor - position_on_lens;
				const double coefficient = pow(dot(normalize(-1.0 * camera.imagesensor_dir), normalize(x0_xI)), 2.0) / x0_xI.length_squared(); // センサ上の点の半球積分を、レンズ上の点の積分に変数変換した時に導入される係数

				image_buffer[y * width + x] = image_buffer[y * width + x] +
					camera.W * L * coefficient / (P_Image * P_lens);
			}
		}

		if (sample % 128 == 0) {
			char buf[256];
			sprintf(buf, "ref_pt_%03d.hdr", sample);
			for (int y = 0; y < height; ++y)
				for (int x = 0; x < width; ++x)
					*hdr.image_ptr(width - x - 1, height - y - 1) = image_buffer[y * width + x] / ((double)width * height * (sample + 1));
			hdr.save(buf);
		}
	}
	return 0;
}

int render_by_pathtracing(const Camera &camera) {
	const int width = camera.image_width_px;
	const int height = camera.image_height_px;
	const int spp = camera.samples_per_pixel;

	HDRImage hdr(width, height);
	std::vector<Color> image_buffer(width * height);
	
	for (int sample = 0; sample < spp; ++sample) {
		std::cout << sample << " " << std::endl;
#pragma omp parallel for schedule(dynamic, 1) num_threads(10)
		for (int y = 0; y < height; ++y) {
			Random random(y * spp + sample);
			for (int x = 0; x < width; ++x) {
				std::vector<Vertex> vertices;
				const PathtracingResult result = generate_vertices_by_pathtracing(camera, x, y, &random, vertices);

				if (result.is_light_hit) {
					if (is_valid_value(result.value)) {
						image_buffer[y * width + x] = image_buffer[y * width + x] + result.value;
					}
				}
			}
		}

		if (sample % 128 == 0) {
			char buf[256];
			sprintf(buf, "pt_%03d.hdr", sample);
			for (int y = 0; y < height; ++y)
				for (int x = 0; x < width; ++x) {
					// std::cout << image_buffer[y * width + x].x << " ";
					*hdr.image_ptr(width - x - 1, height - y - 1) = image_buffer[y * width + x] / ((double)width * height * (sample + 1));
				}
			hdr.save(buf);
		}
	}
	return 0;
}

void render_by_lighttracing(const Camera &camera) {
	const int width = camera.image_width_px;
	const int height = camera.image_height_px;

	HDRImage hdr(width, height);
	const int num_threads = 10;
	std::vector<Color> image_buffer(width * height * num_threads);
	omp_set_num_threads(num_threads);

	int thread_id;
#pragma omp parallel private(thread_id)
	{
		thread_id = omp_get_thread_num();
		Random random(thread_id * 32);

		for (long long sample = 0; sample < (long long)camera.samples_per_pixel * width * height; ++sample) {
			std::vector<Vertex> vertices;
			LighttracingResult result = generate_vertices_by_lighttracing(camera, &random, vertices);
			if (result.is_lens_hit) {
				if (is_valid_value(result.value)) {
					const int x = result.imagebuffer_x, y = result.imagebuffer_y;
					image_buffer[(thread_id * width * height) + y * width + x] = image_buffer[(thread_id * width * height) + y * width + x] + result.value;
				}
			}

			if (sample % (16 * width * height) == 0) {
#pragma omp barrier
				if (thread_id == 0) {
					std::cout << sample << " ";
					char buf[256];
					sprintf(buf, "lt_%03d.hdr", sample / (width * height));
					for (int y = 0; y < height; ++y)
						for (int x = 0; x < width; ++x)
							*hdr.image_ptr(x, y) = Color();
					for (int t = 0; t < num_threads; ++t)
						for (int y = 0; y < height; ++y)
							for (int x = 0; x < width; ++x)
								*hdr.image_ptr(width - x - 1, height - y - 1) = *hdr.image_ptr(width - x - 1, height - y - 1) + image_buffer[(t * width * height) + y * width + x] / (sample + 1) / num_threads;
					hdr.save(buf);
				}
#pragma omp barrier
			}
		}
	}
}


void render_by_bidirectional_pathtracing(const Camera &camera) {
	const int width = camera.image_width_px;
	const int height = camera.image_height_px;

	HDRImage hdr(width, height);
	const int num_threads = 10;
	std::vector<Color> image_buffer(width * height * num_threads);
	omp_set_num_threads(num_threads);

	int thread_id;
#pragma omp parallel private(thread_id)
	{
		thread_id = omp_get_thread_num();
		Random random(thread_id * 32);
		
		for (int sample = 0; sample < camera.samples_per_pixel; ++sample) {
			std::cout << sample << std::endl;
			for (int y = 0; y < height; ++y) {
				for (int x = 0; x < width; ++x) {
					BidirectionalPathtracingResult bpt_result = bidirectional_pathtracing(camera, x, y, &random);

					for (int i = 0; i < bpt_result.data.size(); ++i) {
						const int ix = bpt_result.data[i].imagebuffer_x;
						const int iy = bpt_result.data[i].imagebuffer_y;
						const int idx = (thread_id * width * height) + (iy * width + ix);
						if (is_valid_value(bpt_result.data[i].value))
							image_buffer[idx] = image_buffer[idx] + bpt_result.data[i].value;
					}
				}
			}
			if (1) {
#pragma omp barrier
				if (thread_id == 0) {
					std::cout << sample << " ";
					char buf[256];
					sprintf(buf, "bpt_%03d.hdr", sample);
					for (int y = 0; y < height; ++y)
						for (int x = 0; x < width; ++x)
							*hdr.image_ptr(x, y) = Color();
					for (int t = 0; t < num_threads; ++t)
						for (int y = 0; y < height; ++y)
							for (int x = 0; x < width; ++x) {
								const int idx = (t * width * height) + (y * width + x);
								*hdr.image_ptr(width - x - 1, height - y - 1) = 
									*hdr.image_ptr(width - x - 1, height - y - 1) + image_buffer[idx] / ((double) width * height * (sample + 1)) / num_threads;
							}
					hdr.save(buf);
				}
#pragma omp barrier
			}
		}
	}
}



int main() {
	Camera camera(
		640, 480, 
		Vec(50.0, 40.8, 220.0), normalize(Vec(0.0, 0.0, -1.0)), Vec(0.0, 1.0, 0.0),
		30.0, 30.0, 140.0, 5.0, 
		28.0, 8192);

//	render_by_refernce_pathtracing(camera);
//	render_by_pathtracing(camera);
//	render_by_lighttracing(camera);

	render_by_bidirectional_pathtracing(camera);
	return 0;
}