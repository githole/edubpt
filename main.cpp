#include "render.h"

int main() {
	edubpt::Camera camera(
		640, 480, 
		edubpt::Vec(50.0, 40.8, 220.0), normalize(edubpt::Vec(0.0, 0.0, -1.0)), edubpt::Vec(0.0, 1.0, 0.0),
		30.0, 30.0, 140.0, 5.0, 
		28.0, 
		256); // 256 samples/pixel

	const int num_threads = 10;

//	render_by_refernce_pathtracing(camera, num_threads);
//	render_by_pathtracing(camera, num_threads);
//	render_by_lighttracing(camera, num_threads);
	edubpt::render_by_bidirectional_pathtracing(camera, num_threads);
	return 0;
}