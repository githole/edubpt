#include "render.h"

int main() {
	edubpt::Camera camera(
		// 出力画像の解像度（イメージセンサにマッピングされる）
		640, 480, 
		// イメージセンサの中央部の位置
		edubpt::Vec(50.0, 40.8, 220.0),
		// イメージセンサの向いている方向
		normalize(edubpt::Vec(0.0, 0.0, -1.0)),
		// イメージセンサに対する上方ベクトル
		edubpt::Vec(0.0, 1.0, 0.0),
		// イメージセンサの縦幅（物理的大きさ）
		30.0, 
		// イメージセンサとレンズの距離
		30.0, 
		// 合焦距離（レンズとオブジェクトプレーンの距離）
		140.0, 
		// レンズ半径
		5.0, 
		// センサ感度（最終的な出力画素値に対するスケール値。今回はRGB全て同じとする）
		// 適当に決めた
		28.0, 
		// ピクセルあたりのサンプリング回数
		// 双方向パストレの場合、ピクセルあたりのeyeパス・lightパスの追跡回数となる
		10);

	const int num_threads = 10;

	// リファレンス用のパストレーシングで画像生成（双方向パストレとは独立）
//	edubpt::render_by_refernce_pathtracing(camera, num_threads);

	// カメラ側から頂点生成して画像生成
//	edubpt::render_by_pathtracing(camera, num_threads);

	// 光源側から頂点生成して画像生成
//	edubpt::render_by_lighttracing(camera, num_threads);

	// 双方向パストレで画像生成
	edubpt::render_by_bidirectional_pathtracing(camera, num_threads);
	return 0;
}