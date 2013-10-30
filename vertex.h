#ifndef _VERTEX_H_
#define _VERTEX_H_

#include "vec.h"
#include "material.h"

namespace edubpt {

struct Vertex {
	double total_pdf_A; // この頂点を含めた、この頂点にいたるまでのサンプリング確率密度
	Color throughput; // この頂点に至るまでの、モンテカルロスループット（最終的なモンテカルロ積分の分母に相当したり、レンズ側からならBRDFやGの積になる）
	Vec position;
	int objectID;
	Vec orienting_normal;
	Vec object_normal; // この頂点の位置の物体の元の法線

	enum ObjectType {
		OBJECT_TYPE_LIGHT,
		OBJECT_TYPE_LENS,
		OBJECT_TYPE_SPECULAR,
		OBJECT_TYPE_DIFFUSE,
	};
	ObjectType type;

	Vertex(const Vec &position, const Vec &orienting_normal, const Vec &object_normal, const int objectID, const ObjectType type, const double total_pdf_A, const Color &throughput) :
	position(position), orienting_normal(orienting_normal), object_normal(object_normal), objectID(objectID), total_pdf_A(total_pdf_A), throughput(throughput), type(type) {}
};

};

#endif