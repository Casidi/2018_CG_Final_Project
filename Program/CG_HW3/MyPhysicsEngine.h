#pragma once

#include "glm/glm.h"
#include <vector>

using std::vector;

struct Point3D
{
	float *position;
	float normal[3];
	float velocity[3];
	float force[3];
	int normalCount;
};

struct Spring
{
	int i, j;					// points indexes
	float length;				// rest length
};

class MyPhysicsEngine
{
public:
	float Pressure = 0;
	static const float MASS;
	static const float GY;
	static const float FINAL_PRESSURE;
	static const float KS, KD;

	vector<Spring> allSprings;
	vector<Point3D> allPoints;
	vector<GLMtriangle*> allTriangles;

	MyPhysicsEngine() {};
	~MyPhysicsEngine() {};

	static double length(const float a[3]) {
		return sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
	}

	double distance(const Point3D& a, const Point3D& b) {
		double dx = a.position[0] - b.position[0];
		double dy = a.position[1] - b.position[1];
		double dz = a.position[2] - b.position[2];
		return sqrt(dx*dx + dy*dy + dz*dz);
	}

	float distance(const float* a, const float* b) {
		float dx = a[0] - b[0];
		float dy = a[1] - b[1];
		float dz = a[2] - b[2];
		return sqrt(dx*dx + dy * dy + dz * dz);
	}

	static void cross(const float u[3], const float v[3], float result[3]) {
		result[0] = u[1] * v[2] - u[2] * v[1];
		result[1] = u[2] * v[0] - u[0] * v[2];
		result[2] = u[0] * v[1] - u[1] * v[0];

		double len = length(result);
		result[0] /= len;
		result[1] /= len;
		result[2] /= len;
	}

	void normalize(float* a) {
		double len = length(a);
		a[0] /= len;
		a[1] /= len;
		a[2] /= len;
	}

	void addSpring(GLMmodel *model, int i, int j) {
		Spring newSpring;
		newSpring.i = i;
		newSpring.j = j;
		newSpring.length =
			distance(&(model->vertices[newSpring.i * 3]), &(model->vertices[newSpring.j * 3]));
		allSprings.push_back(newSpring);
	}

	static float triangleArea(float* a, float* b, float* c) {
		float ab[3], ac[3], crossVal[3];
		for (int i = 0; i < 3; ++i) {
			ab[i] = b[i] - a[i];
			ac[i] = c[i] - a[i];
		}
		cross(ab, ac, crossVal);
		return 0.5 * length(crossVal);
	}

	void addSoftBall(GLMmodel *model) {
		GLMgroup *group = model->groups;
		int numVertices = 0;
		while (group) {
			numVertices += group->numtriangles * 3;
			group = group->next;
		}

		allPoints.resize(model->numvertices+1);
		for (int i = 0; i < allPoints.size(); ++i) {
			allPoints[i].position = &(model->vertices[3 * i]);
			allPoints[i].velocity[0] = allPoints[i].velocity[1] = allPoints[i].velocity[2] = 0;
			allPoints[i].normalCount = 0;
		}

		group = model->groups;
		while (group) {
			for (int i = 0; i < group->numtriangles; ++i) {
				GLMtriangle* triangle = &(model->triangles[group->triangles[i]]);
				allTriangles.push_back(triangle);

				addSpring(model, triangle->vindices[0] , triangle->vindices[1]);
				addSpring(model, triangle->vindices[1] , triangle->vindices[2]);
				addSpring(model, triangle->vindices[2] , triangle->vindices[0]);

				for (int j = 0; j < 3; ++j) {
					float* pointNorm = allPoints[triangle->vindices[j]].normal;
					for (int k = 0; k < 3; ++k) {
						pointNorm[k] += model->normals[triangle->nindices[j]*3 + k];
						allPoints[triangle->vindices[j]].normalCount++;
					}
				}
			}

			group = group->next;
		}

		for (int i = 1; i < allPoints.size(); ++i) {
			normalize(allPoints[i].normal);
		}
	}

	void update(float deltaTime) {
		deltaTime *= 4;
		//Compute gravity
		for (int i = 0; i < allPoints.size(); ++i) {
			allPoints[i].force[0] = 0.0f;
			allPoints[i].force[1] = MASS * GY * (Pressure - FINAL_PRESSURE >= 0);
			allPoints[i].force[2] = 0.0f;
		}

		//Spring force
		for (int i = 0; i < allSprings.size(); ++i) {
			Spring& s = allSprings[i];
			Point3D& p1 = allPoints[s.i];
			Point3D& p2 = allPoints[s.j];
			float r12d = distance(p1, p2);

			if (r12d != 0) {
				float v12[3];
				for (int j = 0; j < 3; ++j)
					v12[j] = p1.velocity[j] - p2.velocity[j];

				float f = (r12d - s.length)*KS 
					+ (v12[0]*(p1.position[0] - p2.position[0])
						+ v12[1] * (p1.position[1] - p2.position[1])
						+ v12[2] * (p1.position[2] - p2.position[2])) * KD / r12d;

				float F;
				for (int j = 0; j < 3; ++j) {
					F = ((p1.position[j] - p2.position[j]) / r12d) * f;
					p1.force[j] -= F;
					p2.force[j] += F;
				}
			}
		}

		//Calculate volume
		float max3[3] = { allPoints[1].position[0], allPoints[1].position[1], allPoints[1].position[2] };
		float min3[3] = { allPoints[1].position[0], allPoints[1].position[1], allPoints[1].position[2] };
		for (int i = 1; i < allPoints.size(); ++i) {
			for (int j = 0; j < 3; ++j) {
				if (max3[j] < allPoints[i].position[j])
					max3[j] = allPoints[i].position[j];
				if (min3[j] > allPoints[i].position[j])
					min3[j] = allPoints[i].position[j];
			}
		}
		float volume = (max3[0] - min3[0]) * (max3[1] - min3[1]) * (max3[2] - min3[2]);

		//Calculate pressure force
		for (int i = 0; i < allTriangles.size(); ++i) {
			Point3D& p1 = allPoints[allTriangles[i]->vindices[0]];
			Point3D& p2 = allPoints[allTriangles[i]->vindices[1]];
			Point3D& p3 = allPoints[allTriangles[i]->vindices[2]];
			float area = triangleArea(p1.position, p2.position, p3.position);

			float pressurev = area * Pressure * (1.0f / volume);
			for (int j = 0; j < 3; ++j) {
				p1.force[j] += p1.normal[j] * pressurev;
				p2.force[j] += p2.normal[j] * pressurev;
				p3.force[j] += p3.normal[j] * pressurev;
			}
		}

		//update the position
		for (int i = 1; i < allPoints.size(); ++i) {
			Point3D& p = allPoints[i];

			p.force[0] /= p.normalCount;
			p.force[1] /= p.normalCount;
			p.force[2] /= p.normalCount;

			//x
			p.velocity[0] = p.velocity[0] + (p.force[0] / MASS) * deltaTime;
			float delta_pos = p.velocity[0] * deltaTime;
			p.position[0] += delta_pos;

			p.velocity[1] = p.velocity[1] + (p.force[1] / MASS) * deltaTime;
			delta_pos = p.velocity[1] * deltaTime;
			float Y_BOUND = -3.5f;
			if (p.position[1] + delta_pos < Y_BOUND) {
				p.position[1] = Y_BOUND;
				p.velocity[1] = -0.3 * p.velocity[1];
				p.velocity[0] *= 0.95;
				p.velocity[2] *= 0.95;
			}
			else {
				p.position[1] += delta_pos;
			}

			//z
			p.velocity[2] = p.velocity[2] + (p.force[2] / MASS) * deltaTime;
			delta_pos = p.velocity[2] * deltaTime;
			p.position[2] += delta_pos;
		}

		if(Pressure < FINAL_PRESSURE)
			Pressure += FINAL_PRESSURE / 10.0f;
	}
};

const float MyPhysicsEngine::MASS = 2.0f;
const float MyPhysicsEngine::GY = -3.0f;
const float MyPhysicsEngine::FINAL_PRESSURE = 100.0f;
const float MyPhysicsEngine::KS = 2000.0f;
const float MyPhysicsEngine::KD = 20.0f;
