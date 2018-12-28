#pragma once

#include "glm/glm.h"
#include <vector>

using std::vector;

struct Point3D
{
	float position[3];
	float normal[3];
	float texcoord[2];
	float color[3];
	float velocity[3];
	float force[3];
};

struct Spring
{
	int i, j;					// points indexes
	float length;				// rest length
	float	norm[3];            	// normal vector   
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

	MyPhysicsEngine() {};
	~MyPhysicsEngine() {};

	size_t getAllPointsSize() {
		return sizeof(Point3D) * allPoints.size();
	}

	double distance(const Point3D& a, const Point3D& b) {
		double dx = a.position[0] - b.position[0];
		double dy = a.position[1] - b.position[1];
		double dz = a.position[2] - b.position[2];
		return sqrt(dx*dx + dy*dy + dz*dz);
	}

	void addSoftBall(GLMmodel *model) {
		GLMgroup *group = model->groups;
		int numVertices = 0;
		while (group) {
			numVertices += group->numtriangles * 3;
			group = group->next;
		}

		allPoints.resize(numVertices);
		allSprings.resize(numVertices);
		group = model->groups;
		while (group) {
			for (int i = 0; i < group->numtriangles; ++i) {
				GLMtriangle* triangle = &(model->triangles[group->triangles[i]]);
				for (int j = 0; j < 3; ++j) {
					for (int k = 0; k < 3; ++k) {
						allPoints[i * 3 + j].position[k] = model->vertices[3 * triangle->vindices[j] + k];
						allPoints[i * 3 + j].normal[k] = model->normals[3 * triangle->nindices[j] + k];
						allPoints[i * 3 + j].color[k] = model->materials[triangle->material].diffuse[k];
						allPoints[i * 3 + j].velocity[k] = 0.0f;
						allPoints[i * 3 + j].force[k] = 0.0f;
					}
					for (int k = 0; k < 2; ++k) {
						allPoints[i * 3 + j].texcoord[k] = model->texcoords[2 * triangle->tindices[j] + k];
					}

					allSprings[i * 3 + j].i = i * 3 + j;
					allSprings[i * 3 + j].j = i * 3 + (j+1)%3;
					allSprings[i * 3 + j].length = 
						distance(allPoints[allSprings[i * 3 + j].i], allPoints[allSprings[i * 3 + j].j]);
				}
			}

			group = group->next;
		}
	}

	void addFloor() {

	}

	void update(float deltaTime) {
		//Compute gravity
		for (int i = 0; i < allPoints.size(); ++i) {
			allPoints[i].force[0] = 0.0f;
			allPoints[i].force[1] = MASS * GY * (Pressure - FINAL_PRESSURE >= 0);
			allPoints[i].force[2] = 0.0f;
		}

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

				s.norm[0] = p1.normal[0];
				s.norm[1] = p1.normal[1];
				s.norm[2] = p1.normal[2];
			}
		}

		//Calculate volume
		float volume = 5.0f;
		Point3D center;
		center.position[0] = center.position[1] = center.position[2] = 0;
		for (int i = 0; i < allPoints.size(); ++i) {
			center.position[0] += allPoints[i].position[0];
			center.position[1] += allPoints[i].position[1];
			center.position[2] += allPoints[i].position[2];
		}
		center.position[0] /= allPoints.size();
		center.position[1] /= allPoints.size();
		center.position[2] /= allPoints.size();

		double maxRadius = 0.0;
		for (int i = 0; i < allPoints.size(); ++i) {
			double r = distance(center, allPoints[i]);
			if (r > maxRadius){
				maxRadius = r;
			}
		}
		volume = 4 / 3 * 3.14 * pow(maxRadius, 3);

		//Calculate pressure force
		for (int i = 0; i < allSprings.size(); ++i) {
			Spring& s = allSprings[i];
			Point3D& p1 = allPoints[s.i];
			Point3D& p2 = allPoints[s.j];
			float r12d = distance(p1, p2);

			float pressurev = r12d * Pressure * (1.0f / volume);
			for (int j = 0; j < 3; ++j) {
				p1.force[j] += s.norm[j] * pressurev;
				p2.force[j] += s.norm[j] * pressurev;
			}
		}

		//update the position in x,y,z
		//pos = old_pos + delta_pos 
		//	= old_pos + delta_time*vel
		//	= old_pos + delta_time*(old_vel + delta_time*acc)
		//	= old_pos + delta_time*(old_vel + delta_time*F/m)
		for (int i = 0; i < allPoints.size(); ++i) {
			Point3D& p = allPoints[i];
			for (int j = 0; j < 3; ++j) {
				p.velocity[j] = p.velocity[j] + (p.force[j] / MASS) * deltaTime;
				float delta_pos = p.velocity[j] * deltaTime;
				p.position[j] += delta_pos;
			}
		}

		//update normals for all points
		//the 3 points in the same triangle share one normal
		for (int i = 0; i < allPoints.size(); i += 3) {

		}

		Pressure += FINAL_PRESSURE / 100.0f;
	}
};

const float MyPhysicsEngine::MASS = 1.0f;
const float MyPhysicsEngine::GY = -10.0f;
const float MyPhysicsEngine::FINAL_PRESSURE = 85.0f;
const float MyPhysicsEngine::KS = 5.0f;
const float MyPhysicsEngine::KD = 35.0f;
