#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <algorithm>

#include <glm/gtc/matrix_transform.hpp>

#include "bumper_graph.h"

#include <unordered_set>

using namespace SM;
using namespace SM::Graph;

// ReSharper disable once CppDFATimeOver
float BumperGraph::PUSH_EPSILON = 0.000001f;

void BumperCapsuloid::translate(const glm::vec3 &t)
{
	for (auto& p : neib)
		p.translate(t);

	for (auto& c : cones)
		c.translate(t);
}

void BumperCapsuloid::scale(const float r)
{
	for (auto& p : neib)
		p.scale(r);

	for (auto& c : cones)
		c.scale(r);
}

void BumperCapsuloid::rotate(const glm::mat3 &rot)
{
	for (auto& p : neib)
		p.rotate(rot);

	for (auto& c : cones)
		c.rotate(rot);
}

void BumperCapsuloid::rotateAround(const glm::vec3& axis, const float angle)
{
	const auto rot = glm::mat3(glm::rotate(glm::mat4(1.0f), angle, axis));
	rotate(rot);
}

bool BumperCapsuloid::operator == (const BumperCapsuloid &bc) const
{
	return sphereIndex[0] == bc.sphereIndex[0] && sphereIndex[1] == bc.sphereIndex[1];
}

BumperCapsuloid::BumperCapsuloid (const int indexA, const int indexB, const Sphere &a, const Sphere &b)
{
	sphereIndex[0] = indexA < indexB ? indexA : indexB;
	sphereIndex[1] = indexA < indexB ? indexB : indexA;

	const Sphere &sA = indexA < indexB ? a : b;
	const Sphere &sB = indexA < indexB ? b : a;

	const glm::vec3 d = sB.center - sA.center;
	const float dLengthSquared = dot(d, d);
	dLength = std::sqrt(dLengthSquared);
	dNorm = d / dLength;
	dR = (sB.radius - sA.radius) / dLength;

	const float r0r1 = sB.radius - sA.radius;
	const float dRSquared = r0r1 * r0r1;
	const float diff = dLengthSquared - dRSquared;
	if (diff <= 0.0f || dLength <= 0.0f)
	{
		std::cerr << "Capsuloid degenerate" << std::endl;
		k = 0.0f;
	}
	else
		k = r0r1 / (std::sqrt(diff) * dLength);

	const float diffRadius = fabs(sB.radius - sA.radius);
	float alpha = 0.0f;
	if (dLength > 1e-6f && diffRadius < dLength) alpha = asinf(diffRadius / dLength);

	cones[0].apex = sA.center;
	cones[0].axis = dNorm;
	cones[0].halfAngle = alpha;

	cones[1].apex = sB.center;
	cones[1].axis = -dNorm;
	cones[1].halfAngle = alpha;
}

void BumperSphere::translate(const glm::vec3 &t)
{
	for (auto& cone : cones)
		cone.translate(t);
}

void BumperSphere::scale(float r)
{
	for (auto& cone : cones)
		cone.scale(r);
}

void BumperSphere::rotate(const glm::mat3 &rot)
{
	for (auto& cone : cones)
		cone.rotate(rot);
}

void BumperSphere::rotateAround(const glm::vec3& axis, const float angle)
{
	const auto rot = glm::mat3(glm::rotate(glm::mat4(1.0f), angle, axis));
	rotate(rot);
}

bool BumperPrysmoid::operator==(const BumperPrysmoid &bp) const
{
	return sphereIndex[0] == bp.sphereIndex[0] && sphereIndex[1] == bp.sphereIndex[1] && sphereIndex[2] == bp.sphereIndex[2];
}

void BumperPrysmoid::translate(const glm::vec3 &t)
{
	upperPlane.translate(t);
	midPlane.translate(t);

	m0.translate(t);
	m1.translate(t);
	m2.translate(t);
}

void BumperPrysmoid::scale(const float r)
{
	upperPlane.scale(r);
	midPlane.scale(r);

	m0.scale(r);
	m1.scale(r);
	m2.scale(r);
}

void BumperPrysmoid::rotate(const glm::mat3 &rot)
{
	upperPlane.rotate(rot);
	midPlane.rotate(rot);

	m1.rotate(rot);
	m2.rotate(rot);
	m0.rotate(rot);
}

void BumperPrysmoid::rotateAround(const glm::vec3& axis, const float angle)
{
	const auto rotMatrix = glm::mat3(glm::rotate(glm::mat4(1.0f), angle, axis));
	rotate(rotMatrix);
}

BumperPrysmoid::BumperPrysmoid (
		int indexA,
		int indexB,
		int indexC,
		const Sphere &v0,
		const Sphere &v1,
		const Sphere &v2,
		int direction)
{
	std::vector orderedIndices = {indexA, indexB, indexC};
	std::ranges::sort(orderedIndices);

	sphereIndex[0] = orderedIndices[0];
	sphereIndex[1] = orderedIndices[1];
	sphereIndex[2] = orderedIndices[2];

	upperPlane = computeUpperPlane(v0, v1, v2, direction);

	glm::vec3 u = v1.center - v0.center;
	glm::vec3 v = v2.center - v0.center;

	midPlane = Plane(cross(u, v), v0.center);
	if (direction < 0) midPlane.flip();

	m0 = Plane(cross(v1.center - v2.center, upperPlane.n) * static_cast<float>(direction), v1.center);
	m1 = Plane(cross(v2.center - v0.center, upperPlane.n) * static_cast<float>(direction), v2.center);
	m2 = Plane(cross(v0.center - v1.center, upperPlane.n) * static_cast<float>(direction), v0.center);
}

#ifdef DEBUG
inline void savePrysmoidAsSM(const std::string& filename, const Sphere &sa, const Sphere &sb, const Sphere &sc)
{
	std::ofstream ofs(filename);

	ofs << "Sphere Mesh 2.0" << std::endl;
	ofs << 3;
	ofs << " " << 1;
	ofs << " " << 0 << std::endl;

	ofs << sa.center.x << " " << sa.center.y << " " << sa.center.z << " " << sa.radius << std::endl;
	ofs << sb.center.x << " " << sb.center.y << " " << sb.center.z << " " << sb.radius << std::endl;
	ofs << sc.center.x << " " << sc.center.y << " " << sc.center.z << " " << sc.radius << std::endl;

	ofs << 0 << " " << 1 << " " << 2 << std::endl;
	ofs.close();
}
#endif

Plane BumperPrysmoid::computeUpperPlane(const Sphere &sa, const Sphere &sb, const Sphere &sc, int direction)
{
	glm::vec3 a = sa.center;
	glm::vec3 b = sb.center;
	glm::vec3 c = sc.center;

	auto sign = static_cast<float>(direction);

	glm::vec3 n = sign * cross(b - a, c - a);
	n = normalize(n);
	glm::vec3 startN = n;

	for (int i = 0; i < 1000; i++){
		a = sa.center + n * sa.radius;
		b = sb.center + n * sb.radius;
		c = sc.center + n * sc.radius;

		glm::vec3 new_n = normalize(sign * cross(b - a, c - a));
		if (dot(n, new_n) >= 0.999f)
		{
			if (dot(startN, new_n) < 0)
				std::cerr << "Degenerate Prysmoid detected (Normal flipped)" << std::endl;

			return {new_n, a};
		}
		n = new_n;
	}

	std::cerr << "Degenerate Prysmoid detected (Normal diverged)" << std::endl;

	return {n, a};
}

std::pair<int, float> BumperGraph::signedDistanceFromBumper(const int i, const glm::vec3& p) const
{
	if (bumper[i].shapeType == Bumper::SPHERE)
		return sampleSignedDistanceFromSphere(i, p);
	if (bumper[i].shapeType == Bumper::CAPSULOID)
		return sampleSignedDistanceFromCapsuloid(i, p);
	if (bumper[i].shapeType == Bumper::PRYSMOID)
		return sampleSignedDistanceFromPrysmoid(i, p);

	return sampleSignedDistanceFromQuad(i, p);
}

void BumperGraph::sampleDistanceFromBumperWithPosition(const glm::vec<3, float> &p, const int proposedIndex, glm::vec3 &closestPos) const
{
	glm::vec3 closestNorm;
	if (bumper[proposedIndex].shapeType == Bumper::SPHERE)
		signedDistanceFromSphere(proposedIndex, p, closestPos, closestNorm);
	else if (bumper[proposedIndex].shapeType == Bumper::CAPSULOID)
		signedDistanceFromCapsuloid(proposedIndex, p, closestPos, closestNorm);
	else
		signedDistanceFromPrysmoid(proposedIndex, p, closestPos, closestNorm);
}

void BumperGraph::translate(const glm::vec3 &t)
{
	for (Sphere &s : sphere)
		s.center += t;

	for (auto& bn : bumper)
		if (bn.shapeType == Bumper::PRYSMOID)
		{
			auto &bp = std::get<BumperPrysmoid>(bn.bumper);
			bp.translate(t);
		}
}

void BumperGraph::scale(const float s)
{
	for (Sphere& sp : sphere)
	{
		sp.center *= s;
		sp.radius *= s;
	}

	for (auto& bn : bumper)
		if (bn.shapeType == Bumper::CAPSULOID)
		{
			auto &bc = std::get<BumperCapsuloid>(bn.bumper);
			bc.dLength *= s;

			const float BA = bc.dLength;
			const float r0r1 = sphere[bc.sphereIndex[1]].radius - sphere[bc.sphereIndex[0]].radius;

			const float sq_BA = BA * BA;
			const float sq_r0r1 = r0r1 * r0r1;

			const float L = std::sqrt(sq_BA - sq_r0r1);

			bc.k = r0r1 / (L * BA);
		}
		else if (bn.shapeType == Bumper::PRYSMOID)
		{
			auto &bp = std::get<BumperPrysmoid>(bn.bumper);
			bp.scale(s);
		}
}

void BumperGraph::rotateY(const int angle)
{
	if (angle % 90 != 0)
	{
		std::cerr << "Error: Rotation angle is not a multiple of 90 degrees, not applying rotation." << std::endl;
		return;
	}

	const float angleRadians = glm::radians(static_cast<float>(angle));
	const auto rotationMatrix = glm::mat3(rotate(glm::mat4(1.0f), angleRadians, glm::vec3(0, 1, 0)));

	for (Sphere &s : sphere)
		s.center = rotationMatrix * s.center;

	for (auto& bn : bumper)
		if (bn.shapeType == Bumper::CAPSULOID)
		{
			auto &bc = std::get<BumperCapsuloid>(bn.bumper);
			bc.dNorm = rotationMatrix * bc.dNorm;
		}
		else if (bn.shapeType == Bumper::PRYSMOID)
		{
			auto &bp = std::get<BumperPrysmoid>(bn.bumper);
			bp.rotate(rotationMatrix);
		}
}

bool BumperGraph::pushOutside(glm::vec3 &p, glm::vec3 &n, int &bumperIndex) const
{
	if (bumperIndex == -1) return pushOutsideBruteForce(p, n, bumperIndex);

	// TODO: Check the current bumperIndex, and all the neighbors or the current one
	if (bumper[bumperIndex].shapeType == Bumper::SPHERE)
	{
		if (!pushOutsideSphere(p, n, bumperIndex)) return false;

		for (const Cone& c : std::get<BumperSphere>(bumper[bumperIndex].bumper).cones)
		{
			int idx = c.capsuloidIndex;
			if (c.contains(p) && pushOutsideCapsuloid(p, n, idx))
			{
				bumperIndex = idx;
				return true;
			}
		}
	}
	if (bumper[bumperIndex].shapeType == Bumper::CAPSULOID)
	{
		if(!pushOutsideCapsuloid(p, n, bumperIndex)) return false;

		const int idx = bumperIndex;
		if (!pushOutsideCapsuloid(p, n, bumperIndex))
			return false;

		// Skipping opp planes
		for (int index : std::get<BumperCapsuloid>(bumper[idx].bumper).neibBumpers)
		{
			switch (bumper[index].shapeType)
			{
				case Bumper::PRYSMOID:
					if (pushOutsidePrysmoid(p, n, index))
					{
						bumperIndex = index;
						return true;
					}
					break;
				case Bumper::QUAD:
					if(pushOutsideQuad(p, n, index))
					{
						bumperIndex = index;
						return true;
					}
					break;
				default: break;
			}
		}
	}
	if (bumper[bumperIndex].shapeType == Bumper::PRYSMOID)	return pushOutsidePrysmoid(p, n, bumperIndex);
	if (bumper[bumperIndex].shapeType == Bumper::QUAD)		return pushOutsideQuad(p, n, bumperIndex);

	return false;
}

void BumperGraph::initializeBumperSpheres()
{
	for (int i = 0; i < sphere.size(); i++)
	{
		auto bs = BumperSphere();
		bs.sphereIndex = i;

		Bumper node {};
		node.shapeType = Bumper::SPHERE;
		node.bumper = bs;

		bumper.push_back(node);
	}
}

void BumperGraph::initializeBumperCapsuloids(const SphereMesh &sm, std::vector<std::vector<int>> &capsuloidAdj)
{
	for(const auto &[indices] : sm.capsuloids)
	{
		auto bc = BumperCapsuloid(
			indices[0],
			indices[1],
			sphere[indices[0]],
			sphere[indices[1]]);

		if (capsuloidAdj[indices[0]][indices[1]] != -1 || capsuloidAdj[indices[1]][indices[0]] != -1)
			continue;

		bc.neibExtreme[0] = indices[0];
		bc.neibExtreme[1] = indices[1];

		Cone& cone0 = bc.cones[0];
		Cone& cone1 = bc.cones[1];
		cone0.capsuloidIndex = bumper.size();
		cone1.capsuloidIndex = bumper.size();
		std::get<BumperSphere>(bumper[bc.sphereIndex[0]].bumper).cones.push_back(cone0);
		std::get<BumperSphere>(bumper[bc.sphereIndex[1]].bumper).cones.push_back(cone1);
		bc.hasFather = false;

		Bumper node {};
		node.shapeType = Bumper::CAPSULOID;
		node.bumper = bc;

		bumper.push_back(node);

		capsuloidAdj[indices[0]][indices[1]] = bumper.size() - 1;
		capsuloidAdj[indices[1]][indices[0]] = bumper.size() - 1;
	}
}

void BumperGraph::initializeBumperPrysmoids(const SphereMesh &sm, std::vector<std::vector<int>> &capsuloidAdj)
{
	for (const auto &[indices] : sm.prysmoids)
	{
		auto bpUpper = BumperPrysmoid(
			indices[0],
			indices[1],
			indices[2],
			sphere[indices[0]],
			sphere[indices[1]],
			sphere[indices[2]],
			1
		);

		auto bpLower = BumperPrysmoid(
			indices[0],
			indices[1],
			indices[2],
			sphere[indices[0]],
			sphere[indices[1]],
			sphere[indices[2]],
			-1
		);

		auto bc0 = BumperCapsuloid(
			indices[0],
			indices[1],
			sphere[indices[0]],
			sphere[indices[1]]);

		auto bc1 = BumperCapsuloid(
			indices[0],
			indices[2],
			sphere[indices[0]],
			sphere[indices[2]]);

		auto bc2 = BumperCapsuloid(
			indices[1],
			indices[2],
			sphere[indices[1]],
			sphere[indices[2]]);

		Bumper nodeUpper {};
		nodeUpper.shapeType = Bumper::PRYSMOID;
		bpUpper.neibOpp = bumper.size() + 1;
		nodeUpper.bumper = bpUpper;

		Bumper nodeLower {};
		nodeLower.shapeType = Bumper::PRYSMOID;
		bpLower.neibOpp = bumper.size();
		nodeLower.bumper = bpLower;

		Bumper edge0 {};
		edge0.shapeType = Bumper::CAPSULOID;
		bc0.neibExtreme[0] = indices[0];
		bc0.neibExtreme[1] = indices[1];
		edge0.bumper = bc0;

		Bumper edge1 {};
		edge1.shapeType = Bumper::CAPSULOID;
		bc1.neibExtreme[0] = indices[0];
		bc1.neibExtreme[1] = indices[2];
		edge1.bumper = bc1;

		Bumper edge2 {};
		edge2.shapeType = Bumper::CAPSULOID;
		bc2.neibExtreme[0] = indices[1];
		bc2.neibExtreme[1] = indices[2];
		edge2.bumper = bc2;

		bumper.push_back(nodeUpper);
		bumper.push_back(nodeLower);

		int upperIndex = bumper.size() - 2;
		int lowerIndex = bumper.size() - 1;

		if (capsuloidAdj[indices[0]][indices[1]] == -1 || capsuloidAdj[indices[1]][indices[0]] == -1)
		{
			capsuloidAdj[indices[0]][indices[1]] = bumper.size();
			capsuloidAdj[indices[1]][indices[0]] = bumper.size();

			bumper.push_back(edge0);
		}

		int idx = capsuloidAdj[indices[0]][indices[1]];

		std::get<BumperCapsuloid>(bumper[idx].bumper).neibBumpers.push_back(upperIndex);
		std::get<BumperCapsuloid>(bumper[idx].bumper).neibBumpers.push_back(lowerIndex);

		std::get<BumperPrysmoid>(bumper[upperIndex].bumper).neibSide[2] = idx;
		std::get<BumperPrysmoid>(bumper[lowerIndex].bumper).neibSide[2] = idx;

		if (capsuloidAdj[indices[0]][indices[2]] == -1 || capsuloidAdj[indices[2]][indices[0]] == -1)
		{
			capsuloidAdj[indices[0]][indices[2]] = bumper.size();
			capsuloidAdj[indices[2]][indices[0]] = bumper.size();

			bumper.push_back(edge1);
		}

		idx = capsuloidAdj[indices[0]][indices[2]];

		std::get<BumperCapsuloid>(bumper[idx].bumper).neibBumpers.push_back(upperIndex);
		std::get<BumperCapsuloid>(bumper[idx].bumper).neibBumpers.push_back(lowerIndex);

		std::get<BumperPrysmoid>(bumper[upperIndex].bumper).neibSide[1] = idx;
		std::get<BumperPrysmoid>(bumper[lowerIndex].bumper).neibSide[1] = idx;

		if (capsuloidAdj[indices[1]][indices[2]] == -1 || capsuloidAdj[indices[2]][indices[1]] == -1)
		{
			capsuloidAdj[indices[1]][indices[2]] = bumper.size();
			capsuloidAdj[indices[2]][indices[1]] = bumper.size();

			bumper.push_back(edge2);
		}

		idx = capsuloidAdj[indices[1]][indices[2]];

		std::get<BumperCapsuloid>(bumper[idx].bumper).neibBumpers.push_back(upperIndex);
		std::get<BumperCapsuloid>(bumper[idx].bumper).neibBumpers.push_back(lowerIndex);

		std::get<BumperPrysmoid>(bumper[upperIndex].bumper).neibSide[0] = idx;
		std::get<BumperPrysmoid>(bumper[lowerIndex].bumper).neibSide[0] = idx;
	}
}

void BumperGraph::initializeBumperQuads(const SphereMesh &sm, std::vector<std::vector<int>> &capsuloidAdj)
{
	for (const Quadrilateral& p : sm.quadrilaterals)
	{
		FourSpheres upper = computeFourSpheresFrom(p, 1);
		FourSpheres lower = computeFourSpheresFrom(p, -1);

		auto bqUpper = BumperQuad(upper);
		auto bqLower = BumperQuad(lower);

		auto bc0 = BumperCapsuloid(
			p.indices[0],
			p.indices[1],
			sphere[p.indices[0]],
			sphere[p.indices[1]]);

		auto bc1 = BumperCapsuloid(
			p.indices[1],
			p.indices[2],
			sphere[p.indices[1]],
			sphere[p.indices[2]]);

		auto bc2 = BumperCapsuloid(
			p.indices[2],
			p.indices[3],
			sphere[p.indices[2]],
			sphere[p.indices[3]]);

		auto bc3 = BumperCapsuloid(
			p.indices[3],
			p.indices[0],
			sphere[p.indices[3]],
			sphere[p.indices[0]]);

		Bumper nodeUpper {};
		nodeUpper.shapeType = Bumper::QUAD;
		bqUpper.neibOpp = bumper.size() + 1;
		nodeUpper.bumper = bqUpper;

		Bumper nodeLower {};
		nodeLower.shapeType = Bumper::QUAD;
		bqLower.neibOpp = bumper.size();
		nodeLower.bumper = bqLower;

		Bumper edge0 {};
		edge0.shapeType = Bumper::CAPSULOID;
		bc0.neibExtreme[0] = p.indices[0];
		bc0.neibExtreme[1] = p.indices[1];
		edge0.bumper = bc0;

		Bumper edge3 {};
		edge3.shapeType = Bumper::CAPSULOID;
		bc3.neibExtreme[0] = p.indices[3];
		bc3.neibExtreme[1] = p.indices[0];
		edge3.bumper = bc3;

		Bumper edge1 {};
		edge1.shapeType = Bumper::CAPSULOID;
		bc1.neibExtreme[0] = p.indices[1];
		bc1.neibExtreme[1] = p.indices[2];
		edge1.bumper = bc1;

		Bumper edge2 {};
		edge2.shapeType = Bumper::CAPSULOID;
		bc2.neibExtreme[0] = p.indices[2];
		bc2.neibExtreme[1] = p.indices[3];
		edge2.bumper = bc2;

		bumper.push_back(nodeUpper);
		bumper.push_back(nodeLower);

		int upperIndex = bumper.size() - 2;
		int lowerIndex = bumper.size() - 1;

		if (capsuloidAdj[p.indices[0]][p.indices[1]] == -1 || capsuloidAdj[p.indices[1]][p.indices[0]] == -1)
		{
			capsuloidAdj[p.indices[0]][p.indices[1]] = bumper.size();
			capsuloidAdj[p.indices[1]][p.indices[0]] = bumper.size();

			bumper.push_back(edge0);
		}

		int idx = capsuloidAdj[p.indices[0]][p.indices[1]];

		std::get<BumperCapsuloid>(bumper[idx].bumper).neibBumpers.push_back(upperIndex);
		std::get<BumperCapsuloid>(bumper[idx].bumper).neibBumpers.push_back(lowerIndex);

		std::get<BumperQuad>(bumper[upperIndex].bumper).neibSide[0] = idx;
		std::get<BumperQuad>(bumper[lowerIndex].bumper).neibSide[0] = idx;

		if (capsuloidAdj[p.indices[1]][p.indices[2]] == -1 || capsuloidAdj[p.indices[2]][p.indices[1]] == -1)
		{
			capsuloidAdj[p.indices[1]][p.indices[2]] = bumper.size();
			capsuloidAdj[p.indices[2]][p.indices[1]] = bumper.size();

			bumper.push_back(edge1);
		}

		idx = capsuloidAdj[p.indices[1]][p.indices[2]];

		std::get<BumperCapsuloid>(bumper[idx].bumper).neibBumpers.push_back(upperIndex);
		std::get<BumperCapsuloid>(bumper[idx].bumper).neibBumpers.push_back(lowerIndex);

		std::get<BumperQuad>(bumper[upperIndex].bumper).neibSide[1] = idx;
		std::get<BumperQuad>(bumper[lowerIndex].bumper).neibSide[1] = idx;

		if (capsuloidAdj[p.indices[2]][p.indices[3]] == -1 || capsuloidAdj[p.indices[3]][p.indices[2]] == -1)
		{
			capsuloidAdj[p.indices[2]][p.indices[3]] = bumper.size();
			capsuloidAdj[p.indices[3]][p.indices[2]] = bumper.size();

			bumper.push_back(edge2);
		}

		idx = capsuloidAdj[p.indices[2]][p.indices[3]];

		std::get<BumperCapsuloid>(bumper[idx].bumper).neibBumpers.push_back(upperIndex);
		std::get<BumperCapsuloid>(bumper[idx].bumper).neibBumpers.push_back(lowerIndex);

		std::get<BumperQuad>(bumper[upperIndex].bumper).neibSide[2] = idx;
		std::get<BumperQuad>(bumper[lowerIndex].bumper).neibSide[2] = idx;

		if (capsuloidAdj[p.indices[3]][p.indices[0]] == -1 || capsuloidAdj[p.indices[0]][p.indices[3]] == -1)
		{
			capsuloidAdj[p.indices[3]][p.indices[0]] = bumper.size();
			capsuloidAdj[p.indices[0]][p.indices[3]] = bumper.size();

			bumper.push_back(edge3);
		}

		idx = capsuloidAdj[p.indices[3]][p.indices[0]];

		std::get<BumperCapsuloid>(bumper[idx].bumper).neibBumpers.push_back(upperIndex);
		std::get<BumperCapsuloid>(bumper[idx].bumper).neibBumpers.push_back(lowerIndex);

		std::get<BumperQuad>(bumper[upperIndex].bumper).neibSide[3] = idx;
		std::get<BumperQuad>(bumper[lowerIndex].bumper).neibSide[3] = idx;
	}
}

FourSpheres BumperGraph::computeFourSpheresFrom(const Quadrilateral& quad, int direction) const
{
	Sphere s0 = sphere[quad.indices[0]];
	Sphere s1 = sphere[quad.indices[1]];
	Sphere s2 = sphere[quad.indices[2]];
	Sphere s3 = sphere[quad.indices[3]];

	glm::vec3 diagonalAB = s0.center - s2.center;
	glm::vec3 diagonalCD = s1.center - s3.center;
	glm::vec3 midPlaneNormal = normalize(cross(diagonalAB, diagonalCD));

	Plane midPlane {midPlaneNormal, (s0.center + s1.center + s2.center + s3.center) / 4.0f};

	s0.center = midPlane.project(s0.center);
	s1.center = midPlane.project(s1.center);
	s2.center = midPlane.project(s2.center);
	s3.center = midPlane.project(s3.center);

	if (direction < 0) midPlane.flip();

	Plane upperPlane = commonTangentPlane(s0, s1, s2, s3, direction);

	s0.radius = distance(s0.center, upperPlane.project(s0.center));
	s1.radius = distance(s1.center, upperPlane.project(s1.center));
	s2.radius = distance(s2.center, upperPlane.project(s2.center));
	s3.radius = distance(s3.center, upperPlane.project(s3.center));

	FourSpheres original = {sphere[quad.indices[0]], sphere[quad.indices[1]], sphere[quad.indices[2]], sphere[quad.indices[3]]};
	FourSpheres jq = {s0, s1, s2, s3};

	jq.indices[0] = quad.indices[0];
	jq.indices[1] = quad.indices[1];
	jq.indices[2] = quad.indices[2];
	jq.indices[3] = quad.indices[3];
	jq.midPlane = midPlane;
	jq.upperPlane = upperPlane;

	jq.computeError(original);
	if (!upperPlane.valid)
		jq.error = FLT_MAX;

	return jq;
}

void BumperGraph::initializeBumperNodes(const SphereMesh &sm)
{
	initializeBumperSpheres();
	std::vector capsuloidAdj (sphere.size(), std::vector(sphere.size(), -1));
	initializeBumperCapsuloids(sm, capsuloidAdj);
	initializeBumperPrysmoids(sm, capsuloidAdj);
	initializeBumperQuads(sm, capsuloidAdj);

	sortByType();
}

void BumperGraph::constructFrom (const SphereMesh &sm)
{
	sphere = sm.spheres;

	initializeBumperNodes(sm);
	floodfill();
}

bool BumperGraph::pushOutsideSphere(glm::vec3& p, glm::vec3& n, const int& bumperIndex) const
{
	const Bumper& node = bumper[bumperIndex];
	auto [sphereIndex, cones] = std::get<BumperSphere>(node.bumper);
	const Sphere& s = sphere[sphereIndex];

	const glm::vec3 v = p - s.center;
	if (dot(v, v) > s.radius * s.radius)
		return false;

	n = normalize(v);
	p = s.center + n * (s.radius + PUSH_EPSILON);

	return true;
}

float BumperGraph::signedDistanceFromSphere(const int bumperIndex, const glm::vec3& p, glm::vec3& closestPos, glm::vec3& closestNorm)
const
{
	const Bumper& node = bumper[bumperIndex];
	auto [sphereIndex, cones] = std::get<BumperSphere>(node.bumper);
	const Sphere& s = sphere[sphereIndex];

	const float currDist = length(p - s.center);

	closestNorm = (p - s.center) / currDist;
	closestPos = s.center + closestNorm * s.radius;

	return currDist - s.radius;
}

float BumperGraph::closestSphereOn(const glm::vec3& p, const BumperCapsuloid& bc) const
{
	const Sphere& a = sphere[bc.sphereIndex[0]];

	const glm::vec3 pa = p - a.center;
	float t = dot(pa, bc.dNorm);
	const float pq = std::sqrt(dot(pa, pa) - t * t);

	t +=  bc.k * pq;

	return glm::clamp(t, 0.0f, bc.dLength);
}

Sphere BumperGraph::getInterpolatedSphere(const BumperCapsuloid& bc, const float t) const
{
	Sphere s{};
	const Sphere& a = sphere[bc.sphereIndex[0]];

	s.center = a.center + t * bc.dNorm;
	s.radius = a.radius + t * bc.dR;

	return s;
}

bool BumperGraph::pushOutsideCapsuloid (glm::vec3 &p, glm::vec3 &n, int &bumperIndex) const
{
	const Bumper& node = bumper[bumperIndex];
	const BumperCapsuloid bc = std::get<BumperCapsuloid>(node.bumper);

	const float t = closestSphereOn(p, bc);

	if (t == 0)
	{
		const int index = bc.sphereIndex[0];
		if (!pushOutsideSphere(p, n, index))
			return false;

		bumperIndex = index;
		return true;
	}

	if (t == bc.dLength)
	{
		const int index = bc.sphereIndex[1];
		if(!pushOutsideSphere(p, n, index))
			return false;

		bumperIndex = index;
		return true;
	}

	const Sphere s = getInterpolatedSphere(bc, t);
	const glm::vec3 pqPrime = p - s.center;
	const float sqrdDist = dot(pqPrime, pqPrime);

	if (sqrdDist >= s.radius * s.radius)
		return false;

	n = pqPrime / std::sqrt(sqrdDist);
	p = s.center + n * (s.radius + PUSH_EPSILON);

	return true;
}

float BumperGraph::signedDistanceFromCapsuloid(const int bumperIndex, const glm::vec3& p, glm::vec3& closestPos, glm::vec3&
                                               closestNorm) const
{
	const Bumper& node = bumper[bumperIndex];
	const BumperCapsuloid bc = std::get<BumperCapsuloid>(node.bumper);

	const float t = closestSphereOn(p, bc);
	const Sphere s = getInterpolatedSphere(bc, t);
	const float dist = distance(p, s.center);

	closestNorm = (p - s.center) / dist;
	closestPos = s.center + closestNorm * s.radius;

	return dist - s.radius;
}

bool BumperGraph::pushOutsidePrysmoid (glm::vec3 &p, glm::vec3 &n, int &bumperIndex) const // NOLINT(*-no-recursion)
{
	const Bumper& node = bumper[bumperIndex];
	const BumperPrysmoid bp = std::get<BumperPrysmoid>(node.bumper);

	if (bp.midPlane.isBehind(p))
	{
		bumperIndex = bp.neibOpp;
		return pushOutsidePrysmoid(p, n, bumperIndex);
	}

	if (bp.m1.isBehind(p))
	{
		bumperIndex = bp.neibSide[1];
		return pushOutsideCapsuloid(p, n, bumperIndex);
	}
	if (bp.m2.isBehind(p))
	{
		bumperIndex = bp.neibSide[2];
		return pushOutsideCapsuloid(p, n, bumperIndex);
	}
	if (bp.m0.isBehind(p))
	{
		bumperIndex = bp.neibSide[0];
		return pushOutsideCapsuloid(p, n, bumperIndex);
	}

	const float d = dot(p, normalize(bp.upperPlane.n)) - bp.upperPlane.k;

	if (d > 0) return false;

	n = bp.upperPlane.n;
	p -= n * (d - PUSH_EPSILON);

	return true;
}

bool BumperGraph::pushOutsideQuad(glm::vec3 &p, glm::vec3 &n, int &bumperIndex) const
{
	const Bumper& node = bumper[bumperIndex];
	const BumperQuad bq = std::get<BumperQuad>(node.bumper);

	if (bq.midPlane.isBehind(p))
	{
		bumperIndex = bq.neibOpp;
		return pushOutsidePrysmoid(p, n, bumperIndex);
	}

	if (bq.sidePlanes[0].isBehind(p))
	{
		bumperIndex = bq.neibSide[0];
		return pushOutsideCapsuloid(p, n, bumperIndex);
	}
	if (bq.sidePlanes[1].isBehind(p))
	{
		bumperIndex = bq.neibSide[1];
		return pushOutsideCapsuloid(p, n, bumperIndex);
	}
	if (bq.sidePlanes[2].isBehind(p))
	{
		bumperIndex = bq.neibSide[2];
		return pushOutsideCapsuloid(p, n, bumperIndex);
	}
	if (bq.sidePlanes[3].isBehind(p))
	{
		bumperIndex = bq.neibSide[3];
		return pushOutsideCapsuloid(p, n, bumperIndex);
	}

	const float d = dot(p, normalize(bq.upperPlane.n)) - bq.upperPlane.k;

	if (d > 0) return false;

	n = bq.upperPlane.n;
	p -= n * (d - PUSH_EPSILON);

	return true;
}

float BumperGraph::signedDistanceFromPrysmoid(const int bumperIndex, const glm::vec3& p, glm::vec3& closestPos, glm::vec3& // NOLINT(*-no-recursion)
                                              closestNorm) const
{
	const Bumper& node = bumper[bumperIndex];
	const BumperPrysmoid bp = std::get<BumperPrysmoid>(node.bumper);

	if (bp.midPlane.isBehind(p)) return signedDistanceFromPrysmoid(bp.neibOpp, p, closestPos, closestNorm);

	if (bp.m1.isBehind(p)) return signedDistanceFromCapsuloid(bp.neibSide[1], p, closestPos,  closestNorm);
	if (bp.m2.isBehind(p)) return signedDistanceFromCapsuloid(bp.neibSide[2], p, closestPos,  closestNorm);
	if (bp.m0.isBehind(p)) return signedDistanceFromCapsuloid(bp.neibSide[0], p, closestPos,  closestNorm);

	closestNorm = bp.upperPlane.n;
	closestPos = bp.upperPlane.project(p);

	return bp.upperPlane.distance(p);
}

float BumperGraph::signedDistanceFromQuad(const int bumperIndex, const glm::vec3 &p, glm::vec3 &closestPos,
                                          glm::vec3 &closestNorm) const
{
	const Bumper& node = bumper[bumperIndex];
	const BumperQuad bq = std::get<BumperQuad>(node.bumper);

	if (bq.midPlane.isBehind(p)) return signedDistanceFromPrysmoid(bq.neibOpp, p, closestPos, closestNorm);

	if (bq.sidePlanes[0].isBehind(p)) return signedDistanceFromCapsuloid(bq.neibSide[0], p, closestPos,  closestNorm);
	if (bq.sidePlanes[1].isBehind(p)) return signedDistanceFromCapsuloid(bq.neibSide[1], p, closestPos,  closestNorm);
	if (bq.sidePlanes[2].isBehind(p)) return signedDistanceFromCapsuloid(bq.neibSide[2], p, closestPos,  closestNorm);
	if (bq.sidePlanes[3].isBehind(p)) return signedDistanceFromCapsuloid(bq.neibSide[3], p, closestPos,  closestNorm);

	closestNorm = bq.upperPlane.n;
	closestPos = bq.upperPlane.project(p);

	return bq.upperPlane.distance(p);
}

std::pair<int, float> BumperGraph::sampleSignedDistanceFromSphere (int bumperIndex, const glm::vec3 &p) const
{
	const Bumper& node = bumper[bumperIndex];
	auto [sphereIndex, cones] = std::get<BumperSphere>(node.bumper);
	const Sphere& s = sphere[sphereIndex];

	const float currDist = length(p - s.center);

	return {bumperIndex, currDist - s.radius};
}

std::pair<int, float> BumperGraph::sampleSignedDistanceFromCapsuloid (int bumperIndex, const glm::vec3 &p) const
{
	const Bumper& node = bumper[bumperIndex];
	const BumperCapsuloid bc = std::get<BumperCapsuloid>(node.bumper);

	const float t = closestSphereOn(p, bc);

	if (t == 0)
		return sampleSignedDistanceFromSphere(bc.neibExtreme[0], p);
	if (t == bc.dLength)
		return sampleSignedDistanceFromSphere(bc.neibExtreme[1], p);

	const Sphere s = getInterpolatedSphere(bc, t);
	const float dist = distance(p, s.center);

	return {bumperIndex, dist - s.radius};
}

std::pair<int, float> BumperGraph::sampleSignedDistanceFromPrysmoid (int bumperIndex, const glm::vec3 &p) const
{
	const Bumper& node = bumper[bumperIndex];
	const BumperPrysmoid bp = std::get<BumperPrysmoid>(node.bumper);

	if (bp.midPlane.isBehind(p))
		return {-1, FLT_MAX};

	if (bp.m1.isBehind(p)) return sampleSignedDistanceFromCapsuloid(bp.neibSide[1], p);
	if (bp.m2.isBehind(p)) return sampleSignedDistanceFromCapsuloid(bp.neibSide[2], p);
	if (bp.m0.isBehind(p)) return sampleSignedDistanceFromCapsuloid(bp.neibSide[0], p);

	return {bumperIndex, bp.upperPlane.distance(p)};
}

std::pair<int, float> BumperGraph::sampleSignedDistanceFromQuad(int bumperIndex, const glm::vec3 &p) const
{
	const Bumper& node = bumper[bumperIndex];
	const BumperQuad bq = std::get<BumperQuad>(node.bumper);

	if (bq.midPlane.isBehind(p)) return {-1, FLT_MAX};

	if (bq.sidePlanes[0].isBehind(p)) return sampleSignedDistanceFromCapsuloid(bq.neibSide[0], p);
	if (bq.sidePlanes[1].isBehind(p)) return sampleSignedDistanceFromCapsuloid(bq.neibSide[1], p);
	if (bq.sidePlanes[2].isBehind(p)) return sampleSignedDistanceFromCapsuloid(bq.neibSide[2], p);
	if (bq.sidePlanes[3].isBehind(p)) return sampleSignedDistanceFromCapsuloid(bq.neibSide[2], p);

	return {bumperIndex, bq.upperPlane.distance(p)};
}

std::vector<std::vector<bool>> BumperGraph::getBumperAdjMatrix() const
{
	std::vector sphereAdj (sphere.size(), std::vector(sphere.size(), false));

	for (auto& b : bumper)
	{
		if (b.shapeType == Bumper::SPHERE) continue;

		if (b.shapeType == Bumper::CAPSULOID)
		{
			auto& bc = std::get<BumperCapsuloid>(b.bumper);
			sphereAdj[bc.sphereIndex[0]][bc.sphereIndex[1]] = true;
			sphereAdj[bc.sphereIndex[1]][bc.sphereIndex[0]] = true;
		}

		if (b.shapeType == Bumper::PRYSMOID)
		{
			auto& bp = std::get<BumperPrysmoid>(b.bumper);

			sphereAdj[bp.sphereIndex[0]][bp.sphereIndex[1]] = true;
			sphereAdj[bp.sphereIndex[1]][bp.sphereIndex[0]] = true;

			sphereAdj[bp.sphereIndex[0]][bp.sphereIndex[2]] = true;
			sphereAdj[bp.sphereIndex[2]][bp.sphereIndex[0]] = true;

			sphereAdj[bp.sphereIndex[1]][bp.sphereIndex[2]] = true;
			sphereAdj[bp.sphereIndex[2]][bp.sphereIndex[1]] = true;
		}

		if (b.shapeType == Bumper::QUAD)
		{
			auto& bq = std::get<BumperQuad>(b.bumper);

			sphereAdj[bq.sphereIndex[0]][bq.sphereIndex[1]] = true;
			sphereAdj[bq.sphereIndex[1]][bq.sphereIndex[0]] = true;

			sphereAdj[bq.sphereIndex[0]][bq.sphereIndex[3]] = true;
			sphereAdj[bq.sphereIndex[3]][bq.sphereIndex[0]] = true;

			sphereAdj[bq.sphereIndex[1]][bq.sphereIndex[2]] = true;
			sphereAdj[bq.sphereIndex[2]][bq.sphereIndex[1]] = true;

			sphereAdj[bq.sphereIndex[2]][bq.sphereIndex[3]] = true;
			sphereAdj[bq.sphereIndex[3]][bq.sphereIndex[2]] = true;
		}
	}

	return sphereAdj;
}

void BumperGraph::fillBranchWithTransform(const int value, const std::vector<std::vector<bool>> &adjMatrix,
	std::queue<int>& q, std::unordered_set<int>& visited)
{
	q.push(value);
	visited.insert(value);

	while (!q.empty())
	{
		const int top = q.front();
		q.pop();

		rig[top] = value;

		for (int i = 0; i < sphere.size(); i++)
			if (i != top && adjMatrix[top][i] && !visited.contains(i))
			{
				visited.insert(i);
				q.push(i);
			}
	}
}

void BumperGraph::floodfill()
{
	transformMapper[safeLeft] = {{0.0f, 0.0f, 1.0f}, 0.5f, sphere[safeLeft].center};
	transformMapper[safeLeft] = {{0.0f, 0.0f, 1.0f}, -0.5f, sphere[safeRight].center};
	transformMapper[cancer] = {{0.0f, 0.0f, 1.0f}, 0.0f, sphere[cancer].center};

	// Sphere floodfill
	const auto adjMatrix = getBumperAdjMatrix();

	std::queue<int> q;
	std::unordered_set<int> visited;

	fillBranchWithTransform(cancer, adjMatrix, q, visited);
	fillBranchWithTransform(safeLeft, adjMatrix, q, visited);
	fillBranchWithTransform(safeRight, adjMatrix, q, visited);

	// Bumpers floodfill
	for (int i = 0; i < bumper.size(); i++)
	{
		auto& b = bumper[i];

		if (b.shapeType == Bumper::SPHERE) continue;
		if (b.shapeType == Bumper::CAPSULOID)
		{
			const auto& bc = std::get<BumperCapsuloid>(b.bumper);
			const int index0 = rig[bc.sphereIndex[0]];
			const int index1 = rig[bc.sphereIndex[1]];

			if (index0 == index1 && index0 != cancer) rig[i] = index0;
		}
		if (b.shapeType == Bumper::PRYSMOID)
		{
			const auto& bp = std::get<BumperPrysmoid>(b.bumper);
			const int index0 = rig[bp.sphereIndex[0]];
			const int index1 = rig[bp.sphereIndex[1]];
			const int index2 = rig[bp.sphereIndex[2]];

			if (index0 == index1 == index2 && index0 != cancer) rig[i] = index0;
		}
		if (b.shapeType == Bumper::QUAD)
		{
			const auto& bq = std::get<BumperPrysmoid>(b.bumper);
			const int index0 = rig[bq.sphereIndex[0]];
			const int index1 = rig[bq.sphereIndex[1]];
			const int index2 = rig[bq.sphereIndex[2]];
			const int index3 = rig[bq.sphereIndex[3]];

			if (index0 == index1 == index2 == index3 && index0 != cancer) rig[i] = index0;
		}
	}
}

float BumperGraph::projectOnBruteForce(glm::vec3& p, glm::vec3& n, int& bumperIndex)  const
{
	int counter = 50;

	while (true)
	{
		glm::vec3 bestPos = p;
		glm::vec3 bestNorm = n;
		int bumpIdx = bumperIndex;
		float minDistance = FLT_MAX;

		for (int i = 0; i < bumper.size(); i++)
		{
			glm::vec3 newPos, newNorm;
			float currDist = FLT_MAX;
			bool sideUp = false;
			switch (bumper[i].shapeType)
			{
				case Bumper::SPHERE:
					currDist = signedDistanceFromSphere(i, p, newPos, newNorm);
					break;
				case Bumper::CAPSULOID:
					currDist = signedDistanceFromCapsuloid(i, p, newPos, newNorm);
					break;
				case Bumper::PRYSMOID:
					sideUp = isPointOverPrysmoid(i, p);
					if (!sideUp) i++; // Skipping the upper plane
					currDist = signedDistanceFromPrysmoid(i, p, newPos, newNorm);
					if (sideUp) i++; // Skipping the lower plane
					break;
				case Bumper::QUAD:
					sideUp = isPointOverQuad(i, p);
					if (!sideUp) i++; // Skipping the upper plane
					currDist = signedDistanceFromQuad(i, p, newPos, newNorm);
					if (sideUp) i++; // Skipping the lower plane
					break;
			}
			if (currDist < minDistance)
			{
				bestPos = newPos;
				bestNorm = newNorm;
				bumpIdx = i;
				minDistance = currDist;
			}
		}

		p = bestPos;
		n = bestNorm;
		bumperIndex = bumpIdx;

		if (counter-- <= 0 || minDistance > -0.001f)
			return minDistance;
	}
}

bool BumperGraph::pushOutsideBruteForce(glm::vec3& p, glm::vec3& n, const int& bumperIndex) const
{
	glm::vec3 tmpPos = p;
	glm::vec3 tmpNorm = n;
	int tmpBumperIndex = bumperIndex;
	projectOnBruteForce(tmpPos, tmpNorm, tmpBumperIndex);

	if (bumper[tmpBumperIndex].shapeType == Bumper::SPHERE)	return pushOutsideSphere(p, n, tmpBumperIndex);
	if (bumper[tmpBumperIndex].shapeType == Bumper::CAPSULOID) return pushOutsideCapsuloid(p, n, tmpBumperIndex);
	if (bumper[tmpBumperIndex].shapeType == Bumper::PRYSMOID)	return pushOutsidePrysmoid(p, n, tmpBumperIndex);
	if (bumper[tmpBumperIndex].shapeType == Bumper::QUAD)		return pushOutsideQuad(p, n, tmpBumperIndex);

	return false;
}

inline int getRandomInVector(const std::vector<int> & v, const int randSeed){
	if (v.empty()) return -1;
	return v[ randSeed % v.size() ];
}

bool BumperGraph::projectOn (glm::vec3 &p, glm::vec3 &n, int &bumperIndex) const
{
	projectOnBruteForce(p, n, bumperIndex);
	return true;
}

void BumperGraph::animateRig(int iterations)
{
	for (int it = 0; it < iterations; ++it)
	{
		for (int i = 0; i < bumper.size(); i++)
		{
			auto& b = bumper[i];

			if (b.shapeType == Bumper::SPHERE && rig[i] != cancer)
			{
				auto [axis, angle, pivot] = transformMapper[rig[i]];

				glm::vec3 translatedPos = sphere[i].center - pivot;
				glm::mat3 rotationMatrix = glm::mat3(glm::rotate(glm::mat4(1.0f), angle, axis));
				glm::vec3 rotatedPos = rotationMatrix * translatedPos;

				sphere[i].center = rotatedPos + pivot;

				auto& bn = std::get<BumperSphere>(b.bumper);

				bn.translate(-pivot);
				bn.rotate(rotationMatrix);
				bn.translate(pivot);
			}
			if (b.shapeType == Bumper::CAPSULOID && rig[i] != cancer)
			{
				auto [axis, angle, pivot] = transformMapper[rig[i]];

				glm::vec3 translatedPos = sphere[i].center - pivot;
				glm::mat3 rotationMatrix = glm::mat3(glm::rotate(glm::mat4(1.0f), angle, axis));
				glm::vec3 rotatedPos = rotationMatrix * translatedPos;

				auto& bc = std::get<BumperCapsuloid>(b.bumper);

				bc.translate(-pivot);
				bc.rotate(rotationMatrix);
				bc.translate(pivot);
			}
			if (b.shapeType == Bumper::PRYSMOID && rig[i] != cancer)
			{
				auto [axis, angle, pivot] = transformMapper[rig[i]];

				glm::vec3 translatedPos = sphere[i].center - pivot;
				glm::mat3 rotationMatrix = glm::mat3(glm::rotate(glm::mat4(1.0f), angle, axis));
				glm::vec3 rotatedPos = rotationMatrix * translatedPos;

				auto& bp = std::get<BumperPrysmoid>(b.bumper);

				bp.translate(-pivot);
				bp.rotate(rotationMatrix);
				bp.translate(pivot);
			}
			if (b.shapeType == Bumper::QUAD && rig[i] != cancer)
			{
				auto [axis, angle, pivot] = transformMapper[rig[i]];

				glm::vec3 translatedPos = sphere[i].center - pivot;
				glm::mat3 rotationMatrix = glm::mat3(glm::rotate(glm::mat4(1.0f), angle, axis));
				glm::vec3 rotatedPos = rotationMatrix * translatedPos;

				auto& bq = std::get<BumperQuad>(b.bumper);

				bq.translate(-pivot);
				bq.rotate(rotationMatrix);
				bq.translate(pivot);
			}
		}
	}
}

bool BumperGraph::isPointOverPrysmoid (const int bumperIndex, const glm::vec3 &p) const
{
	const Bumper& node = bumper[bumperIndex];
	const BumperPrysmoid bp = std::get<BumperPrysmoid>(node.bumper);

	const Sphere& a = sphere[bp.sphereIndex[0]];
	const glm::vec3 pa = p - a.center;

	return dot(bp.upperPlane.n, pa) > 0; // If the point is over the prysmoid
}

bool BumperGraph::isPointOverQuad(const int bumperIndex, const glm::vec3 &p) const
{
	const Bumper& node = bumper[bumperIndex];
	const BumperQuad bq = std::get<BumperQuad>(node.bumper);

	const Sphere& a = sphere[bq.sphereIndex[0]];
	const glm::vec3 pa = p - a.center;

	return dot(bq.upperPlane.n, pa) > 0; // If the point is over the quad
}

void BumperQuad::translate(const glm::vec3 &t)
{
	upperPlane.translate(t);
	midPlane.translate(t);

	for (auto& sideP : sidePlanes)
		sideP.translate(t);
}

void BumperQuad::scale(float r)
{
	upperPlane.scale(r);
	midPlane.scale(r);

	for (auto& sideP : sidePlanes)
		sideP.scale(r);
}

void BumperQuad::rotate(const glm::mat3 &rot)
{
	upperPlane.rotate(rot);
	midPlane.rotate(rot);

	for (auto& sideP : sidePlanes)
		sideP.rotate(rot);
}

void BumperQuad::rotateAround(const glm::vec3& axis, const float angle)
{
	const auto rot = glm::mat3(glm::rotate(glm::mat4(1.0f), angle, axis));
	rotate(rot);
}

BumperQuad::BumperQuad(const FourSpheres &fs)
{
	upperPlane = fs.upperPlane;
	midPlane = fs.midPlane;

	sphereIndex[0] = fs.indices[0];
	sphereIndex[1] = fs.indices[1];
	sphereIndex[2] = fs.indices[2];
	sphereIndex[3] = fs.indices[3];

	auto bp = BumperPrysmoid(sphereIndex[0], sphereIndex[1], sphereIndex[2], fs.spheres[0], fs.spheres[1], fs.spheres[2], 1);
	auto bp1 = BumperPrysmoid(sphereIndex[0], sphereIndex[2], sphereIndex[3], fs.spheres[0], fs.spheres[2], fs.spheres[3], 1);

	getSidePlanes(bp, bp1);
}

void BumperQuad::getSidePlanes(const BumperPrysmoid &bp, const BumperPrysmoid &bp1)
{
	sidePlanes[0] = bp.m2;
	sidePlanes[1] = bp.m0;
	sidePlanes[2] = bp1.m0;
	sidePlanes[3] = bp1.m1;
}

bool Bumper::hasAParent () const
{
	if (shapeType == SPHERE)
		return true;
	if (shapeType == CAPSULOID)
		return std::get<BumperCapsuloid>(bumper).hasFather;

	return false;
}

std::string Bumper::serialize() const
{
	std::ostringstream oss;

	if (shapeType == SPHERE)
		return "";

	if (shapeType == CAPSULOID)
	{
		auto bc = std::get<BumperCapsuloid>(bumper);
		oss << bc.sphereIndex[0]
		<< " " << bc.sphereIndex[1]
		<< " " << bc.dNorm.x
		<< " " << bc.dNorm.y
		<< " " << bc.dNorm.z
		<< " " << bc.dLength
		<< " " << bc.dR
		<< " " << bc.k;
	}
	else if (shapeType == PRYSMOID)
	{
		auto bp = std::get<BumperPrysmoid>(bumper);
		oss << bp.sphereIndex[0]
		<< " " << bp.sphereIndex[1]
		<< " " << bp.sphereIndex[2]
		<< " " << bp.neibSide[0]
		<< " " << bp.neibSide[1]
		<< " " << bp.neibSide[2]
		<< " " << bp.neibOpp
		<< " " << bp.m1.serialize()
		<< " " << bp.m2.serialize()
		<< " " << bp.m0.serialize()
		<< " " << bp.midPlane.serialize()
		<< " " << bp.upperPlane.serialize();
	}
	else
		throw std::runtime_error("Unknown shape type encountered during serialization.");

	return oss.str();
}

void Bumper::translate(const glm::vec3 &t)
{
	if (shapeType == SPHERE) std::get<BumperSphere>(bumper).translate(t);
	if (shapeType == CAPSULOID) std::get<BumperCapsuloid>(bumper).translate(t);
	if (shapeType == PRYSMOID) std::get<BumperPrysmoid>(bumper).translate(t);
	if (shapeType == QUAD) std::get<BumperQuad>(bumper).translate(t);
}

void Bumper::scale(float r)
{
	if (shapeType == SPHERE) std::get<BumperSphere>(bumper).scale(r);
	if (shapeType == CAPSULOID) std::get<BumperCapsuloid>(bumper).scale(r);
	if (shapeType == PRYSMOID) std::get<BumperPrysmoid>(bumper).scale(r);
	if (shapeType == QUAD) std::get<BumperQuad>(bumper).scale(r);
}

void Bumper::rotate(const glm::mat3 &rot)
{
	if (shapeType == SPHERE) std::get<BumperSphere>(bumper).rotate(rot);
	if (shapeType == CAPSULOID) std::get<BumperCapsuloid>(bumper).rotate(rot);
	if (shapeType == PRYSMOID) std::get<BumperPrysmoid>(bumper).rotate(rot);
	if (shapeType == QUAD) std::get<BumperQuad>(bumper).rotate(rot);
}

void Bumper::rotateAround(const glm::vec3& axis, const float angle)
{
	if (shapeType == SPHERE) std::get<BumperSphere>(bumper).rotateAround(axis, angle);
	if (shapeType == CAPSULOID) std::get<BumperCapsuloid>(bumper).rotateAround(axis, angle);
	if (shapeType == PRYSMOID) std::get<BumperPrysmoid>(bumper).rotateAround(axis, angle);
	if (shapeType == QUAD) std::get<BumperQuad>(bumper).rotateAround(axis, angle);
}

bool Bumper::operator == (const Bumper& bn) const
{
	if (shapeType != bn.shapeType) return false;

	switch (shapeType)
	{
		case SPHERE: return std::get<BumperSphere>(bumper).sphereIndex == std::get<BumperSphere>(bn.bumper).sphereIndex;
		case CAPSULOID: return std::get<BumperCapsuloid>(bumper) == std::get<BumperCapsuloid>(bn.bumper);
		case PRYSMOID: return std::get<BumperPrysmoid>(bumper) == std::get<BumperPrysmoid>(bn.bumper);
		default: return false;
	}

	return false;
}

void BumperGraph::sortByType()
{
	std::vector<int> permutation (bumper.size());

	for (int i = 0; i < bumper.size(); i++)
		permutation[i] = i;

	std::ranges::sort(permutation, [this](const int a, const int b) {
		if (bumper[a].shapeType < bumper[b].shapeType) return true;
		if (bumper[a].shapeType > bumper[b].shapeType) return false;
		return a < b;
	});

	std::vector<int> inversePermutation (bumper.size());
	for (int i = 0; i < bumper.size(); i++)
		inversePermutation[permutation[i]] = i;

	for (auto & bn : bumper)
	{
		if (bn.shapeType == Bumper::CAPSULOID)
		{
			auto& bc = std::get<BumperCapsuloid>(bn.bumper);

			bc.neibExtreme[0] = inversePermutation[bc.neibExtreme[0]];
			bc.neibExtreme[1] = inversePermutation[bc.neibExtreme[1]];
		}
		else if (bn.shapeType == Bumper::PRYSMOID)
		{
			auto& bp = std::get<BumperPrysmoid>(bn.bumper);

			bp.neibSide[0] = inversePermutation[bp.neibSide[0]];
			bp.neibSide[1] = inversePermutation[bp.neibSide[1]];
			bp.neibSide[2] = inversePermutation[bp.neibSide[2]];

			bp.neibOpp = inversePermutation[bp.neibOpp];
		}
	}

	const std::vector<Bumper> oldBumperNode = bumper;
	for (const int i : permutation)
		bumper[i] = oldBumperNode[permutation[i]];
}

void BumperGraphMesh::appendSubMesh(SubMesh &dest, const SubMesh &src) {
	const int offset = dest.vertices.size();

    for (const auto &v : src.vertices)
        dest.vertices.push_back(v);
    for (const auto &f : src.faces)
        dest.faces.push_back({ f[0] + offset, f[1] + offset, f[2] + offset });
}

void SubMesh::addQuad(const int a, const int b, const int c, const int d)
{
	addTriangle(a, b, c);
	addTriangle(c, d, a);
}

void SubMesh::addTriangle(const int a, const int b, const int c)
{
	if (a == b || b == c || a == c) return;

	faces.push_back({a, b, c});
}

SubMesh BumperGraphMesh::createSphereSubMesh(const Sphere &s, const int resolution)
{
    SubMesh sub;
	const int n = resolution;
	const int m = resolution * 2;
    sub.name = "sphere";

    for (int i = 1; i < n; ++i) {
	    const float parallel = glm::pi<float>() * i / n;
        for (int j = 0; j < m; ++j) {
	        const float meridian = 2.0f * glm::pi<float>() * j / m;
	        const float x = sin(parallel) * cos(meridian);
	        const float y = cos(parallel);
	        const float z = sin(parallel) * sin(meridian);
            sub.vertices.push_back(s.center + s.radius * glm::vec3(x, y, z));
        }
    }
	sub.vertices.push_back(s.center + s.radius * glm::vec3(0, 1, 0));
	sub.vertices.push_back(s.center + s.radius * glm::vec3(0, -1, 0));
    const int nPole = sub.vertices.size() - 2, sPole = sub.vertices.size() - 1;

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
        	int a = (i - 1 + 0) * m + (j + 0) % m;
        	int b = (i - 1 + 1) * m + (j + 0) % m;
        	int c = (i - 1 + 1) * m + (j + 1) % m;
        	int d = (i - 1 + 0) * m + (j + 1) % m;

        	if (i == 0)
        		a = d = nPole;
        	else if (i == n - 1)
        		b = c = sPole;

        	sub.addQuad(c, b, a, d);
        }
    }

    return sub;
}

SubMesh BumperGraphMesh::createCapsuleSubMesh(const Sphere &s0, const Sphere &s1, int resolution)
{
    SubMesh sub;
    sub.name = "capsule";

    glm::vec3 v0 = s0.center, v1 = s1.center;
    float r0 = s0.radius, r1 = s1.radius;

    if (r1 < r0) std::swap(v0, v1); std::swap(r0, r1);

    glm::vec3 d = v0 - v1;
    float dLength = length(d);
    if (dLength < 1e-6f) return sub;

    float diff = r1 - r0;
    float l = sqrt(dLength * dLength - diff * diff);
    float r0Bis = r0 * (l / dLength);
    float r1Bis = r1 * (l / dLength);
    float d0 = diff * (r0 / dLength);
    float d1 = diff * (r1 / dLength);

    glm::vec3 dir = normalize(d);
    glm::vec3 v0Bis = v0 + dir * d0;
    glm::vec3 v1Bis = v1 + dir * d1;
    glm::vec3 arbitrary = fabs(dir.x) < 0.99f ? glm::vec3(1, 0, 0) : glm::vec3(0, 1, 0);
    glm::vec3 right = normalize(cross(dir, arbitrary));
    glm::vec3 up = normalize(cross(right, dir));
    std::vector<glm::vec3> circle1(resolution), circle2(resolution);

    for (int i = 0; i < resolution; ++i) {
        float theta = 2.0f * glm::pi<float>() * i / resolution;
        glm::vec3 offset1 = (right * cos(theta) + up * sin(theta)) * r0Bis;
        circle1[i] = v0Bis + offset1;
        glm::vec3 offset2 = (right * cos(theta) + up * sin(theta)) * r1Bis;
        circle2[i] = v1Bis + offset2;
    }

    for (int i = 0; i < resolution; ++i) {
        int next = (i + 1) % resolution;
        int i1 = sub.vertices.size() + 1; sub.vertices.push_back(circle1[i]);
        int i2 = sub.vertices.size() + 1; sub.vertices.push_back(circle2[i]);
        int i3 = sub.vertices.size() + 1; sub.vertices.push_back(circle1[next]);
        sub.faces.push_back({ i1, i2, i3 });
        int i4 = sub.vertices.size() + 1; sub.vertices.push_back(circle2[next]);
        sub.faces.push_back({ i2, i4, i3 });
    }
    return sub;
}

SubMesh BumperGraphMesh::createClosedCapsuleSubMesh(const Sphere &s0, const Sphere &s1, int resolution) {
    Sphere sphereA = s0;
    Sphere sphereB = s1;
    glm::vec3 vA = s0.center, vB = s1.center;
    float rA = s0.radius, rB = s1.radius;

    if (rB < rA) {
        std::swap(sphereA, sphereB);
        std::swap(vA, vB);
        std::swap(rA, rB);
    }

    glm::vec3 d_centers = vA - vB;
    float dLength = length(d_centers);
    if (dLength < 1e-6f) return createCapsuleSubMesh(s0, s1, resolution);
    float diff = rB - rA;
    float l = sqrt(dLength * dLength - diff * diff);
    float rABis = rA * (l / dLength);
    float rBBis = rB * (l / dLength);
    float dA = diff * (rA / dLength);
    float dB = diff * (rB / dLength);
    glm::vec3 dir = normalize(d_centers);
    glm::vec3 endpointA = vA + dir * dA;
    glm::vec3 endpointB = vB + dir * dB;

    glm::vec3 arbitrary = fabs(dir.x) < 0.99f ? glm::vec3(1,0,0) : glm::vec3(0,1,0);
    glm::vec3 right = normalize(cross(dir, arbitrary));
    glm::vec3 up = normalize(cross(right, dir));

	SubMesh capsule;
    for (int i = 0; i < resolution; ++i) {
        float theta = 2.0f * glm::pi<float>() * i / resolution;
    	glm::vec3 circAPoint = endpointA + rABis * (right * cos(theta) + up * sin(theta));
    	glm::vec3 circBPoint = endpointB + rBBis * (right * cos(theta) + up * sin(theta));
    	capsule.vertices.push_back(circAPoint);
    	capsule.vertices.push_back(circBPoint);
    }
	capsule.vertices.push_back(sphereA.center);
	capsule.vertices.push_back(sphereB.center);

	int apexA = capsule.vertices.size() - 2, apexB = capsule.vertices.size() - 1;
	for (int i = 0; i < resolution; i++)
	{
		int j = (i + 1) % resolution;
		int a = 2*i;
		int b = 2*i + 1;
		int c = 2*j + 1;
		int d = 2*j;
		capsule.addQuad(a, b, c, d);
		capsule.addTriangle(a, d, apexA);
		capsule.addTriangle(c, b, apexB);
	}

    return capsule;
}

SubMesh BumperGraphMesh::createPlaneSubMesh(const Plane &p, float width, float height)
{
    SubMesh sub;
    sub.name = "plane";
    glm::vec3 n = normalize(p.n);
    glm::vec3 center = n * p.k;
    glm::vec3 tangent = normalize(cross(n, glm::vec3(0,1,0)));
    if (length(tangent) < 1e-3f)
        tangent = normalize(cross(n, glm::vec3(1,0,0)));
    glm::vec3 bitangent = normalize(cross(n, tangent));
    glm::vec3 v0 = center + tangent*(width/2.0f) + bitangent*(height/2.0f);
    glm::vec3 v1 = center - tangent*(width/2.0f) + bitangent*(height/2.0f);
    glm::vec3 v2 = center - tangent*(width/2.0f) - bitangent*(height/2.0f);
    glm::vec3 v3 = center + tangent*(width/2.0f) - bitangent*(height/2.0f);
    glm::vec3 computedNormal = normalize(cross(v1-v0, v2-v0));

    if (dot(computedNormal, n) < 0.0f) {
        sub.vertices.push_back(v0);
        sub.vertices.push_back(v2);
        sub.vertices.push_back(v1);
        sub.vertices.push_back(v3);
        sub.faces.push_back({1,2,3});
        sub.faces.push_back({1,4,2});
    } else {
        sub.vertices.push_back(v0);
        sub.vertices.push_back(v1);
        sub.vertices.push_back(v2);
        sub.vertices.push_back(v3);
        sub.faces.push_back({1,2,3});
        sub.faces.push_back({1,3,4});
    }

    return sub;
}

SubMesh BumperGraphMesh::createSubMeshFromBumper(const Bumper &b, const std::vector<Sphere>& spheres) {
    SubMesh sub;

    std::visit([&]<typename T0>(T0 &&arg) {
        using T = std::decay_t<T0>;
        if constexpr (std::is_same_v<T, BumperSphere>) {
	        const int idx = arg.sphereIndex;

            if (idx >= 0 && idx < static_cast<int>(spheres.size())) {
                sub = createSphereSubMesh(spheres[idx], 128);
                sub.name = "sphere_" + std::to_string(idx);
            }
        } else if constexpr (std::is_same_v<T, BumperCapsuloid>) {
	        const int idx0 = arg.neibExtreme[0] >= 0 ? arg.neibExtreme[0] : arg.sphereIndex[0];
	        const int idx1 = arg.neibExtreme[1] >= 0 ? arg.neibExtreme[1] : arg.sphereIndex[1];

            if (idx0 >= 0 && idx0 < static_cast<int>(spheres.size()) &&
                idx1 >= 0 && idx1 < static_cast<int>(spheres.size())) {
                sub = createClosedCapsuleSubMesh(spheres[idx0], spheres[idx1], 128);
                sub.name = "capsuloid_" + std::to_string(idx0) + "_" + std::to_string(idx1);
            }
        } else if constexpr (std::is_same_v<T, BumperPrysmoid>) {
            sub.name = "prysmoid_pending";
        } else if constexpr (std::is_same_v<T, BumperQuad>) {
            sub = createPlaneSubMesh(arg.upperPlane, 1.0f, 1.0f);
            sub.name = "quad";
        }
    }, b.bumper);

    return sub;
}

SubMesh BumperGraphMesh::createPrysmoidSubMesh(const BumperPrysmoid &prys, const std::vector<Sphere> &spheres)
{
    SubMesh sub;
    sub.name = "prysmoid";

    int i0 = prys.sphereIndex[0], i1 = prys.sphereIndex[1], i2 = prys.sphereIndex[2];
    if (i0 < 0 || i1 < 0 || i2 < 0 ||
        i0 >= spheres.size() || i1 >= spheres.size() || i2 >= spheres.size())
        return sub;

    glm::vec3 U0 = prys.upperPlane.project(spheres[i0].center);
    glm::vec3 U1 = prys.upperPlane.project(spheres[i1].center);
    glm::vec3 U2 = prys.upperPlane.project(spheres[i2].center);

    glm::vec3 M0 = prys.midPlane.project(spheres[i0].center);
    glm::vec3 M1 = prys.midPlane.project(spheres[i1].center);
    glm::vec3 M2 = prys.midPlane.project(spheres[i2].center);

    int vU0 = sub.vertices.size(); sub.vertices.push_back(U0);
    int vU1 = sub.vertices.size(); sub.vertices.push_back(U1);
    int vU2 = sub.vertices.size(); sub.vertices.push_back(U2);
    sub.faces.push_back({ vU0, vU1, vU2 });

    int vM0 = sub.vertices.size(); sub.vertices.push_back(M0);
    int vM1 = sub.vertices.size(); sub.vertices.push_back(M1);
    int vM2 = sub.vertices.size(); sub.vertices.push_back(M2);

    sub.faces.push_back({ vM0, vM2, vM1 });

    sub.faces.push_back({ vU0, vU1, vM1 });
    sub.faces.push_back({ vU0, vM1, vM0 });

    sub.faces.push_back({ vU1, vU2, vM2 });
    sub.faces.push_back({ vU1, vM2, vM1 });

    sub.faces.push_back({ vU2, vU0, vM0 });
    sub.faces.push_back({ vU2, vM0, vM2 });

    return sub;
}

SubMesh BumperGraphMesh::createExtraLateralCutPlaneForPrysmoidSideWithElongation(
    const BumperPrysmoid &prys,
    const std::vector<Sphere> &spheres,
    int side,
    float topOffset,
    float sideElongation)
{
    SubMesh sub;
    sub.name = "prysmoid_extra_cut_side";

    int i0 = prys.sphereIndex[0], i1 = prys.sphereIndex[1], i2 = prys.sphereIndex[2];
    if (i0 < 0 || i1 < 0 || i2 < 0 ||
        i0 >= spheres.size() || i1 >= spheres.size() || i2 >= spheres.size())
        return sub;

    glm::vec3 U0 = prys.upperPlane.project(spheres[i0].center);
    glm::vec3 U1 = prys.upperPlane.project(spheres[i1].center);
    glm::vec3 U2 = prys.upperPlane.project(spheres[i2].center);

    glm::vec3 M0 = prys.midPlane.project(spheres[i0].center);
    glm::vec3 M1 = prys.midPlane.project(spheres[i1].center);
    glm::vec3 M2 = prys.midPlane.project(spheres[i2].center);

    glm::vec3 topLeft, topRight, bottomLeft, bottomRight;
    if (side == 0) {
        topLeft = U0;   topRight = U1;
        bottomLeft = M0; bottomRight = M1;
    } else if (side == 1) {
        topLeft = U1;   topRight = U2;
        bottomLeft = M1; bottomRight = M2;
    } else if (side == 2) {
        topLeft = U2;   topRight = U0;
        bottomLeft = M2; bottomRight = M0;
    } else
        return sub;

    glm::vec3 lateralNormal = normalize(cross(topRight - topLeft, bottomLeft - topLeft));

    glm::vec3 globalUp(0,1,0);
    if (dot(globalUp, prys.upperPlane.n) <= 0) globalUp = glm::vec3(0,-1,0);

    glm::vec3 T = globalUp - dot(globalUp, lateralNormal) * lateralNormal;

    if (length(T) > 1e-6f)
        T = normalize(T) * topOffset;
    else
        T = glm::vec3(0);

    glm::vec3 topRightVec = normalize(topRight - topLeft);
    glm::vec3 bottomRightVec = normalize(bottomRight - bottomLeft);

    glm::vec3 dispTopLeft = -topRightVec * sideElongation;
    glm::vec3 dispTopRight = topRightVec * sideElongation;
    glm::vec3 dispBottomLeft = -bottomRightVec * sideElongation;
    glm::vec3 dispBottomRight = bottomRightVec * sideElongation;

    glm::vec3 newTopLeft  = topLeft  + T + dispTopLeft;
    glm::vec3 newTopRight = topRight + T + dispTopRight;

    glm::vec3 newBottomLeft  = bottomLeft + dispBottomLeft;
    glm::vec3 newBottomRight = bottomRight + dispBottomRight;

    int v0 = sub.vertices.size(); sub.vertices.push_back(newTopLeft);
    int v1 = sub.vertices.size(); sub.vertices.push_back(newTopRight);
    int v2 = sub.vertices.size(); sub.vertices.push_back(newBottomRight);
    int v3 = sub.vertices.size(); sub.vertices.push_back(newBottomLeft);
    sub.addQuad(v0, v1, v2, v3);

    return sub;
}

BumperGraphMesh::BumperGraphMesh(const BumperGraph &bg) {
    submeshes.clear();
	std::unordered_map<std::pair<int,int>, int, hash_pair> capsuleMap;
    for (size_t i = 0; i < bg.bumper.size(); i++) {
        const Bumper &b = bg.bumper[i];
        std::visit([&]<typename T0>(T0 &&arg) {
            using T = std::decay_t<T0>;
            if constexpr (std::is_same_v<T, BumperCapsuloid>) {
	            const int idx0 = arg.neibExtreme[0] >= 0 ? arg.neibExtreme[0] : arg.sphereIndex[0];
	            const int idx1 = arg.neibExtreme[1] >= 0 ? arg.neibExtreme[1] : arg.sphereIndex[1];
            	
                if (idx0 >= 0 && idx1 >= 0) {
	                const std::pair<int,int> key = std::minmax(idx0, idx1);
                    SubMesh cap = createClosedCapsuleSubMesh(bg.sphere[key.first], bg.sphere[key.second], 128);
                    cap.name = "capsuloid_" + std::to_string(key.first) + "_" + std::to_string(key.second);
                    capsuleMap[key] = submeshes.size();
                    submeshes.push_back(cap);
                }
            }
        }, b.bumper);
    }

    for (size_t i = 0; i < bg.bumper.size(); i++) {
        const Bumper &b = bg.bumper[i];
        std::visit([&]<typename T0>(T0 &&arg) {
            using T = std::decay_t<T0>;
            if constexpr (!std::is_same_v<T, BumperCapsuloid>) {
            	if constexpr (std::is_same_v<T, BumperPrysmoid>) {
		            const SubMesh sub = createPrysmoidSubMesh(arg, bg.sphere);
            		submeshes.push_back(sub);
            		float topOffset = 2.0f;
					float sideElongation = 3.0f;

					for (int side = 0; side < 3; ++side) {
						SubMesh extraSide = createExtraLateralCutPlaneForPrysmoidSideWithElongation(arg,
							bg.sphere, side, topOffset, sideElongation);
						extraSide.name = "prysmoid_extra_side_" + std::to_string(side);
						submeshes.push_back(extraSide);
					}
				} else {
		            const SubMesh sm = createSubMeshFromBumper(b, bg.sphere);
                    submeshes.push_back(sm);
                }
            }
        }, b.bumper);
    }
}

void BumperGraphMesh::exportObj(const std::string &filename) const {
    glm::vec3 globalCentroid(0.0f);
    size_t totalVertices = 0;
    for (const auto &[name, vertices, faces] : submeshes) {
        for (const auto &v : vertices) {
            globalCentroid += v;
            ++totalVertices;
        }
    }

    if (totalVertices > 0) globalCentroid /= static_cast<float>(totalVertices);

    float explosionFactor = 0.0f;

    std::ofstream out(filename);
    if (!out)
        throw std::runtime_error("Failed to open file: " + filename);

    int globalVertexOffset = 0;
    for (const auto &[name, vertices, faces] : submeshes) {
        glm::vec3 subCentroid(0.0f);
        for (const auto &v : vertices)
            subCentroid += v;
        if (!vertices.empty())
            subCentroid /= static_cast<float>(vertices.size());
        glm::vec3 offset = (subCentroid - globalCentroid) * explosionFactor;

        glm::vec3 color(1.0f); // white
        if (name.find("sphere") != std::string::npos)
            color = glm::vec3(1, 0, 0); // red
        else if (name.find("capsuloid") != std::string::npos)
            color = glm::vec3(0, 1, 0); // green
        else if (name.find("prysmoid") != std::string::npos)
            color = glm::vec3(0, 0, 1); // blue

        out << "o " << name << "\n";
        for (const auto &v : vertices) {
            out << "v "
        		<< v.x + offset.x << " "
                << v.y + offset.y << " "
                << v.z + offset.z << " "
                << color.r << " " << color.g << " " << color.b << "\n";
        }

        for (const auto &f : faces) {
            out << "f "
                << f[0] + 1 + globalVertexOffset << " "
                << f[1] + 1 + globalVertexOffset << " "
                << f[2] + 1 + globalVertexOffset << "\n";
        }
        globalVertexOffset += vertices.size();
    }
    out.close();
}

void BumperGraphMesh::exportPly(const std::string &filename) const {
    glm::vec3 globalCentroid(0.0f);
    size_t totalVertices = 0;
    for (const auto &[name, vertices, faces] : submeshes) {
        for (const auto &v : vertices) {
            globalCentroid += v;
            ++totalVertices;
        }
    }

    if (totalVertices > 0) globalCentroid /= static_cast<float>(totalVertices);

    float explosionFactor = 0.0f;

    std::vector<VertexColor> combinedVertices;
    std::vector<std::array<int, 3>> combinedFaces;
    int vertexOffset = 0;

    for (const auto &[name, vertices, faces] : submeshes) {
        glm::vec3 color(1.0f); // white
        if (name.find("sphere") != std::string::npos)
            color = glm::vec3(1, 0, 0); // red
        else if (name.find("capsuloid") != std::string::npos)
            color = glm::vec3(0, 1, 0); // green
        else if (name.find("prysmoid") != std::string::npos)
            color = glm::vec3(0, 0, 1); // blue
        else if (name.find("quad") != std::string::npos)
            color = glm::vec3(128/255.0f, 128/255.0f, 128/255.0f); // gray

        glm::vec3 subCentroid(0.0f);
        for (const auto &v : vertices)
            subCentroid += v;
        if (!vertices.empty())
            subCentroid /= static_cast<float>(vertices.size());
        glm::vec3 offset = (subCentroid - globalCentroid) * explosionFactor;

        for (const auto &v : vertices) {
            VertexColor vc {};
            vc.pos = v + offset;

            vc.r = static_cast<unsigned char>(glm::clamp(color.r, 0.0f, 1.0f) * 255);
            vc.g = static_cast<unsigned char>(glm::clamp(color.g, 0.0f, 1.0f) * 255);
            vc.b = static_cast<unsigned char>(glm::clamp(color.b, 0.0f, 1.0f) * 255);

            combinedVertices.push_back(vc);
        }

        for (const auto &f : faces)
            combinedFaces.push_back({ f[0] + vertexOffset, f[1] + vertexOffset, f[2] + vertexOffset });

        vertexOffset += vertices.size();
    }

    std::ofstream out(filename);
    if (!out)
        throw std::runtime_error("Failed to open file: " + filename);

    out << "ply\n";
    out << "format ascii 1.0\n";
    out << "element vertex " << combinedVertices.size() << "\n";
    out << "property float x\n";
    out << "property float y\n";
    out << "property float z\n";
    out << "property uchar red\n";
    out << "property uchar green\n";
    out << "property uchar blue\n";
    out << "element face " << combinedFaces.size() << "\n";
    out << "property list uchar int vertex_index\n";
    out << "end_header\n";

    for (const auto &[pos, r, g, b] : combinedVertices)
        out << pos.x << " " << pos.y << " " << pos.z << " "
            << static_cast<int>(r) << " " << static_cast<int>(g) << " " << static_cast<int>(b) << "\n";

    for (const auto &f : combinedFaces)
        out << "3 " << f[0] - 1 << " " << f[1] - 1 << " " << f[2] - 1 << "\n";

    out.close();
}