#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <algorithm>

#include <glm/gtc/matrix_transform.hpp>

#include "bumper_graph.h"

using namespace SM;
using namespace SM::Graph;

float BumperGraph::PUSH_EPSILON = 0.000001f;

bool BumperCapsuloid::operator == (const BumperCapsuloid &bc) const
{
	return sphereIndex[0] == bc.sphereIndex[0] && sphereIndex[1] == bc.sphereIndex[1];
}

BumperCapsuloid::BumperCapsuloid (const int indexA, const int indexB, const Sphere &s1, const Sphere &s2)
{
	sphereIndex[0] = indexA < indexB ? indexA : indexB;
	sphereIndex[1] = indexA < indexB ? indexB : indexA;

	const Sphere &sA = indexA < indexB ? s1 : s2;
	const Sphere &sB = indexA < indexB ? s2 : s1;

	const glm::vec3 d = sB.center - sA.center;
	const float dLengthSquared = glm::dot(d, d);
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

void BumperPrysmoid::scale(float s)
{
	upperPlane.scale(s);
	midPlane.scale(s);

	m0.scale(s);
	m1.scale(s);
	m2.scale(s);
}

void BumperPrysmoid::rotate(const glm::mat3 &rot)
{
	upperPlane.rotate(rot);
	midPlane.rotate(rot);

	m1.rotate(rot);
	m2.rotate(rot);
	m0.rotate(rot);
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
	std::vector<int> orderedIndices = {indexA, indexB, indexC};
	std::ranges::sort(orderedIndices);

	sphereIndex[0] = orderedIndices[0];
	sphereIndex[1] = orderedIndices[1];
	sphereIndex[2] = orderedIndices[2];

	upperPlane = computeUpperPlane(v0, v1, v2, direction);

	glm::vec3 u = v1.center - v0.center;
	glm::vec3 v = v2.center - v0.center;

	midPlane = Plane(glm::cross(u, v), v0.center);
	if (direction < 0)
		midPlane.flip();

	m0 = Plane(glm::cross(v1.center - v2.center, upperPlane.n) * static_cast<float>(direction), v1.center);
	m1 = Plane(glm::cross(v2.center - v0.center, upperPlane.n) * static_cast<float>(direction), v2.center);
	m2 = Plane(glm::cross(v0.center - v1.center, upperPlane.n) * static_cast<float>(direction), v0.center);
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

	glm::vec3 n = sign * glm::cross(b - a, c - a);
	n = glm::normalize(n);
	glm::vec3 startN = n;

	for (int i = 0; i < 1000; i++){
		a = sa.center + n * sa.radius;
		b = sb.center + n * sb.radius;
		c = sc.center + n * sc.radius;

		glm::vec3 new_n = glm::normalize(sign * glm::cross(b - a, c - a));
		if (glm::dot(n, new_n) >= 0.999f)
		{
			if (glm::dot(startN, new_n) < 0)
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
	const auto rotationMatrix = glm::mat3(glm::rotate(glm::mat4(1.0f), angleRadians, glm::vec3(0, 1, 0)));

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
	if (bumper[bumperIndex].shapeType == Bumper::SPHERE)	return pushOutsideSphere(p, n, bumperIndex);
	if (bumper[bumperIndex].shapeType == Bumper::CAPSULOID) return pushOutsideCapsuloid(p, n, bumperIndex);
	if (bumper[bumperIndex].shapeType == Bumper::PRYSMOID)	return pushOutsidePrysmoid(p, n, bumperIndex);
	if (bumper[bumperIndex].shapeType == Bumper::QUAD)		return pushOutsideQuad(p, n, bumperIndex);

	return false;
}

void BumperGraph::initializeBumperSpheres()
{
	for (int i = 0; i < sphere.size(); i++)
	{
		BumperSphere bs = BumperSphere();
		bs.sphereIndex = i;

		Bumper node {};
		node.shapeType = Bumper::SPHERE;
		node.bumper = bs;

		bumper.push_back(node);
	}
}

void BumperGraph::initializeBumperCapsuloids(const SphereMesh &sm, std::vector<std::vector<int>> &capsuloidAdj)
{
	for(const Capsuloid &c : sm.capsuloids)
	{
		BumperCapsuloid bc = BumperCapsuloid(
			c.indices[0],
			c.indices[1],
			sphere[c.indices[0]],
			sphere[c.indices[1]]);

		if (capsuloidAdj[c.indices[0]][c.indices[1]] != -1 || capsuloidAdj[c.indices[1]][c.indices[0]] != -1)
			continue;

		bc.neibExtreme[0] = c.indices[0];
		bc.neibExtreme[1] = c.indices[1];
		bc.hasFather = false;

		Bumper node {};
		node.shapeType = Bumper::CAPSULOID;
		node.bumper = bc;

		bumper.push_back(node);

		capsuloidAdj[c.indices[0]][c.indices[1]] = bumper.size() - 1;
		capsuloidAdj[c.indices[1]][c.indices[0]] = bumper.size() - 1;
	}
}

void BumperGraph::initializeBumperPrysmoids(const SphereMesh &sm, std::vector<std::vector<int>> &capsuloidAdj)
{
	for (const Prysmoid &p : sm.prysmoids)
	{
		BumperPrysmoid bpUpper = BumperPrysmoid(
			p.indices[0],
			p.indices[1],
			p.indices[2],
			sphere[p.indices[0]],
			sphere[p.indices[1]],
			sphere[p.indices[2]],
			1
		);

		BumperPrysmoid bpLower = BumperPrysmoid(
			p.indices[0],
			p.indices[1],
			p.indices[2],
			sphere[p.indices[0]],
			sphere[p.indices[1]],
			sphere[p.indices[2]],
			-1
		);

		BumperCapsuloid bc0 = BumperCapsuloid(
			p.indices[0],
			p.indices[1],
			sphere[p.indices[0]],
			sphere[p.indices[1]]);

		BumperCapsuloid bc1 = BumperCapsuloid(
			p.indices[0],
			p.indices[2],
			sphere[p.indices[0]],
			sphere[p.indices[2]]);

		BumperCapsuloid bc2 = BumperCapsuloid(
			p.indices[1],
			p.indices[2],
			sphere[p.indices[1]],
			sphere[p.indices[2]]);

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
		bc0.neibExtreme[0] = p.indices[0];
		bc0.neibExtreme[1] = p.indices[1];
		edge0.bumper = bc0;

		Bumper edge1 {};
		edge1.shapeType = Bumper::CAPSULOID;
		bc1.neibExtreme[0] = p.indices[0];
		bc1.neibExtreme[1] = p.indices[2];
		edge1.bumper = bc1;

		Bumper edge2 {};
		edge2.shapeType = Bumper::CAPSULOID;
		bc2.neibExtreme[0] = p.indices[1];
		bc2.neibExtreme[1] = p.indices[2];
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

		std::get<BumperPrysmoid>(bumper[upperIndex].bumper).neibSide[2] = idx;
		std::get<BumperPrysmoid>(bumper[lowerIndex].bumper).neibSide[2] = idx;

		if (capsuloidAdj[p.indices[0]][p.indices[2]] == -1 || capsuloidAdj[p.indices[2]][p.indices[0]] == -1)
		{
			capsuloidAdj[p.indices[0]][p.indices[2]] = bumper.size();
			capsuloidAdj[p.indices[2]][p.indices[0]] = bumper.size();

			bumper.push_back(edge1);
		}

		idx = capsuloidAdj[p.indices[0]][p.indices[2]];

		std::get<BumperPrysmoid>(bumper[upperIndex].bumper).neibSide[1] = idx;
		std::get<BumperPrysmoid>(bumper[lowerIndex].bumper).neibSide[1] = idx;

		if (capsuloidAdj[p.indices[1]][p.indices[2]] == -1 || capsuloidAdj[p.indices[2]][p.indices[1]] == -1)
		{
			capsuloidAdj[p.indices[1]][p.indices[2]] = bumper.size();
			capsuloidAdj[p.indices[2]][p.indices[1]] = bumper.size();

			bumper.push_back(edge2);
		}

		idx = capsuloidAdj[p.indices[1]][p.indices[2]];

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

		BumperQuad bqUpper = BumperQuad(upper);
		BumperQuad bqLower = BumperQuad(lower);

		BumperCapsuloid bc0 = BumperCapsuloid(
			p.indices[0],
			p.indices[1],
			sphere[p.indices[0]],
			sphere[p.indices[1]]);

		BumperCapsuloid bc1 = BumperCapsuloid(
			p.indices[1],
			p.indices[2],
			sphere[p.indices[1]],
			sphere[p.indices[2]]);

		BumperCapsuloid bc2 = BumperCapsuloid(
			p.indices[2],
			p.indices[3],
			sphere[p.indices[2]],
			sphere[p.indices[3]]);

		BumperCapsuloid bc3 = BumperCapsuloid(
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

		std::get<BumperQuad>(bumper[upperIndex].bumper).neibSide[0] = idx;
		std::get<BumperQuad>(bumper[lowerIndex].bumper).neibSide[0] = idx;

		if (capsuloidAdj[p.indices[1]][p.indices[2]] == -1 || capsuloidAdj[p.indices[2]][p.indices[1]] == -1)
		{
			capsuloidAdj[p.indices[1]][p.indices[2]] = bumper.size();
			capsuloidAdj[p.indices[2]][p.indices[1]] = bumper.size();

			bumper.push_back(edge1);
		}

		idx = capsuloidAdj[p.indices[1]][p.indices[2]];

		std::get<BumperQuad>(bumper[upperIndex].bumper).neibSide[1] = idx;
		std::get<BumperQuad>(bumper[lowerIndex].bumper).neibSide[1] = idx;

		if (capsuloidAdj[p.indices[2]][p.indices[3]] == -1 || capsuloidAdj[p.indices[3]][p.indices[2]] == -1)
		{
			capsuloidAdj[p.indices[2]][p.indices[3]] = bumper.size();
			capsuloidAdj[p.indices[3]][p.indices[2]] = bumper.size();

			bumper.push_back(edge2);
		}

		idx = capsuloidAdj[p.indices[2]][p.indices[3]];

		std::get<BumperQuad>(bumper[upperIndex].bumper).neibSide[2] = idx;
		std::get<BumperQuad>(bumper[lowerIndex].bumper).neibSide[2] = idx;

		if (capsuloidAdj[p.indices[3]][p.indices[0]] == -1 || capsuloidAdj[p.indices[0]][p.indices[3]] == -1)
		{
			capsuloidAdj[p.indices[3]][p.indices[0]] = bumper.size();
			capsuloidAdj[p.indices[0]][p.indices[3]] = bumper.size();

			bumper.push_back(edge3);
		}

		idx = capsuloidAdj[p.indices[3]][p.indices[0]];

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
	glm::vec3 midPlaneNormal = glm::normalize(glm::cross(diagonalAB, diagonalCD));

	Plane midPlane {midPlaneNormal, (s0.center + s1.center + s2.center + s3.center) / 4.0f};

	s0.center = midPlane.project(s0.center);
	s1.center = midPlane.project(s1.center);
	s2.center = midPlane.project(s2.center);
	s3.center = midPlane.project(s3.center);

	if (direction < 0) midPlane.flip();

	Plane upperPlane = commonTangentPlane(s0, s1, s2, s3, direction);

	s0.radius = glm::distance(s0.center, upperPlane.project(s0.center));
	s1.radius = glm::distance(s1.center, upperPlane.project(s1.center));
	s2.radius = glm::distance(s2.center, upperPlane.project(s2.center));
	s3.radius = glm::distance(s3.center, upperPlane.project(s3.center));

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
	std::vector<std::vector<int>> capsuloidAdj (sphere.size(), std::vector<int>(sphere.size(), -1));
	initializeBumperCapsuloids(sm, capsuloidAdj);
	initializeBumperPrysmoids(sm, capsuloidAdj);
	initializeBumperQuads(sm, capsuloidAdj);

	sortByType();
}

void BumperGraph::constructFrom (const SphereMesh &sm)
{
	sphere = sm.spheres;

	initializeBumperNodes(sm);
}

bool BumperGraph::pushOutsideSphere(glm::vec3& p, glm::vec3& n, int& bumperIndex) const
{
	const Bumper& node = bumper[bumperIndex];
	auto [sphereIndex, cones] = std::get<BumperSphere>(node.bumper);
	const Sphere& s = sphere[sphereIndex];

	const glm::vec3 v = p - s.center;
	if (glm::dot(v, v) > s.radius * s.radius)
		return false;

	n = glm::normalize(v);
	p = s.center + n * (s.radius + PUSH_EPSILON);

	return true;
}

float BumperGraph::signedDistanceFromSphere(const int bumperIndex, const glm::vec3& p, glm::vec3& closestPos, glm::vec3& closestNorm)
const
{
	const Bumper& node = bumper[bumperIndex];
	auto [sphereIndex, cones] = std::get<BumperSphere>(node.bumper);
	const Sphere& s = sphere[sphereIndex];

	const float currDist = glm::length(p - s.center);

	closestNorm = (p - s.center) / currDist;
	closestPos = s.center + closestNorm * s.radius;

	return currDist - s.radius;
}

float BumperGraph::closestSphereOn(const glm::vec3& p, const BumperCapsuloid& bc) const
{
	const Sphere& a = sphere[bc.sphereIndex[0]];

	const glm::vec3 pa = p - a.center;
	float t = glm::dot(pa, bc.dNorm);
	const float pq = std::sqrt(glm::dot(pa, pa) - (t * t));

	t +=  bc.k * pq;

	return glm::clamp(t, 0.0f, bc.dLength);
}

Sphere BumperGraph::getInterpolatedSphere(const BumperCapsuloid& bc, float t) const
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
		int index = bc.sphereIndex[0];
		if (!pushOutsideSphere(p, n, index))
			return false;

		bumperIndex = index;
		return true;
	}

	if (t == bc.dLength)
	{
		int index = bc.sphereIndex[1];
		if(!pushOutsideSphere(p, n, index))
			return false;

		bumperIndex = index;
		return true;
	}

	const Sphere s = getInterpolatedSphere(bc, t);
	const glm::vec3 pqPrime = p - s.center;
	const float sqrdDist = glm::dot(pqPrime, pqPrime);

	if (sqrdDist >= s.radius * s.radius)
		return false;

	n = pqPrime / std::sqrt(sqrdDist);
	p = s.center + n * (s.radius + PUSH_EPSILON);

	return true;
}

float BumperGraph::signedDistanceFromCapsuloid(int bumperIndex, const glm::vec3& p, glm::vec3& closestPos, glm::vec3&
closestNorm) const
{
	const Bumper& node = bumper[bumperIndex];
	const BumperCapsuloid bc = std::get<BumperCapsuloid>(node.bumper);

	const float t = closestSphereOn(p, bc);
	const Sphere s = getInterpolatedSphere(bc, t);
	const float dist = glm::distance(p, s.center);

	closestNorm = (p - s.center) / dist;
	closestPos = s.center + closestNorm * s.radius;

	return dist - s.radius;
}

bool BumperGraph::pushOutsidePrysmoid (glm::vec3 &p, glm::vec3 &n, int &bumperIndex) const
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

	const float d = glm::dot(p, glm::normalize(bp.upperPlane.n)) - bp.upperPlane.k;

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

	const float d = glm::dot(p, glm::normalize(bq.upperPlane.n)) - bq.upperPlane.k;

	if (d > 0) return false;

	n = bq.upperPlane.n;
	p -= n * (d - PUSH_EPSILON);

	return true;
}

float BumperGraph::signedDistanceFromPrysmoid(int bumperIndex, const glm::vec3& p, glm::vec3& closestPos, glm::vec3&
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

	const float currDist = glm::length(p - s.center);

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
	const float dist = glm::distance(p, s.center);

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

bool BumperGraph::pushOutsideBruteForce(glm::vec3& p, glm::vec3& n, int& bumperIndex) const
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

bool BumperGraph::isPointOverPrysmoid (const int bumperIndex, const glm::vec3 &p) const
{
	const Bumper& node = bumper[bumperIndex];
	const BumperPrysmoid bp = std::get<BumperPrysmoid>(node.bumper);

	const Sphere& a = sphere[bp.sphereIndex[0]];
	const glm::vec3 pa = p - a.center;

	return glm::dot(bp.upperPlane.n, pa) > 0; // If the point is over the prysmoid
}

bool BumperGraph::isPointOverQuad(int bumperIndex, const glm::vec3 &p) const
{
	const Bumper& node = bumper[bumperIndex];
	const BumperQuad bq = std::get<BumperQuad>(node.bumper);

	const Sphere& a = sphere[bq.sphereIndex[0]];
	const glm::vec3 pa = p - a.center;

	return glm::dot(bq.upperPlane.n, pa) > 0; // If the point is over the quad
}

BumperQuad::BumperQuad(const FourSpheres &fs)
{
	upperPlane = fs.upperPlane;
	midPlane = fs.midPlane;

	sphereIndex[0] = fs.indices[0];
	sphereIndex[1] = fs.indices[1];
	sphereIndex[2] = fs.indices[2];
	sphereIndex[3] = fs.indices[3];

	BumperPrysmoid bp = BumperPrysmoid(sphereIndex[0], sphereIndex[1], sphereIndex[2], fs.spheres[0], fs.spheres[1], fs.spheres[2], 1);
	BumperPrysmoid bp1 = BumperPrysmoid(sphereIndex[0], sphereIndex[2], sphereIndex[3], fs.spheres[0], fs.spheres[2], fs.spheres[3], 1);

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
