#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <algorithm>

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtc/matrix_transform.hpp>

#include "bumper_grid.h"

#include "glm/gtx/string_cast.hpp"

using namespace SM;
using namespace SM::Grid;

float BumperGrid::PUSH_EPSILON = 0.000001f;

bool BumperCapsuloid::operator == (const BumperCapsuloid &bc) const
{
	return sphereIndex[0] == bc.sphereIndex[0] && sphereIndex[1] == bc.sphereIndex[1];
}

BumperCapsuloid::BumperCapsuloid (const int indexA, const int indexB, const Sphere &s1, const Sphere &s2)
{
	sphereIndex[0] = indexA;
	sphereIndex[1] = indexB;

	const Sphere &sA = s1;
	const Sphere &sB = s2;

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
		std::cerr << "Capsuloid " << indexA << ", " << indexB << " degenerate" << std::endl;
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

void BumperPrysmoid::scale(const float s)
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
	std::vector orderedIndices = {indexA, indexB, indexC};
	std::ranges::sort(orderedIndices);

	sphereIndex[0] = orderedIndices[0];
	sphereIndex[1] = orderedIndices[1];
	sphereIndex[2] = orderedIndices[2];

	upperPlane = computeUpperPlane(v0, v1, v2, direction);

	glm::vec3 u = v1.center - v0.center;
	glm::vec3 v = v2.center - v0.center;

	midPlane = Plane(cross(u, v), v0.center);
	if (direction < 0)
		midPlane.flip();

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

std::pair<int, float> BumperGrid::signedDistanceFromBumper(const int i, const glm::vec3& p) const
{
	if (bumper[i].shapeType == Bumper::SPHERE)
		return sampleSignedDistanceFromSphere(i, p);
	if (bumper[i].shapeType == Bumper::CAPSULOID)
		return sampleSignedDistanceFromCapsuloid(i, p);
	if (bumper[i].shapeType == Bumper::PRYSMOID)
		return sampleSignedDistanceFromPrysmoid(i, p);

	return sampleSignedDistanceFromQuad(i, p);
}

void BumperGrid::sampleDistanceFromBumperWithPosition(const glm::vec<3, float> &p, const int proposedIndex, glm::vec3 &closestPos) const
{
	glm::vec3 closestNorm;
	if (bumper[proposedIndex].shapeType == Bumper::SPHERE)
		signedDistanceFromSphere(proposedIndex, p, closestPos, closestNorm);
	else if (bumper[proposedIndex].shapeType == Bumper::CAPSULOID)
		signedDistanceFromCapsuloid(proposedIndex, p, closestPos, closestNorm);
	else if (bumper[proposedIndex].shapeType == Bumper::PRYSMOID)
		signedDistanceFromPrysmoid(proposedIndex, p, closestPos, closestNorm);
	else
		signedDistanceFromQuad(proposedIndex, p, closestPos, closestNorm);
}

void BumperGrid::createGridSamplesFrom(const std::vector<glm::vec3>& pos, const float skinThickness)
{
	gridSamples.clear();
	std::unordered_map<int, int> pushedOutSpheres;
	const float maxDistOffset = skinThickness;

	for (const auto& p : pos)
	{
		int proposedIndex = -1;
		float proposedDist = FLT_MAX;
		glm::vec3 proposedDest = {};
		bool rejected = false;

		for (int i = 0; i < bumper.size(); i++)
		{
			if (bumper[i].hasAParent()) continue;

			auto [idx, dst] = signedDistanceFromBumper(i, p);

			if (dst <= 0.0f && dst >= -maxDistOffset && dst <= proposedDist)
			{
				proposedIndex = idx;
				proposedDist = dst;
				sampleDistanceFromBumperWithPosition(p, proposedIndex, proposedDest);
			}
			else if (dst < -maxDistOffset)
				rejected = true;
		}

		if (proposedIndex == -1 || rejected)
			continue;

		bool inside = false;
		for (int i = 0; i < bumper.size(); i++)
		{
			if (bumper[i].hasAParent()) continue;

			auto [idx, dst] = signedDistanceFromBumper(i, proposedDest);

			if (dst <= 0.0f)
			{
				inside = true;
				break;
			}
		}

		if (!inside)
			continue;

		gridSamples.push_back({proposedDest, proposedDest, proposedIndex});
		gridSamples.push_back({p, proposedDest, proposedIndex});
	}
}

void BumperGrid::translate(const glm::vec3 &t)
{
	for (Sphere &s : sphere)
		s.center += t;

	for (auto& bn : bumper)
		if (bn.shapeType == Bumper::PRYSMOID)
		{
			auto &bp = std::get<BumperPrysmoid>(bn.bumper);
			bp.translate(t);
		}

	grid.translateGrid(t);
}

void BumperGrid::scale(const float s)
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

	grid.scaleGrid(s);
}

void BumperGrid::rotateY(const int angle)
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

	grid.rotateYGrid(angle);
}

bool BumperGrid::serialize(const std::string &filepath) const
{
	std::ofstream file(filepath);

	if (!file.is_open()) {
		std::cerr << "Error: Could not open file for writing: " << filepath << std::endl;
		return false;
	}

	file << "BumperGrid 2\n";
	int capsuloids = 0, prysmoids = 0, quads = 0, composite = 0;
	for (auto& bn : bumper)
	{
		if (bn.shapeType == Bumper::SPHERE) continue;
		if (bn.shapeType == Bumper::CAPSULOID) ++capsuloids;
		else if (bn.shapeType == Bumper::PRYSMOID) ++prysmoids;
		else if (bn.shapeType == Bumper::QUAD) ++quads;
		else ++composite;
	}

	file << sphere.size() << " " << capsuloids << " " << prysmoids << " " << quads << " " << composite << "\n";
	for (auto& s : sphere)
		file << s.center.x << " " << s.center.y << " " << s.center.z << " " << s.radius << "\n";

	for (auto& bn : bumper)
	{
		if (bn.shapeType == Bumper::SPHERE) continue;

		file << bn.serialize() << "\n";
	}

	file << grid.serialize();

	file.close();
	return true;
}

void BumperGrid::deserializeCapsuloids(std::ifstream &file, const int numSpheres, const int numCapsuloids)
{
	file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	for (int i = 0; i < numCapsuloids; ++i) {
		std::string bumperNodeData;
		std::getline(file, bumperNodeData);

		std::istringstream iss(bumperNodeData);
		BumperCapsuloid bc;
		if (!(iss >> bc.sphereIndex[0] >> bc.sphereIndex[1] >>
		      bc.dNorm.x >> bc.dNorm.y >> bc.dNorm.z >>
		      bc.dLength >> bc.dR >> bc.k))
			throw std::runtime_error("Error deserializing BumperCapsuloid data.");

		bumper[numSpheres + i].shapeType = Bumper::CAPSULOID;
		bumper[numSpheres + i].bumper = bc;
	}
}

void BumperGrid::deserializePrysmoids(std::ifstream &file, int numSpheres, int numCapsuloids, int numPrysmoids)
{
	for (int i = 0; i < numPrysmoids; i++)
	{
		std::string bumperNodeData;
		std::getline(file, bumperNodeData);

		std::istringstream iss(bumperNodeData);
		BumperPrysmoid bp;

		if (!(iss >> bp.sphereIndex[0] >> bp.sphereIndex[1] >> bp.sphereIndex[2] >>
		      bp.neibSide[0] >> bp.neibSide[1] >> bp.neibSide[2] >> bp.neibOpp >>
		      bp.m1 >> bp.m2 >> bp.m0 >> bp.midPlane >> bp.upperPlane))
			throw std::runtime_error("Error deserializing BumperPrysmoid data.");

		bumper[numSpheres + numCapsuloids + i].shapeType = Bumper::PRYSMOID;
		bumper[numSpheres + numCapsuloids + i].bumper = bp;
	}
}

void BumperGrid::deserializeQuads(std::ifstream &file, int numSpheres, int numCapsuloids, int numPrysmoids,
	int numQuads)
{
	for (int i = 0; i < numQuads; i++)
	{
		std::string bumperNodeData;
		std::getline(file, bumperNodeData);

		std::istringstream iss(bumperNodeData);
		BumperQuad bq;

		if (!(iss >> bq.sphereIndex[0] >> bq.sphereIndex[1] >> bq.sphereIndex[2] >> bq.sphereIndex[3] >>
			bq.neibSide[0] >> bq.neibSide[1] >> bq.neibSide[2] >> bq.neibSide[3] >>
			bq.neibOpp >>
			bq.sidePlanes[0] >> bq.sidePlanes[1] >> bq.sidePlanes[2] >> bq.sidePlanes[3] >> bq.midPlane >> bq.upperPlane))
			throw std::runtime_error("Error deserializing BumperQuad data.");

		bumper[numSpheres + numCapsuloids + numPrysmoids + i].shapeType = Bumper::QUAD;
		bumper[numSpheres + numCapsuloids + numPrysmoids + i].bumper = bq;
	}
}

void BumperGrid::deserializeComposite(std::ifstream &file, int numSpheres, int numCapsuloids, int numPrysmoids, int numQuads, int numCompositeBumpers)
{
	for (int i = 0; i < numCompositeBumpers; i++)
	{
		std::string bumperNodeData;
		std::getline(file, bumperNodeData);

		std::istringstream iss(bumperNodeData);
		CompositeBumper cb;

		int numBumpers = 0;
		if (!(iss >> numBumpers))
			throw std::runtime_error("Error deserializing CompositeBumper data.");

		for (int j = 0; j < numBumpers; j++)
		{
			std::string bumperIndex;
			if (!(iss >> bumperIndex))
				throw std::runtime_error("Error deserializing CompositeBumper data.");

			if (bumperIndex.find('|') != std::string::npos)
			{
				std::string flag = bumperIndex.substr(bumperIndex.find('|') + 1);
				bumperIndex = bumperIndex.substr(0, bumperIndex.find('|'));

				cb.flags[j] = std::stoi(flag, nullptr, 16);
			}

			cb.bumperIndices[j] = std::stoi(bumperIndex);
		}

		bumper[numSpheres + numCapsuloids + numPrysmoids + numQuads + i].shapeType = Bumper::COMPOSITE;
		bumper[numSpheres + numCapsuloids + numPrysmoids + numQuads + i].bumper = cb;
	}
}

void BumperGrid::deserializeBumpers(std::ifstream& file, const int numSpheres, const int numCapsuloids,
                                    const int numPrysmoids, const int numQuads, const int numCompositeBumpers)
{
	deserializeCapsuloids(file, numSpheres, numCapsuloids);
	deserializePrysmoids(file, numSpheres, numCapsuloids, numPrysmoids);
	deserializeQuads(file, numSpheres, numCapsuloids, numPrysmoids, numQuads);
	deserializeComposite(file, numSpheres, numCapsuloids, numPrysmoids, numQuads, numCompositeBumpers);
}

void BumperGrid::deserialize(const std::string &filepath)
{
	std::ifstream file(filepath);
	if (!file.is_open())
		throw std::runtime_error("Error: Could not open file for reading: " + filepath);

	std::string header;
	int version;
	file >> header >> version;
	if (header != "BumperGrid" || version != 2)
		throw std::runtime_error("Error: Unsupported file format or version in file: " + filepath);

	// Read the number of spheres and bumper nodes
	int numSpheres, numCapsuloids, numPrysmoids, numQuads, numCompositeBumpers;
	file >> numSpheres >> numCapsuloids >> numPrysmoids >> numQuads >> numCompositeBumpers;

	// Resize the vectors to accommodate the data
	sphere.clear();
	bumper.clear();

	sphere.resize(numSpheres);
	bumper.resize(numSpheres + numCapsuloids + numPrysmoids + numQuads + numCompositeBumpers);

	for (int i = 0; i < numSpheres; ++i) {
		float x, y, z, radius;
		file >> x >> y >> z >> radius;
		sphere[i] = Sphere(glm::vec3(x, y, z), radius);

		BumperSphere bs;
		bs.sphereIndex = i;

		bumper[i].shapeType = Bumper::SPHERE;
		bumper[i].bumper = bs;
	}

	deserializeBumpers(file, numSpheres, numCapsuloids, numPrysmoids, numQuads, numCompositeBumpers);

	std::string gridData((std::istreambuf_iterator(file)), std::istreambuf_iterator<char>());
	grid.deserialize(gridData);

	file.close();
}

void BumperGrid::clearAllFlagsForTesting()
{
	for (int& cell : grid.cells)
	{
		if (cell < 0) continue;

		cell &= ~0xF;
		const int idx = cell >> 4;
		if (bumper[idx].shapeType == Bumper::COMPOSITE)
		{
			auto& cb = std::get<CompositeBumper>(bumper[idx].bumper);
			for (short & flag : cb.flags)
				flag = 0;
		}
	}
}

void BumperGrid::populateGrid (WIPGrid& wipGrid) const
{
	wipGrid.populate(gridSamples);
}

std::vector<std::pair<Bumper, int>> BumperGrid::sortCompositeIndices(CompositeBumper &cb)
{
	std::vector<std::pair<Bumper, int>> nodes;
	for (const int & bumperIdx : cb.bumperIndices)
		if (bumperIdx > 0)
			nodes.emplace_back(bumper[bumperIdx], bumperIdx);

	std::ranges::sort(nodes, [](const std::pair<Bumper, int>& a, const std::pair<Bumper, int>& b)
	{
		if (a.first.shapeType != b.first.shapeType)
			return a.first.shapeType < b.first.shapeType;
		return a.second < b.second;
	});

	// reverse nodes
	std::ranges::reverse(nodes);
	return nodes;
}

void BumperGrid::addDontFollowFlagsPrysmoidToPrysmoid(CompositeBumper &cb,
	const std::vector<std::pair<Bumper,int>>& nodes,
	const int i) const
{
	for (int j = i + 1; j < nodes.size(); j++)
	{
		if (nodes[j].first.shapeType != Bumper::PRYSMOID)
			continue;

		const BumperPrysmoid& bp1 = std::get<BumperPrysmoid>(nodes[i].first.bumper);
		const BumperPrysmoid& bp2 = std::get<BumperPrysmoid>(nodes[j].first.bumper);

		auto& [v0_a, v1_a, v2_a] = bp1.sphereIndex;
		auto& [v0_b, v1_b, v2_b] = bp2.sphereIndex;

		if (v1_a == v1_b && v2_a == v2_b)
			cb.flags[i] |= DONT_FOLLOW_SIDE_0;
		if (v0_a == v0_b && v2_a == v2_b)
			cb.flags[i] |= DONT_FOLLOW_SIDE_1;
		if (v0_a == v0_b && v1_a == v1_b)
			cb.flags[i] |= DONT_FOLLOW_SIDE_2;
	}
}

void BumperGrid::addDontFollowFlagsPrysmoidToQuad(
    CompositeBumper &cb,
    const std::vector<std::pair<Bumper, int>>& nodes,
    const int i) const
{
    for (int j = i + 1; j < nodes.size(); j++)
    {
        if (nodes[j].first.shapeType != Bumper::QUAD)
            continue;

        const BumperPrysmoid &bp1 = std::get<BumperPrysmoid>(nodes[i].first.bumper);
        const BumperQuad &bq     = std::get<BumperQuad>(nodes[j].first.bumper);

        auto &[v0_a, v1_a, v2_a] = bp1.sphereIndex;
        auto &[v0_b, v1_b, v2_b, v3_b] = bq.sphereIndex;

        std::array ring = { v0_b, v1_b, v2_b, v3_b };

        auto edgesMatch = [](const int a, const int b, const int c, const int d) {
            return (a == c && b == d) || (a == d && b == c);
        };

        // A helper that checks if a prysm side (p,q) matches ANY of the ring edges
        auto sideMatchesQuadPerimeter = [&](const int p, const int q) {
            // ring: (r[0]->r[1]), (r[1]->r[2]), (r[2]->r[3]), (r[3]->r[0])
            for (int k = 0; k < 4; ++k)
            {
	            const int r0 = ring[k];
	            const int r1 = ring[(k+1) % 4];
                if (edgesMatch(p, q, r0, r1))
                    return true;
            }
            return false;
        };

        if (sideMatchesQuadPerimeter(v1_a, v2_a))
            cb.flags[i] |= DONT_FOLLOW_SIDE_0;
        if (sideMatchesQuadPerimeter(v0_a, v2_a))
            cb.flags[i] |= DONT_FOLLOW_SIDE_1;
        if (sideMatchesQuadPerimeter(v0_a, v1_a))
            cb.flags[i] |= DONT_FOLLOW_SIDE_2;
    }
}

void BumperGrid::addDontFollowFlagsQuadToPrysmoid(
    CompositeBumper &cb,
    const std::vector<std::pair<Bumper, int>> &nodes,
    const int i) const
{
    if (nodes[i].first.shapeType != Bumper::QUAD)
        return;

    const BumperQuad &bq = std::get<BumperQuad>(nodes[i].first.bumper);
    auto &[v0_b, v1_b, v2_b, v3_b] = bq.sphereIndex;

    std::array ring = { v0_b, v1_b, v2_b, v3_b };

    auto edgesMatch = [](const int a, const int b, const int c, const int d) {
        return (a == c && b == d) || (a == d && b == c);
    };

    for (int j = i + 1; j < static_cast<int>(nodes.size()); ++j)
    {
        if (nodes[j].first.shapeType != Bumper::PRYSMOID)
            continue;

        const BumperPrysmoid &bp = std::get<BumperPrysmoid>(nodes[j].first.bumper);
        auto &[v0_a, v1_a, v2_a] = bp.sphereIndex;

        std::array quadSides = {
            std::make_pair(ring[0], ring[1]),
            std::make_pair(ring[1], ring[2]),
            std::make_pair(ring[2], ring[3]),
            std::make_pair(ring[3], ring[0])
        };

        // A helper to check if any of the 3 prysm sides matches a specific Quad side
        auto quadSideMatchesAnyPrysmoidSide = [&](const int qb0, const int qb1) {
            // Compare (qb0, qb1) with each prysm side
            if (edgesMatch(qb0, qb1, v1_a, v2_a)) return true;
            if (edgesMatch(qb0, qb1, v0_a, v2_a)) return true;
            if (edgesMatch(qb0, qb1, v0_a, v1_a)) return true;
            return false;
        };

        if (quadSideMatchesAnyPrysmoidSide(quadSides[0].first, quadSides[0].second))
            cb.flags[i] |= DONT_FOLLOW_SIDE_0;
        if (quadSideMatchesAnyPrysmoidSide(quadSides[1].first, quadSides[1].second))
            cb.flags[i] |= DONT_FOLLOW_SIDE_1;
        if (quadSideMatchesAnyPrysmoidSide(quadSides[2].first, quadSides[2].second))
            cb.flags[i] |= DONT_FOLLOW_SIDE_2;
        if (quadSideMatchesAnyPrysmoidSide(quadSides[3].first, quadSides[3].second))
            cb.flags[i] |= DONT_FOLLOW_SIDE_3;
    }
}

void BumperGrid::addDontFollowFlagsQuadToQuad(
    CompositeBumper &cb,
    const std::vector<std::pair<Bumper, int>> &nodes,
    const int i) const
{
    if (nodes[i].first.shapeType != Bumper::QUAD)
        return;

    const BumperQuad &bq_i = std::get<BumperQuad>(nodes[i].first.bumper);
    auto &[v0_i, v1_i, v2_i, v3_i] = bq_i.sphereIndex;

    std::array ringI = { v0_i, v1_i, v2_i, v3_i };

    std::array quadISides = {
        std::make_pair(ringI[0], ringI[1]),
        std::make_pair(ringI[1], ringI[2]),
        std::make_pair(ringI[2], ringI[3]),
        std::make_pair(ringI[3], ringI[0])
    };

    auto edgesMatch = [](const int a, const int b, const int c, const int d) {
        return (a == c && b == d) || (a == d && b == c);
    };

    for (int j = i + 1; j < static_cast<int>(nodes.size()); ++j)
    {
        if (nodes[j].first.shapeType != Bumper::QUAD)
            continue;

        const BumperQuad &bq_j = std::get<BumperQuad>(nodes[j].first.bumper);
        auto &[v0_j, v1_j, v2_j, v3_j] = bq_j.sphereIndex;

        std::array ringJ = { v0_j, v1_j, v2_j, v3_j };

        std::array quadJSides = {
            std::make_pair(ringJ[0], ringJ[1]),
            std::make_pair(ringJ[1], ringJ[2]),
            std::make_pair(ringJ[2], ringJ[3]),
            std::make_pair(ringJ[3], ringJ[0])
        };

        for (auto &e : quadJSides)
        {
            if (edgesMatch(quadISides[0].first, quadISides[0].second,
                           e.first, e.second))
            {
                cb.flags[i] |= DONT_FOLLOW_SIDE_0;
                break;
            }
        }
        for (auto &e : quadJSides)
        {
            if (edgesMatch(quadISides[1].first, quadISides[1].second,
                           e.first, e.second))
            {
                cb.flags[i] |= DONT_FOLLOW_SIDE_1;
                break;
            }
        }
        for (auto &e : quadJSides)
        {
            if (edgesMatch(quadISides[2].first, quadISides[2].second,
                           e.first, e.second))
            {
                cb.flags[i] |= DONT_FOLLOW_SIDE_2;
                break;
            }
        }
        for (auto &e : quadJSides)
        {
            if (edgesMatch(quadISides[3].first, quadISides[3].second,
                           e.first, e.second))
            {
                cb.flags[i] |= DONT_FOLLOW_SIDE_3;
                break;
            }
        }
    }
}


void BumperGrid::processCompositeBumperFlags(CompositeBumper& cb, const std::set<int>& orig, int cellIdx)
{
	const std::vector<std::pair<Bumper, int>> nodes = sortCompositeIndices(cb);

	for (int i = 0; i < nodes.size(); i++)
	{
		if (nodes[i].first.shapeType == Bumper::PRYSMOID)
		{
			cb.flags[i] = computePrysmoidFlagsFromSet(nodes[i].second, orig);
			// cb.flags[i] = computePrysmoidFlagsFromCell(nodes[i].second, cellIdx);
			addDontFollowFlagsPrysmoidToPrysmoid(cb, nodes, i);
			addDontFollowFlagsPrysmoidToQuad(cb, nodes, i);
		}
		else if (nodes[i].first.shapeType == Bumper::QUAD)
		{
			cb.flags[i] = computeQuadFlagsFromSet(nodes[i].second, orig);
			addDontFollowFlagsQuadToPrysmoid(cb, nodes, i);
			addDontFollowFlagsQuadToQuad(cb, nodes, i);
		}
	}

	for (int i = 0; i < nodes.size(); i++)
		cb.bumperIndices[i] = nodes[i].second;
}

Bumper BumperGrid::makeCompositeBumper(const std::set<int> &orig, const std::set<int> &clean, const int cellIdx)
{
	CompositeBumper cb;

	if (clean.size() > MAX_GRID_BUMPERS)
		std::cout << "Warning: More than " << MAX_GRID_BUMPERS << " bumpers in a cell, baking at most " <<
				MAX_GRID_BUMPERS << "!" << std::endl;

	const int range = std::min(MAX_GRID_BUMPERS, static_cast<int>(clean.size()));
	int i = 0;
	for (const int bumperIdx : clean)
		if (i < range)
			cb.bumperIndices[i++] = bumperIdx;

	processCompositeBumperFlags(cb, orig, cellIdx);

	Bumper node;
	node.shapeType = Bumper::COMPOSITE;
	node.bumper = cb;

	return node;
}

int BumperGrid::makeOrFindCompositeBumper(const std::set<int> &orig, const std::set<int> &clean, const int cellIdx)
{
	const Bumper newBumper = makeCompositeBumper(orig, clean, cellIdx);

	for (int i = 0; i < bumper.size(); i++)
		if (bumper[i] == newBumper)
			return i;

	bumper.push_back(newBumper);
	return static_cast<int>(bumper.size()) - 1;
}

void BumperGrid::computeFinalGrid()
{
	grid.cells.resize(gridWIP.size());
	for (int i = 0; i < gridWIP.size(); i++)
		grid.cells[i] = gridIndexOf(gridWIP[i], i);
}

int BumperGrid::gridIndexOf(const std::set<int> &orig, const int cellIdx)
{
	std::set<int> clean = orig;
	clean = compressChildren(clean);
	clean = compressSiblings(clean);

	int index = 0;
	short flags = 0;
	switch (clean.size())
	{
		case 0:
			return -1;
		case 1:
			index = *clean.begin();
			flags = computePrysmoidFlagsFromSet(index, orig) | computeQuadFlagsFromSet(index, orig);
			// flags = computePrysmoidFlagsFromCell(index, cellIdx);
			break;
		default:
			index = makeOrFindCompositeBumper(orig, clean, cellIdx);
			break;
	}

	return index << 4 | flags;
}

void operator += (std::set<int>& a, const std::set<int>& b)
{
	std::ranges::set_union(a, b, std::inserter(a, a.end()));
}

void operator -= (std::set<int>& a, const std::set<int>& b)
{
	for (const int j : b)
		if (a.contains(j))
			a.erase(j);
}

std::set<int> BumperGrid::compressChildren (const std::set<int>& cell) const
{
	std::set<int> tmp = cell;
	std::set<int> toRemove;
	for (const int j : cell)
		toRemove += descendents(j);
	tmp -= toRemove;
	return tmp;
}

int numberOfIntersections(const std::set<int>& a, const std::set<int>& b)
{
	int count = 0;
	for (int j : a)
		if (b.contains(j))
			count++;
	return count;
}

std::set<int> BumperGrid::compressSiblings (const std::set<int>& cell) const
{
	std::set<int> tmp = cell;
	bool done = false;
	while (!done)
	{
		done = true;
		for (int j = 0; j < bumper.size(); j++)
		{
			std::set<int> d = descendents(j);
			const int intersections = numberOfIntersections(tmp, d);
			if (intersections > 1)
			{
				tmp -= d;
				tmp.insert(j);
				done = false;

				break;
			}
		}
	}

	return tmp;
}

short BumperGrid::computePrysmoidFlagsFromCell(const int element, const int cellIdx) const
{
	constexpr int FLAG_IGNORE_SIDE_0 = 1 << 0;
	constexpr int FLAG_IGNORE_SIDE_1 = 1 << 1;
	constexpr int FLAG_IGNORE_SIDE_2 = 1 << 2;
	constexpr int FLAG_IGNORE_OPP = 1 << 3;

	const std::vector<glm::vec3> cellCoords = grid.getCellCorners(cellIdx);

	if (bumper[element].shapeType != Bumper::PRYSMOID) return 0;

	short flags = 0;
	const auto bp = std::get<BumperPrysmoid>(bumper[element].bumper);
	bp.m0.cellTest(cellCoords) == 1 ? flags |= FLAG_IGNORE_SIDE_0 : 0;
	bp.m1.cellTest(cellCoords) == 1 ? flags |= FLAG_IGNORE_SIDE_1 : 0;
	bp.m2.cellTest(cellCoords) == 1 ? flags |= FLAG_IGNORE_SIDE_2 : 0;
	bp.upperPlane.cellTest(cellCoords) == 1 ? flags |= FLAG_IGNORE_OPP : 0;

	return flags;
}

short BumperGrid::computeQuadFlagsFromSet(const int i, const std::set<int> &set) const
{
	constexpr int FLAG_IGNORE_SIDE_0 = 1 << 0;
	constexpr int FLAG_IGNORE_SIDE_1 = 1 << 1;
	constexpr int FLAG_IGNORE_SIDE_2 = 1 << 2;
	constexpr int FLAG_IGNORE_SIDE_3 = 1 << 7;
	constexpr int FLAG_IGNORE_OPP = 1 << 3;

	const Bumper bn = bumper[i];
	if (bn.shapeType != Bumper::QUAD)
		return 0;

	std::set<int> tmp = set;
	for (const int elem : set)
		if (bumper[elem].shapeType == Bumper::QUAD)
		{
			const BumperQuad bq = std::get<BumperQuad>(bumper[elem].bumper);
			tmp.insert(bq.neibSide[0]);
			tmp.insert(bq.neibSide[1]);
			tmp.insert(bq.neibSide[2]);
			tmp.insert(bq.neibSide[3]);
		}

	const BumperQuad bq = std::get<BumperQuad>(bn.bumper);

	short result = 0;
	if (!set.contains(bq.neibSide[0])) result |= FLAG_IGNORE_SIDE_0;
	if (!set.contains(bq.neibSide[1])) result |= FLAG_IGNORE_SIDE_1;
	if (!set.contains(bq.neibSide[2])) result |= FLAG_IGNORE_SIDE_2;
	if (!set.contains(bq.neibSide[3])) result |= FLAG_IGNORE_SIDE_3;
	if (!set.contains(bq.neibOpp))	   result |= FLAG_IGNORE_OPP;

	return result;
}

short BumperGrid::computeQuadFlagsFromCell(const int element, const int cellIdx) const
{
	constexpr int FLAG_IGNORE_SIDE_0 = 1 << 0;
	constexpr int FLAG_IGNORE_SIDE_1 = 1 << 1;
	constexpr int FLAG_IGNORE_SIDE_2 = 1 << 2;
	constexpr int FLAG_IGNORE_SIDE_3 = 1 << 7;
	constexpr int FLAG_IGNORE_OPP = 1 << 3;

	const std::vector<glm::vec3> cellCoords = grid.getCellCorners(cellIdx);

	if (bumper[element].shapeType != Bumper::QUAD) return 0;

	short flags = 0;
	const auto bq = std::get<BumperQuad>(bumper[element].bumper);
	bq.sidePlanes[0].cellTest(cellCoords) == 1 ? flags |= FLAG_IGNORE_SIDE_0 : 0;
	bq.sidePlanes[1].cellTest(cellCoords) == 1 ? flags |= FLAG_IGNORE_SIDE_1 : 0;
	bq.sidePlanes[2].cellTest(cellCoords) == 1 ? flags |= FLAG_IGNORE_SIDE_2 : 0;
	bq.sidePlanes[3].cellTest(cellCoords) == 1 ? flags |= FLAG_IGNORE_SIDE_3 : 0;
	bq.upperPlane.cellTest(cellCoords) == 1 ? flags |= FLAG_IGNORE_OPP : 0;

	return flags;
}

void BumperGrid::initializeBumperSpheres()
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

void BumperGrid::initializeBumperCapsuloids(const SphereMesh &sm, std::vector<std::vector<int>> &capsuloidAdj)
{
	for(const Capsuloid &c : sm.capsuloids)
	{
		auto bc = BumperCapsuloid(
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

void BumperGrid::initializeBumperPrysmoids(const SphereMesh &sm, std::vector<std::vector<int>> &capsuloidAdj)
{
	for (const Prysmoid &p : sm.prysmoids)
	{
		auto bpUpper = BumperPrysmoid(
			p.indices[0],
			p.indices[1],
			p.indices[2],
			sphere[p.indices[0]],
			sphere[p.indices[1]],
			sphere[p.indices[2]],
			1
		);

		auto bpLower = BumperPrysmoid(
			p.indices[0],
			p.indices[1],
			p.indices[2],
			sphere[p.indices[0]],
			sphere[p.indices[1]],
			sphere[p.indices[2]],
			-1
		);

		auto bc0 = BumperCapsuloid(
			p.indices[0],
			p.indices[1],
			sphere[p.indices[0]],
			sphere[p.indices[1]]);

		auto bc1 = BumperCapsuloid(
			p.indices[0],
			p.indices[2],
			sphere[p.indices[0]],
			sphere[p.indices[2]]);

		auto bc2 = BumperCapsuloid(
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

Quadrilateral BumperGrid::reorderVerticesClockwise(const Quadrilateral& quad) const {
	glm::vec3 v0 = sphere[quad.indices[0]].center;
	glm::vec3 v1 = sphere[quad.indices[1]].center;
	glm::vec3 v2 = sphere[quad.indices[2]].center;
	glm::vec3 v3 = sphere[quad.indices[3]].center;

	glm::vec3 diagonalAB = v0 - v2;
	glm::vec3 diagonalCD = v1 - v3;

	glm::vec3 normal = cross(diagonalAB, diagonalCD);

	// Compute center of the quad
	glm::vec3 center = (v0 + v1 + v2 + v3) / 4.0f;

	// Determine if the vertices are clockwise by comparing normal with center direction
	if (dot(normal, v0 - center) < 0.0f)
		return {quad.indices[0], quad.indices[3], quad.indices[2], quad.indices[1]};

	return quad;
}

FourSpheres BumperGrid::computeFourSpheresFrom(const Quadrilateral& quad, const int direction) const
{
	FourSpheres jq = {sphere[quad.indices[0]], sphere[quad.indices[1]], sphere[quad.indices[2]], sphere[quad.indices[3]]};

	jq.midPlane = jq.computeCommonMidPlane(direction);
	jq.upperPlane = jq.computeCommonTangentPlane(direction);

	jq.indices[0] = quad.indices[0];
	jq.indices[1] = quad.indices[1];
	jq.indices[2] = quad.indices[2];
	jq.indices[3] = quad.indices[3];

	jq.spheres[0] = sphere[quad.indices[0]];
	jq.spheres[1] = sphere[quad.indices[1]];
	jq.spheres[2] = sphere[quad.indices[2]];
	jq.spheres[3] = sphere[quad.indices[3]];

	return jq;
}

void BumperGrid::initializeBumperQuads(const SphereMesh &sm, std::vector<std::vector<int>> &capsuloidAdj)
{
    for (const Quadrilateral& p : sm.quadrilaterals)
    {
        FourSpheres upper = computeFourSpheresFrom(p, 1);
        FourSpheres lower = computeFourSpheresFrom(p, -1);

        auto bqUpper = BumperQuad(upper, 1);
        auto bqLower = BumperQuad(lower, -1);
        std::vector capsuloids = {
            BumperCapsuloid(p.indices[0], p.indices[1], sphere[p.indices[0]], sphere[p.indices[1]]),
            BumperCapsuloid(p.indices[1], p.indices[2], sphere[p.indices[1]], sphere[p.indices[2]]),
            BumperCapsuloid(p.indices[2], p.indices[3], sphere[p.indices[2]], sphere[p.indices[3]]),
            BumperCapsuloid(p.indices[3], p.indices[0], sphere[p.indices[3]], sphere[p.indices[0]])
        };

        Bumper nodeUpper{};
        nodeUpper.shapeType = Bumper::QUAD;
        bqUpper.neibOpp = bumper.size() + 1;
        nodeUpper.bumper = bqUpper;

        Bumper nodeLower{};
        nodeLower.shapeType = Bumper::QUAD;
        bqLower.neibOpp = bumper.size();
        nodeLower.bumper = bqLower;

        bumper.push_back(nodeUpper);
        bumper.push_back(nodeLower);

        int upperIndex = bumper.size() - 2;
        int lowerIndex = bumper.size() - 1;

        // Iterate over each edge to assign capsuloids and neighbors
        for (int edge = 0; edge < 4; ++edge)
        {
            int idx1 = p.indices[edge];
            int idx2 = p.indices[(edge + 1) % 4];

            if (capsuloidAdj[idx1][idx2] == -1 || capsuloidAdj[idx2][idx1] == -1)
            {
                capsuloidAdj[idx1][idx2] = bumper.size();
                capsuloidAdj[idx2][idx1] = bumper.size();

                Bumper edgeBumper{};
                edgeBumper.shapeType = Bumper::CAPSULOID;
                capsuloids[edge].neibExtreme[0] = idx1;
                capsuloids[edge].neibExtreme[1] = idx2;
                edgeBumper.bumper = capsuloids[edge];

                bumper.push_back(edgeBumper);
            }

            int idx = capsuloidAdj[idx1][idx2];

            std::get<BumperCapsuloid>(bumper[idx].bumper).hasFather = true;

            std::get<BumperQuad>(bumper[upperIndex].bumper).neibSide[edge] = idx;
            std::get<BumperQuad>(bumper[lowerIndex].bumper).neibSide[edge] = idx;
        }
    }
}

void BumperGrid::initializeBumperNodes(const SphereMesh &sm)
{
	initializeBumperSpheres();
	std::vector capsuloidAdj (sphere.size(), std::vector(sphere.size(), -1));
	initializeBumperCapsuloids(sm, capsuloidAdj);
	initializeBumperPrysmoids(sm, capsuloidAdj);
	initializeBumperQuads(sm, capsuloidAdj);

	sortByType();
}

void BumperGrid::constructFrom (const SphereMesh &sm)
{
	sphere = sm.spheres;
	bbox = sm.bbox;

	initializeBumperNodes(sm);
}

bool BumperGrid::pushOutsideSphere(glm::vec3& p, glm::vec3& n, const int& bumperIndex) const
{
	const Bumper& node = bumper[bumperIndex];
	auto [sphereIndex] = std::get<BumperSphere>(node.bumper);
	const Sphere& s = sphere[sphereIndex];

	const glm::vec3 v = p - s.center;
	if (dot(v, v) > s.radius * s.radius)
		return false;

	n = normalize(v);
	p = s.center + n * (s.radius + PUSH_EPSILON);

	return true;
}

bool BumperGrid::pushOutsideCapsuloid(glm::vec3 &p, glm::vec3 &n, const int bumperIndex) const
{
	const Bumper& node = bumper[bumperIndex];
	const BumperCapsuloid bc = std::get<BumperCapsuloid>(node.bumper);

	const Sphere& a = sphere[bc.sphereIndex[0]];

	const glm::vec3 pa = p - a.center;
	float t = dot(pa, bc.dNorm);

	t +=  bc.k * std::sqrt(dot(pa, pa) - t * t);
	t = glm::clamp(t, 0.0f, bc.dLength);

	const auto s = Sphere(a.center + t * bc.dNorm, a.radius + t * bc.dR);
	const glm::vec3 pq = p - s.center;
	const float sqrdDist = dot(pq, pq);

	if (sqrdDist >= s.radius * s.radius)
		return false;

	n = pq / std::sqrt(sqrdDist);
	p = s.center + n * (s.radius + PUSH_EPSILON);

	return true;
}

short reverseBits(const short b) {
	short reversed = 0;
	for (int i = 0; i < 8; ++i) {
		if (b & 1 << i) {
			reversed |= 1 << 7 - i;
		}
	}
	return reversed;
}

bool BumperGrid::pushOutsidePrysmoid(glm::vec3 &p, glm::vec3 &n, const int bumperIndex, const short flags) const
{
#ifdef BOOK_KEEP
	static int counter = 0;
	static int totCounter = 0;

	++totCounter;
#endif

	const Bumper& node = bumper[bumperIndex];
	const BumperPrysmoid bp = std::get<BumperPrysmoid>(node.bumper);

	if (!(flags & IGNORE_OPP_PLANE))
	{
		if (bp.midPlane.isBehind(p))
			return pushOutsidePrysmoid(p, n, bp.neibOpp, flags | IGNORE_OPP_PLANE);
	}

	const float m0Distance = bp.m0.distance(p), m1Distance = bp.m1.distance(p), m2Distance = bp.m2.distance(p);
	const float minDistance = std::min({m0Distance, m1Distance, m2Distance});

	if (!(flags & IGNORE_OPP_SIDE_1))
	{
		if (minDistance < 0 && minDistance == m1Distance)
		{
			if (flags & DONT_FOLLOW_SIDE_1)
			{
#ifdef BOOK_KEEP
				counter++;
#endif

				return false;
			}
			return pushOutsideCapsuloid(p, n, bp.neibSide[1]);
		}
	}
	if (!(flags & IGNORE_OPP_SIDE_2))
	{
		if (minDistance < 0 && minDistance == m2Distance)
		{
			if (flags & DONT_FOLLOW_SIDE_2)
			{
#ifdef BOOK_KEEP
				counter++;
#endif
				return false;
			}
			return pushOutsideCapsuloid(p, n, bp.neibSide[2]);
		}
	}
	if (!(flags & IGNORE_OPP_SIDE_0))
	{
		if (minDistance < 0 && minDistance == m0Distance)
		{
			if (flags & DONT_FOLLOW_SIDE_0)
			{
#ifdef BOOK_KEEP
				counter++;
#endif
				return false;
			}
			return pushOutsideCapsuloid(p, n, bp.neibSide[0]);
		}
	}

	const float d = dot(p, bp.upperPlane.n) - bp.upperPlane.k;
	if (d > 0) return false;

	n = bp.upperPlane.n;
	p += n * (-d + PUSH_EPSILON);

	return true;
}

bool BumperGrid::pushOutsideQuad(glm::vec3 &p, glm::vec3 &n, const int bumperIndex, const short flags) const
{
	const Bumper& node = bumper[bumperIndex];
	const BumperQuad bq = std::get<BumperQuad>(node.bumper);

	if (!(flags & IGNORE_OPP_PLANE))
	{
		if (bq.midPlane.isBehind(p))
			return pushOutsideQuad(p, n, bq.neibOpp, flags | IGNORE_OPP_PLANE);
	}

	const float m0Distance = bq.sidePlanes[0].distance(p), m1Distance = bq.sidePlanes[1].distance(p),
	m2Distance = bq.sidePlanes[2].distance(p), m3Distance = bq.sidePlanes[3].distance(p);
	const float minDistance = std::min({m0Distance, m1Distance, m2Distance, m3Distance});

	if (!(flags & IGNORE_OPP_SIDE_0))
	{
		if (minDistance < 0 && minDistance == m0Distance)
		{
			if (flags & DONT_FOLLOW_SIDE_0)
			{
#ifdef BOOK_KEEP
				counter++;
#endif

				return false;
			}
			return pushOutsideCapsuloid(p, n, bq.neibSide[0]);
		}
	}
	if (!(flags & IGNORE_OPP_SIDE_1))
	{
		if (minDistance < 0 && minDistance == m1Distance)
		{
			if (flags & DONT_FOLLOW_SIDE_1)
			{
#ifdef BOOK_KEEP
				counter++;
#endif
				return false;
			}
			return pushOutsideCapsuloid(p, n, bq.neibSide[1]);
		}
	}
	if (!(flags & IGNORE_OPP_SIDE_2))
	{
		if (minDistance < 0 && minDistance == m2Distance)
		{
			if (flags & DONT_FOLLOW_SIDE_2)
			{
#ifdef BOOK_KEEP
				counter++;
#endif
				return false;
			}
			return pushOutsideCapsuloid(p, n, bq.neibSide[2]);
		}
	}
	if (!(flags & IGNORE_OPP_SIDE_3))
	{
		if (minDistance < 0 && minDistance == m3Distance)
		{
			if (flags & DONT_FOLLOW_SIDE_3)
			{
#ifdef BOOK_KEEP
				counter++;
#endif
				return false;
			}
			return pushOutsideCapsuloid(p, n, bq.neibSide[3]);
		}
	}

	const float d = dot(p, bq.upperPlane.n) - bq.upperPlane.k;
	if (d > 0) return false;

	n = bq.upperPlane.n;
	p += n * (-d + PUSH_EPSILON);

	return true;
}

float BumperGrid::signedDistanceFromSphere(const int bumperIndex, const glm::vec3& p, glm::vec3& closestPos, glm::vec3& closestNorm)
const
{
	const Bumper& node = bumper[bumperIndex];
	auto [sphereIndex] = std::get<BumperSphere>(node.bumper);
	const Sphere& s = sphere[sphereIndex];

	const float currDist = length(p - s.center);

	closestNorm = (p - s.center) / currDist;
	closestPos = s.center + closestNorm * s.radius;

	return currDist - s.radius;
}

float BumperGrid::closestSphereOn(const glm::vec3& p, const BumperCapsuloid& bc) const
{
	const Sphere& a = sphere[bc.sphereIndex[0]];

	const glm::vec3 pa = p - a.center;
	float t = dot(pa, bc.dNorm);
	const float pq = std::sqrt(dot(pa, pa) - t * t);

	t +=  bc.k * pq;

	return glm::clamp(t, 0.0f, bc.dLength);
}

Sphere BumperGrid::getInterpolatedSphere(const BumperCapsuloid& bc, const float t) const
{
	Sphere s{};
	const Sphere& a = sphere[bc.sphereIndex[0]];

	s.center = a.center + t * bc.dNorm;
	s.radius = a.radius + t * bc.dR;

	return s;
}

float BumperGrid::signedDistanceFromCapsuloid(const int bumperIndex, const glm::vec3& p, glm::vec3& closestPos, glm::vec3&
                                              closestNorm) const
{
	const Bumper& node = bumper[bumperIndex];
	const BumperCapsuloid bc = std::get<BumperCapsuloid>(node.bumper);

	const float t = closestSphereOn(p, bc);

	if (t == 0)
		return signedDistanceFromSphere(bc.neibExtreme[0], p, closestPos, closestNorm);
	if (t == bc.dLength)
		return signedDistanceFromSphere(bc.neibExtreme[1], p, closestPos, closestNorm);

	const Sphere s = getInterpolatedSphere(bc, t);
	const float dist = distance(p, s.center);

	closestNorm = (p - s.center) / dist;
	closestPos = s.center + closestNorm * s.radius;

	return dist - s.radius;
}

float BumperGrid::signedDistanceFromPrysmoid(const int bumperIndex, const glm::vec3& p, glm::vec3& closestPos, glm::vec3&
                                             closestNorm) const
{
	const Bumper& node = bumper[bumperIndex];
	const BumperPrysmoid bp = std::get<BumperPrysmoid>(node.bumper);

	if (bp.midPlane.isBehind(p)) return signedDistanceFromPrysmoid(bp.neibOpp, p, closestPos, closestNorm);

	const float m0Distance = bp.m0.distance(p), m1Distance = bp.m1.distance(p), m2Distance = bp.m2.distance(p);
	const float minDistance = std::min({m0Distance, m1Distance, m2Distance});

	if (minDistance < 0 && minDistance == m0Distance) return signedDistanceFromCapsuloid(bp.neibSide[0], p, closestPos,  closestNorm);
	if (minDistance < 0 && minDistance == m1Distance) return signedDistanceFromCapsuloid(bp.neibSide[1], p, closestPos,  closestNorm);
	if (minDistance < 0 && minDistance == m2Distance) return signedDistanceFromCapsuloid(bp.neibSide[2], p, closestPos,  closestNorm);

	closestNorm = bp.upperPlane.n;
	closestPos = bp.upperPlane.project(p);

	return bp.upperPlane.distance(p);
}

float BumperGrid::signedDistanceFromQuad(const int bumperIndex, const glm::vec3 &p, glm::vec3 &closestPos,
                                         glm::vec3 &closestNorm) const
{
	const Bumper& node = bumper[bumperIndex];
	const BumperQuad bq = std::get<BumperQuad>(node.bumper);

	if (bq.midPlane.isBehind(p)) return signedDistanceFromQuad(bq.neibOpp, p, closestPos, closestNorm);

	const float m0Distance = bq.sidePlanes[0].distance(p), m1Distance = bq.sidePlanes[1].distance(p),
	m2Distance = bq.sidePlanes[2].distance(p), m3Distance = bq.sidePlanes[3].distance(p);
	const float minDistance = std::min({m0Distance, m1Distance, m2Distance, m3Distance});

	if (minDistance < 0 && minDistance == m0Distance) return signedDistanceFromCapsuloid(bq.neibSide[0], p, closestPos,  closestNorm);
	if (minDistance < 0 && minDistance == m1Distance) return signedDistanceFromCapsuloid(bq.neibSide[1], p, closestPos,  closestNorm);
	if (minDistance < 0 && minDistance == m2Distance) return signedDistanceFromCapsuloid(bq.neibSide[2], p, closestPos,  closestNorm);
	if (minDistance < 0 && minDistance == m3Distance) return signedDistanceFromCapsuloid(bq.neibSide[3], p, closestPos,  closestNorm);

	closestNorm = bq.upperPlane.n;
	closestPos = bq.upperPlane.project(p);

	return bq.upperPlane.distance(p);
}

std::pair<int, float> BumperGrid::sampleSignedDistanceFromSphere (int bumperIndex, const glm::vec3 &p) const
{
	const Bumper& node = bumper[bumperIndex];
	auto [sphereIndex] = std::get<BumperSphere>(node.bumper);
	const Sphere& s = sphere[sphereIndex];

	const float currDist = length(p - s.center);

	return {bumperIndex, currDist - s.radius};
}

std::pair<int, float> BumperGrid::sampleSignedDistanceFromCapsuloid (int bumperIndex, const glm::vec3 &p) const
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

std::pair<int, float> BumperGrid::sampleSignedDistanceFromPrysmoid (int bumperIndex, const glm::vec3 &p) const
{
	const Bumper& node = bumper[bumperIndex];
	const BumperPrysmoid bp = std::get<BumperPrysmoid>(node.bumper);

	if (bp.midPlane.isBehind(p)) return {-1, FLT_MAX};

	const float m0Distance = bp.m0.distance(p), m1Distance = bp.m1.distance(p), m2Distance = bp.m2.distance(p);
	const float minDistance = std::min({m0Distance, m1Distance, m2Distance});

	if (minDistance < 0 && minDistance == m0Distance) return sampleSignedDistanceFromCapsuloid(bp.neibSide[0], p);
	if (minDistance < 0 && minDistance == m1Distance) return sampleSignedDistanceFromCapsuloid(bp.neibSide[1], p);
	if (minDistance < 0 && minDistance == m2Distance) return sampleSignedDistanceFromCapsuloid(bp.neibSide[2], p);

	return {bumperIndex, bp.upperPlane.distance(p)};
}

std::pair<int, float> BumperGrid::sampleSignedDistanceFromQuad(int bumperIndex, const glm::vec3 &p) const
{
	const Bumper& node = bumper[bumperIndex];
	const BumperQuad bq = std::get<BumperQuad>(node.bumper);

	if (bq.midPlane.isBehind(p)) return {-1, FLT_MAX};

	const float m0Distance = bq.sidePlanes[0].distance(p), m1Distance = bq.sidePlanes[1].distance(p),
	m2Distance = bq.sidePlanes[2].distance(p), m3Distance = bq.sidePlanes[3].distance(p);
	const float minDistance = std::min({m0Distance, m1Distance, m2Distance, m3Distance});

	if (minDistance < 0 && minDistance == m0Distance) return sampleSignedDistanceFromCapsuloid(bq.neibSide[0], p);
	if (minDistance < 0 && minDistance == m1Distance) return sampleSignedDistanceFromCapsuloid(bq.neibSide[1], p);
	if (minDistance < 0 && minDistance == m2Distance) return sampleSignedDistanceFromCapsuloid(bq.neibSide[2], p);
	if (minDistance < 0 && minDistance == m3Distance) return sampleSignedDistanceFromCapsuloid(bq.neibSide[3], p);

	return {bumperIndex, bq.upperPlane.distance(p)};
}

float BumperGrid::projectOnBruteForce(glm::vec3& p, glm::vec3& n, int& bumperIndex)  const
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
				default:
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

bool BumperGrid::pushOutsideBruteForce(glm::vec3& p, glm::vec3& n, int& bumperIndex) const
{
	bool hasCollided = false;

	for (int i = 0; i < bumper.size(); i++)
	{
		int index = i;

		if (bumper[i].shapeType == Bumper::SPHERE)
		{
			if (pushOutsideSphere(p, n, index))
			{
				bumperIndex = index;
				hasCollided = true;
			}
		}
		else if (bumper[i].shapeType == Bumper::CAPSULOID)
		{
			if (pushOutsideCapsuloid(p, n, index))
			{
				bumperIndex = index;
				hasCollided = true;
			}
		}
		else if (bumper[i].shapeType == Bumper::PRYSMOID)
		{
			if (pushOutsidePrysmoid(p, n, index))
			{
				bumperIndex = index;
				hasCollided = true;
			}
		}
		else if (bumper[i].shapeType == Bumper::QUAD)
		{
			if (pushOutsideQuad(p, n, index))
			{
				bumperIndex = index;
				hasCollided = true;
			}
		}
	}

	return hasCollided;
}

inline int getRandomInVector(const std::vector<int> & v, const int randSeed){
	if (v.empty()) return -1;
	return v[ randSeed % v.size() ];
}

bool BumperGrid::projectOn (glm::vec3 &p, glm::vec3 &n, int &bumperIndex) const
{
	projectOnBruteForce(p, n, bumperIndex);
	return true;
}

bool BumperGrid::isPointOverPrysmoid (const int bumperIndex, const glm::vec3 &p) const
{
	const Bumper& node = bumper[bumperIndex];
	const BumperPrysmoid bp = std::get<BumperPrysmoid>(node.bumper);

	const Sphere& a = sphere[bp.sphereIndex[0]];
	const glm::vec3 pa = p - a.center;

	return dot(bp.upperPlane.n, pa) > 0; // If the point is over the prysmoid
}

bool BumperGrid::isPointOverQuad(const int bumperIndex, const glm::vec3 &p) const
{
	const Bumper& node = bumper[bumperIndex];
	const BumperQuad bq = std::get<BumperQuad>(node.bumper);

	const Sphere& a = sphere[bq.sphereIndex[0]];
	const glm::vec3 pa = p - a.center;

	return dot(bq.upperPlane.n, pa) > 0; // If the point is over the quad
}

bool BumperGrid::pushOutsideComposite(glm::vec3& p, glm::vec3 n, const int bumperIndex) const
{
	const Bumper& node = bumper[bumperIndex];
	const CompositeBumper& cb = std::get<CompositeBumper>(node.bumper);

	for (int i = 0; i < MAX_GRID_BUMPERS; i++)
	{
		int index = cb.bumperIndices[i];
		if (index == -1)
			break;

		const auto bumperType = bumper[index].shapeType;
		if (bumperType == Bumper::SPHERE)
			pushOutsideSphere(p, n, index);
		if (bumperType == Bumper::CAPSULOID)
			pushOutsideCapsuloid(p, n, index);
		if (bumperType == Bumper::PRYSMOID)
			pushOutsidePrysmoid(p, n, index, cb.flags[i]);
	}

	return false;
}

bool BumperGrid::pushOutsideGrid (glm::vec3 &p, glm::vec3 &n) const
{
	int bumperAtPos = grid.indicesAt(p);
	if (bumperAtPos < 0)
		return false;

	constexpr int FLAGS_MASK = 0xF;
	short flags = bumperAtPos & FLAGS_MASK;
	bumperAtPos >>= 4;

	if (bumper[bumperAtPos].shapeType == Bumper::SPHERE)
	{
		if (pushOutsideSphere(p, n, bumperAtPos)) return true;

		return false;
	}

	if (bumper[bumperAtPos].shapeType == Bumper::CAPSULOID)
		return pushOutsideCapsuloid(p, n, bumperAtPos);
	if (bumper[bumperAtPos].shapeType == Bumper::PRYSMOID)
		return pushOutsidePrysmoid(p, n, bumperAtPos, flags);
	if (bumper[bumperAtPos].shapeType == Bumper::QUAD)
		return pushOutsideQuad(p, n, bumperAtPos);
	if (bumper[bumperAtPos].shapeType == Bumper::COMPOSITE)
		return pushOutsideComposite(p, n, bumperAtPos);

	return true;
}

void BumperGrid::constructGrid (const int numberOfCells)
{
	WIPGrid wipGrid;
	wipGrid.format(bbox, numberOfCells);
	populateGrid(wipGrid);

	wipGrid.bake(*this);
	computeFinalGrid();
}

int SpatialGrid::indicesAt(const glm::vec3& position) const
{
	if (!isIndexValid(getCellCoords(position)))
		return -1;

	return cells[getIndex(position)];
}

void WIPGrid::populate(const std::vector<GridSample>& samples)
{
	cellsSet.resize(numCells());
	for (const auto& [position, dest, index] : samples)
	{
		const int idx = getIndex(position);
		cellsSet[idx].insert(index);
	}

	for (int i = 0; i < cells.size(); i++)
	{
		const std::vector c(cellsSet[i].begin(), cellsSet[i].end());
		cells[i] = c;
	}
}

void WIPGrid::bake(BumperGrid& bg) const
{
	bg.grid.cells.resize(cells.size());
	bg.grid.resolution = resolution;
	bg.grid.offset = offset;
	bg.grid.scale = scale;
	bg.gridWIP = cellsSet;
}

std::string SpatialGrid::serialize() const
{
	std::ostringstream oss;

	oss << resolution.x << " " << resolution.y << " " << resolution.z << "\n";
	oss << offset.x << " " << offset.y << " " << offset.z << "\n";
	oss << scale << "\n";

	const int totalCells = cells.size();

	int lastIdx = cells[0];
	int count = 1;
	for (int idx = 1; idx < totalCells; ++idx)
	{
		const auto& cell = cells[idx];
		if (cell == lastIdx)
			count++;
		else
		{
			if (lastIdx == -1)
			{
				oss << count << " 0 ";
				lastIdx = cell;
				count = 1;

				continue;
			}

			const int realIdx = lastIdx >> 4;
			const short flags = lastIdx & 0xF;

			std::ostringstream hexStream;
			hexStream << std::hex << std::uppercase << std::setw(1) << std::setfill('0') << static_cast<int>(flags);
			std::string hexStr = hexStream.str();

			if (flags > 0)
				oss << count << " " << realIdx + 1 << "|" << hexStr << " ";
			else
				oss << count << " " << realIdx + 1 << " ";

			lastIdx = cell;
			count = 1;
		}
	}

	if (lastIdx == -1)
		oss << count << " 0 ";
	else
	{
		const int realIdx = lastIdx >> 4;
		const short flags = lastIdx & 0xF;

		std::ostringstream hexStream;
		hexStream << std::hex << std::uppercase << std::setw(1) << std::setfill('0') << static_cast<int>(flags);
		std::string hexStr = hexStream.str();

		if (flags > 0)
			oss << count << " " << realIdx + 1 << "|" << hexStr << " ";
		else
			oss << count << " " << realIdx + 1 << " ";
	}

	return oss.str();
}


void SpatialGrid::deserialize(const std::string& data)
{
	std::istringstream iss(data);

	if (!(iss >> resolution.x >> resolution.y >> resolution.z))
		throw std::runtime_error("Error reading resolution values.");

	if (!(iss >> offset.x >> offset.y >> offset.z))
		throw std::runtime_error("Error reading offset values.");

	if (!(iss >> scale))
		throw std::runtime_error("Error reading scale value.");

	const int size = resolution.x * resolution.y * resolution.z;
	cells.clear();
	cells.resize(size);

	int count;
	std::string idx;

	int position = 0;
	while (position < size)
	{
		if (!(iss >> count >> idx))
			throw std::runtime_error("Error reading cell data.");

		for (int j = 0; j < count; ++j)
		{
			if (position >= size)
				throw std::runtime_error("Deserialization error: Exceeded expected cell count.");

			const size_t pos = idx.find('|');
			if (pos != std::string::npos)
			{
				const int realIdx = std::stoi(idx.substr(0, pos)) - 1;
				const short flags = std::stoi(idx.substr(pos + 1), nullptr, 16);

				cells[position++] = realIdx << 4 | flags;
			}
			else
				cells[position++] = std::stoi(idx) - 1 << 4;
		}
	}
}

void WIPGrid::format(const AABB& boundingBox, const int requiredNumCells)
{
    const float bdd = length(boundingBox.maxCorner - boundingBox.minCorner);

    float minScale = bdd / requiredNumCells;
    float maxScale = bdd / (requiredNumCells / 2.0f);

    int minVal = howManyCellsForScale(boundingBox, minScale);
    int maxVal = howManyCellsForScale(boundingBox, maxScale);

    while (minVal > requiredNumCells) {
        minScale *= 0.5f;
        minVal = howManyCellsForScale(boundingBox, minScale);
    }

    while (maxVal < requiredNumCells) {
        maxScale *= 2.0f;
        maxVal = howManyCellsForScale(boundingBox, maxScale);
    }

    // Perform binary search to find the optimal scale
    for (int i = 0; i < 100; i++)
    {
        const float midScale = (minScale + maxScale) / 2.0f;
        const int midVal = howManyCellsForScale(boundingBox, midScale);

        if (midVal < requiredNumCells)
        {
            minScale = midScale;
            minVal = midVal;
        }
        else
        {
            maxScale = midScale;
            maxVal = midVal;
        }
    }

    // Choose the closest value to the required number of cells
    if (std::abs(minVal - requiredNumCells) < std::abs(maxVal - requiredNumCells))
        setSize(boundingBox, minScale);
    else
        setSize(boundingBox, maxScale);

#ifdef DEBUG
    std::cout << "Total cells: " << numCells() << " (Requested: " << requiredNumCells << ")" << std::endl;
    std::cout << "Resolution: " << resolution.x << "x" << resolution.y << "x" << resolution.z << std::endl;
    cells.resize(resolution.x * resolution.y * resolution.z);
#endif
}

void WIPGrid::setSize(const AABB& boundingBox, const float requiredScale)
{
	offset = {0.0f, 0.0f, 0.0f};
	scale = requiredScale;

	const glm::ivec3 minCell = getCellCoords(boundingBox.minCorner);
	const glm::ivec3 maxCell = getCellCoords(boundingBox.maxCorner);

	offset = {-static_cast<float>(minCell.x), -static_cast<float>(minCell.y), -static_cast<float>(minCell.z)};
	resolution = maxCell - minCell + glm::ivec3(1);
}

int WIPGrid::numCells () const
{
	return resolution.x * resolution.y * resolution.z;
}

bool SpatialGrid::isIndexValid (const glm::ivec3 &p) const
{
	return p.x >= 0 && p.x < resolution.x &&
	       p.y >= 0 && p.y < resolution.y &&
	       p.z >= 0 && p.z < resolution.z;
}

int SpatialGrid::getFlatIndex(const int i, const int j, const int k) const
{
	return i + resolution.x * j + k * resolution.y * resolution.x;
}

int SpatialGrid::getIndex (const glm::vec3 &p) const
{
	const glm::vec3 i = getCellCoords(p);
	return i.x + resolution.x * i.y + i.z * resolution.y * resolution.x;
}

int WIPGrid::getIndex (const glm::vec3 &p) const
{
	const glm::vec3 i = getCellCoords(p);
	return i.x + resolution.x * i.y + i.z * resolution.y * resolution.x;
}

glm::ivec3 SpatialGrid::getCellCoords(const glm::vec3& p) const
{
	return floor(p * scale + offset);
}

std::vector<glm::vec3> SpatialGrid::getCellCorners(const int i) const {
	std::vector<glm::vec3> corners;

	const auto cellCoords = glm::ivec3(
		i % resolution.x,
		i / resolution.x % resolution.y,
		i / (resolution.x * resolution.y)
	);

	// Converting to world space
	const glm::vec3 basePosition = (glm::vec3(cellCoords) - offset) / scale;

	auto cellSize = glm::vec3(1.0f / scale);
	const glm::vec3 offsets[8] = {
		{0, 0, 0},
		{cellSize.x, 0, 0},
		{cellSize.x, cellSize.y, 0},
		{0, cellSize.y, 0},
		{0, 0, cellSize.z},
		{cellSize.x, 0, cellSize.z},
		{cellSize.x, cellSize.y, cellSize.z},
		{0, cellSize.y, cellSize.z}
	};

	for (auto offset : offsets)
		corners.push_back(basePosition + offset);

	return corners;
}

void SpatialGrid::translateGrid(const glm::vec3 &t)
{
	const glm::vec3 adjustedTranslation = t * scale;
	offset -= adjustedTranslation;
}

void SpatialGrid::scaleGrid(const float s)
{
	const float newScale = scale * (1 / s);
	scale = newScale;
}

void SpatialGrid::rotateYGrid(const int angle)
{
    int rotation = static_cast<int>(angle) % 360;
    if (rotation < 0)
        rotation += 360;

    auto temp = SpatialGrid(*this);

    glm::ivec3 newResolution = resolution;
    glm::vec3 newOffset = offset;

    switch (rotation) {
        case 90:
            newResolution = { resolution.z, resolution.y, resolution.x };
            newOffset = { offset.z, offset.y, resolution.x - offset.x };
            break;

        case 180:
            newOffset = { -offset.x + resolution.x, offset.y, -offset.z + resolution.z };
            break;

        case 270:
            newResolution = { resolution.z, resolution.y, resolution.x };
            newOffset = { resolution.z - offset.z, offset.y, offset.x };
            break;

        case 0:
        default:
            return;
    }

    resolution = newResolution;
    offset = newOffset;

    std::vector newCells(cells.size(), 0);

    for (int x = 0; x < temp.resolution.x; ++x)
        for (int y = 0; y < temp.resolution.y; ++y)
            for (int z = 0; z < temp.resolution.z; ++z) {
                glm::ivec3 originalCoord = { x, y, z };
                glm::ivec3 newCoord;

                switch (rotation) {
                    case 90:
                        newCoord = { z, y, temp.resolution.x - 1 - x };
                        break;

                    case 180:
                        newCoord = { temp.resolution.x - 1 - x, y, temp.resolution.z - 1 - z };
                        break;

                    case 270:
                        newCoord = { temp.resolution.z - 1 - z, y, x };
                        break;
                	default:
						std::cerr << "Invalid rotation angle: " << rotation << std::endl;
						return;
                }

                int oldIndex = temp.getFlatIndex(originalCoord.x, originalCoord.y, originalCoord.z);
                int newIndex = getFlatIndex(newCoord.x, newCoord.y, newCoord.z);

                newCells[newIndex] = temp.cells[oldIndex];
            }

    cells = std::move(newCells);
}

glm::ivec3 WIPGrid::getCellCoords(const glm::vec3& p) const
{
	return floor(p * scale + offset);
}

int WIPGrid::howManyCellsForScale (const AABB &bbox, const float scale)
{
	WIPGrid grid;
	grid.setSize(bbox, scale);
	return grid.numCells();
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
	else if (shapeType == QUAD)
	{
		auto bq = std::get<BumperQuad>(bumper);
		oss << bq.sphereIndex[0]
		<< " " << bq.sphereIndex[1]
		<< " " << bq.sphereIndex[2]
		<< " " << bq.sphereIndex[3]
		<< " " << bq.neibSide[0]
		<< " " << bq.neibSide[1]
		<< " " << bq.neibSide[2]
		<< " " << bq.neibSide[3]
		<< " " << bq.neibOpp
		<< " " << bq.sidePlanes[0].serialize()
		<< " " << bq.sidePlanes[1].serialize()
		<< " " << bq.sidePlanes[2].serialize()
		<< " " << bq.sidePlanes[3].serialize()
		<< " " << bq.midPlane.serialize()
		<< " " << bq.upperPlane.serialize();
	}
	else if (shapeType == COMPOSITE)
	{
		auto cb = std::get<CompositeBumper>(bumper);
		oss << cb.serialize();
	}
	else
		throw std::runtime_error("Unknown shape type encountered during serialization.");

	return oss.str();
}

short BumperGrid::computePrysmoidFlagsFromSet(int i, const std::set<int> &set) const
{
	constexpr int FLAG_IGNORE_SIDE_0 = 1 << 0;
	constexpr int FLAG_IGNORE_SIDE_1 = 1 << 1;
	constexpr int FLAG_IGNORE_SIDE_2 = 1 << 2;
	constexpr int FLAG_IGNORE_OPP = 1 << 3;

	const Bumper bn = bumper[i];
	if (bn.shapeType != Bumper::PRYSMOID)
		return 0;

	std::set<int> tmp = set;
	for (const int elem : set)
		if (bumper[elem].shapeType == Bumper::PRYSMOID)
		{
			const BumperPrysmoid bp = std::get<BumperPrysmoid>(bumper[elem].bumper);
			tmp.insert(bp.neibSide[0]);
			tmp.insert(bp.neibSide[1]);
			tmp.insert(bp.neibSide[2]);
		}

	const BumperPrysmoid bp = std::get<BumperPrysmoid>(bn.bumper);

	short result = 0;
	if (!set.contains(bp.neibSide[0])) result |= FLAG_IGNORE_SIDE_0;
	if (!set.contains(bp.neibSide[1])) result |= FLAG_IGNORE_SIDE_1;
	if (!set.contains(bp.neibSide[2])) result |= FLAG_IGNORE_SIDE_2;
	if (!set.contains(bp.neibOpp))	   result |= FLAG_IGNORE_OPP;

	return result;
}

bool Bumper::operator == (const Bumper& bn) const
{
	if (shapeType != bn.shapeType) return false;

	switch (shapeType)
	{
		case SPHERE: return std::get<BumperSphere>(bumper).sphereIndex == std::get<BumperSphere>(bn.bumper).sphereIndex;
		case CAPSULOID: return std::get<BumperCapsuloid>(bumper) == std::get<BumperCapsuloid>(bn.bumper);
		case PRYSMOID: return std::get<BumperPrysmoid>(bumper) == std::get<BumperPrysmoid>(bn.bumper);
		case COMPOSITE: return std::get<CompositeBumper>(bumper) == std::get<CompositeBumper>(bn.bumper);
		default: return false;
	}

	return false;
}

void printSphere(const Sphere& s)
{
	std::cout << "Center: " << s.center.x << " " << s.center.y << " " << s.center.z << " " << s.radius << std::endl;
}

BumperQuad::BumperQuad(const FourSpheres &fs, int direction)
{
	upperPlane = fs.upperPlane;
	midPlane = fs.midPlane;

	sphereIndex[0] = fs.indices[0];
	sphereIndex[1] = fs.indices[1];
	sphereIndex[2] = fs.indices[2];
	sphereIndex[3] = fs.indices[3];

	std::cout << "Prysm1 " << std::endl;
	printSphere(fs.spheres[0]);
	printSphere(fs.spheres[1]);
	printSphere(fs.spheres[2]);
	std::cout << "Prysm2 "  << std::endl;
	printSphere(fs.spheres[0]);
	printSphere(fs.spheres[2]);
	printSphere(fs.spheres[3]);

	auto bp = BumperPrysmoid(sphereIndex[0], sphereIndex[1], sphereIndex[2], fs.spheres[0], fs.spheres[1], fs.spheres[2], direction);
	auto bp1 = BumperPrysmoid(sphereIndex[0], sphereIndex[2], sphereIndex[3], fs.spheres[0], fs.spheres[2], fs.spheres[3], direction);

	getSidePlanes(bp, bp1, direction);
}

void BumperQuad::getSidePlanes(const BumperPrysmoid &bp, const BumperPrysmoid &bp1, int direction)
{
	sidePlanes[0] = bp.m2;
	sidePlanes[1] = bp.m0;
	sidePlanes[2] = bp1.m0;
	sidePlanes[3] = bp1.m1;
}

bool CompositeBumper::operator == (const CompositeBumper &bn) const
{
	for (int i = 0; i < MAX_GRID_BUMPERS; i++)
		if (bumperIndices[i] != bn.bumperIndices[i] || flags[i] != bn.flags[i])
			return false;

	return true;
}

std::string CompositeBumper::serialize() const
{
	std::ostringstream oss;

	int size = 0;
	for (int bumperIndice : bumperIndices)
		if (bumperIndice != -1)
			size++;

	oss << size << " ";
	for (int i = 0; i < MAX_GRID_BUMPERS; i ++)
	{
		const int bumperIdx = bumperIndices[i];
		const short flag = flags[i];

		if (bumperIdx == -1)
			break;

		std::ostringstream hexStream;
		hexStream << std::hex << std::uppercase << std::setw(2) << std::setfill('0') << static_cast<int>(flag);
		std::string hexStr = hexStream.str();

		if (flag == 0)
			oss << bumperIdx << " ";
		else
			oss << bumperIdx << "|" << hexStr << " ";
	}

	return oss.str();
}

std::set<int> BumperGrid::descendents (int i) const
{
	std::set<int> res;

	auto bn = bumper[i];

	if (bn.shapeType == Bumper::CAPSULOID)
	{
		auto bc = std::get<BumperCapsuloid>(bn.bumper);
		res.insert(bc.neibExtreme[0]);
		res.insert(bc.neibExtreme[1]);
	}
	else if (bn.shapeType == Bumper::PRYSMOID)
	{
		auto bp = std::get<BumperPrysmoid>(bn.bumper);

		for (int k : bp.neibSide)
		{
			res.insert(k);

			auto bc = std::get<BumperCapsuloid>(bumper[k].bumper);

			res.insert(bc.neibExtreme[0]);
			res.insert(bc.neibExtreme[1]);
		}
	}
	else if (bn.shapeType == Bumper::QUAD)
	{
		auto bq = std::get<BumperQuad>(bn.bumper);

		for (int k : bq.neibSide)
		{
			res.insert(k);

			auto bc = std::get<BumperCapsuloid>(bumper[k].bumper);

			res.insert(bc.neibExtreme[0]);
			res.insert(bc.neibExtreme[1]);
		}
	}

	return res;
}

void BumperGrid::sortByType()
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
		else if (bn.shapeType == Bumper::QUAD)
		{
			auto& bq = std::get<BumperQuad>(bn.bumper);

			bq.neibSide[0] = inversePermutation[bq.neibSide[0]];
			bq.neibSide[1] = inversePermutation[bq.neibSide[1]];
			bq.neibSide[2] = inversePermutation[bq.neibSide[2]];
			bq.neibSide[3] = inversePermutation[bq.neibSide[3]];

			bq.neibOpp = inversePermutation[bq.neibOpp];
		}
	}

	const std::vector<Bumper> oldBumperNode = bumper;
	for (const int i : permutation)
		bumper[i] = oldBumperNode[permutation[i]];
}
