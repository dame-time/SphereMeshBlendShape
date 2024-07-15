#include "SphereMeshBlendShape.hpp"

#include <fstream>
#include <iostream>
#include <string>

float BumperGraph::PUSH_EPSILON = 0.000001f;

void SphereMesh::getBBox ()
{
	bbox.minCorner = glm::vec3(FLT_MAX, FLT_MAX, FLT_MAX);
	bbox.maxCorner = glm::vec3(-FLT_MAX, -FLT_MAX, -FLT_MAX);
	
	for (const auto& s : spheres)
	{
		bbox.addPoint(s.c + glm::vec3(s.r, s.r, s.r));
		bbox.addPoint(s.c - glm::vec3(s.r, s.r, s.r));
	}
}

bool SphereMesh::loadFromFile (const char *filename)
{
	std::ifstream inputFile(filename);
	
	if (!inputFile.is_open()) {
		std::cout << "Error opening the file!" << std::endl;
		return false;
	}
	
	std::string line;
	std::vector<std::string> content;
	
	while (std::getline(inputFile, line))
		content.push_back(line);
	
	inputFile.close();
	
	if (content.size() < 5)
	{
		std::cout << "Invalid file format!" << std::endl;
		return false;
	}
	
	std::string structCounter = content[2];
	std::string delimiter = " ";
	int nSpheres = 0, nCapsuloids = 0, nPrysmoids = 0;
	for (int i = 0; i < 3; i++) {
		std::string token = structCounter.substr(0, structCounter.find(delimiter));
		structCounter.erase(0, structCounter.find(delimiter) + delimiter.length());
		
		if (i == 0)
			nSpheres = std::stoi(token);
		else if (i == 1)
			nPrysmoids = std::stoi(token);
		else
			nCapsuloids = std::stoi(token);
	}
	
	for (int i = 0; i < nSpheres; i++)
	{
		Sphere s = extractSphereFromString(content[4 + i]);
		spheres.push_back(s);
	}
	
	for (int i = 0; i < nPrysmoids; i++)
	{
		Prysmoid p = extractPrysmoidFromString(content[4 + nSpheres + i]);
		prysmoids.push_back(p);
	}
	
	for (int i = 0; i < nCapsuloids; i++)
	{
		Capsuloid c = extractCapsuloidFromString(content[4 + nSpheres + nPrysmoids + i]);
		capsuloids.push_back(c);
	}
	
	getBBox();
	
	return true;
}

Sphere SphereMesh::extractSphereFromString(const std::string &sphereString)
{
	std::string s = sphereString;
	
	Sphere sphere{};
	
	std::string delimiter = " ";
	std::string token = s.substr(0, s.find(delimiter));
	sphere.c.x = std::stof(token);
	s.erase(0, s.find(delimiter) + delimiter.length());
	
	token = s.substr(0, s.find(delimiter));
	sphere.c.y = std::stof(token);
	s.erase(0, s.find(delimiter) + delimiter.length());
	
	token = s.substr(0, s.find(delimiter));
	sphere.c.z = std::stof(token);
	s.erase(0, s.find(delimiter) + delimiter.length());
	
	token = s.substr(0, s.find(delimiter));
	sphere.r = std::stof(token);
	
	return sphere;
}

Capsuloid SphereMesh::extractCapsuloidFromString(const std::string &capsuloidString)
{
	std::string s = capsuloidString;
	
	std::string delimiter = " ";
	std::string token = s.substr(0, s.find(delimiter));
	int sphereA = std::stoi(token);
	s.erase(0, s.find(delimiter) + delimiter.length());
	
	token = s.substr(0, s.find(delimiter));
	int sphereB = std::stoi(token);
	
	Capsuloid capsuloid;
	capsuloid.indices[0] = sphereA;
	capsuloid.indices[1] = sphereB;
	
	return capsuloid;
}

Prysmoid SphereMesh::extractPrysmoidFromString(const std::string &prysmoidString)
{
	std::string s = prysmoidString;
	
	
	std::string delimiter = " ";
	std::string token = s.substr(0, s.find(delimiter));
	int prysmoidA = std::stoi(token);
	s.erase(0, s.find(delimiter) + delimiter.length());
	
	token = s.substr(0, s.find(delimiter));
	int prysmoidB = std::stoi(token);
	s.erase(0, s.find(delimiter) + delimiter.length());
	
	token = s.substr(0, s.find(delimiter));
	int prysmoidC = std::stoi(token);
	
	Prysmoid prysmoid;
	
	prysmoid.indices[0] = prysmoidA;
	prysmoid.indices[1] = prysmoidB;
	prysmoid.indices[2] = prysmoidC;
	
	return prysmoid;
}

bool SphereMeshBlendShape::loadFromFile (const char *filename)
{
	SphereMesh sm;
	
	if (!sm.loadFromFile(filename))
		return false;
	
	shapes.push_back(sm);
	
	return true;
}

BumperCapsuloid::BumperCapsuloid (int indexA, int indexB, const Sphere &a, const Sphere &b)
{
	sphereIndex[0] = indexA;
	sphereIndex[1] = indexB;
	
	float BA = glm::length(b.c - a.c);
	float r0r1 = b.r - a.r;
	
	float sq_BA = BA * BA;
	float sq_r0r1 = r0r1 * r0r1;
	
	float l = std::sqrt(sq_BA - sq_r0r1);
	
	k = r0r1 / (l * BA);
}

BumperPrysmoid::BumperPrysmoid (
		int indexA,
		int indexB,
		int indexC,
		const Sphere &a,
		const Sphere &b,
		const Sphere&c,
		int direction)
{
	sphereIndex[0] = indexA;
	sphereIndex[1] = indexB;
	sphereIndex[2] = indexC;
	
	normal = getNormal(a, b, c, direction);
	invM = glm::inverse(getM(a, b, c));
	basePlane = getPlane(a, b, c, direction);
}

glm::vec3 BumperPrysmoid::getNormal (const Sphere &a, const Sphere &b, const Sphere& c, int direction) const
{
	Triangle tri {};
	
	tri.a = a.c;
	tri.b = b.c;
	tri.c = c.c;
	
	float sign = (direction > 0) ? 1 : -1;
	
	tri.n = sign * glm::cross(b.c - a.c, c.c - a.c);
	tri.n = glm::normalize(tri.n);
	
	for (int i = 0; i < 10; i++){
		tri.a = sign * (tri.a + tri.n * a.r);
		tri.b = sign * (tri.b + tri.n * b.r);
		tri.c = sign * (tri.c + tri.n * c.r);
		
		tri.n = sign * glm::cross(b.c - a.c, c.c - a.c);
		tri.n = glm::normalize(tri.n);
	}
	
	return tri.n;
}

glm::mat3x3 BumperPrysmoid::getM (const Sphere& a, const Sphere& b, const Sphere& c)
{
	glm::mat3x3 M;
	
	glm::vec3 col0, col1, col2;
	
	col0 = b.c - a.c;
	col1 = c.c - a.c;
	col2 = normal;
	
	M[0] = col0;
	M[1] = col1;
	M[2] = col2;
	
	return M;
}

glm::vec4 BumperPrysmoid::getPlane (const Sphere &a, const Sphere &b, const Sphere &c, int direction)
{
	glm::vec3 u = b.c - a.c;
	glm::vec3 v = c.c - a.c;
	
	glm::vec3 n = (direction > 0) ? glm::normalize(glm::cross(u, v)) : -glm::normalize(glm::cross(u, v));
	
	float d = -glm::dot(n, a.c);
	
	return {n, d};
}

void BumperGraph::constructFrom (const SphereMesh &sm)
{
	sphere = sm.spheres;
	
	std::vector<std::vector<int>> capsuloidAdj (sphere.size(), std::vector<int>(sphere.size(), -1));
	
	for (int i = 0; i < sphere.size(); i++)
	{
		BumperSphere bs = BumperSphere();
		bs.sphereIndex = i;
		
		BumperNode node {};
		node.shapeType = BumperNode::SPHERE;
		node.bumper = bs;
		
		bumperNode.push_back(node);
	}
	
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
		
		std::get<BumperSphere>(bumperNode[c.indices[0]].bumper).neibCapsules.push_back(bumperNode.size());
		std::get<BumperSphere>(bumperNode[c.indices[1]].bumper).neibCapsules.push_back(bumperNode.size());
		
		BumperNode node {};
		node.shapeType = BumperNode::CAPSULOID;
		node.bumper = bc;
		
		bumperNode.push_back(node);
		
		capsuloidAdj[c.indices[0]][c.indices[1]] = bumperNode.size() - 1;
		capsuloidAdj[c.indices[1]][c.indices[0]] = bumperNode.size() - 1;
	}
	
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
		
		BumperNode nodeUpper {};
		nodeUpper.shapeType = BumperNode::PRYSMOID;
		bpUpper.niebOpp = bumperNode.size() + 1;
		nodeUpper.bumper = bpUpper;
		
		BumperNode nodeLower {};
		nodeLower.shapeType = BumperNode::PRYSMOID;
		bpLower.niebOpp = bumperNode.size();
		nodeLower.bumper = bpLower;
		
		BumperNode edge0 {};
		edge0.shapeType = BumperNode::CAPSULOID;
		bc0.neibExtreme[0] = p.indices[0];
		bc0.neibExtreme[1] = p.indices[1];
		edge0.bumper = bc0;
		
		BumperNode edge1 {};
		edge1.shapeType = BumperNode::CAPSULOID;
		bc1.neibExtreme[0] = p.indices[0];
		bc1.neibExtreme[1] = p.indices[2];
		edge1.bumper = bc1;
		
		BumperNode edge2 {};
		edge2.shapeType = BumperNode::CAPSULOID;
		bc2.neibExtreme[0] = p.indices[1];
		bc2.neibExtreme[1] = p.indices[2];
		edge2.bumper = bc2;
		
		bumperNode.push_back(nodeUpper);
		bumperNode.push_back(nodeLower);
		
		int upperIndex = bumperNode.size() - 2;
		int lowerIndex = bumperNode.size() - 1;
		
		if (capsuloidAdj[p.indices[0]][p.indices[1]] == -1 || capsuloidAdj[p.indices[1]][p.indices[0]] == -1)
		{
			capsuloidAdj[p.indices[0]][p.indices[1]] = bumperNode.size();
			capsuloidAdj[p.indices[1]][p.indices[0]] = bumperNode.size();
			
			std::get<BumperSphere>(bumperNode[p.indices[0]].bumper).neibCapsules.push_back(bumperNode.size());
			std::get<BumperSphere>(bumperNode[p.indices[1]].bumper).neibCapsules.push_back(bumperNode.size());
			
			bumperNode.push_back(edge0);
		}
		
		int idx = capsuloidAdj[p.indices[0]][p.indices[1]];
		
		std::get<BumperPrysmoid>(bumperNode[upperIndex].bumper).neibSide[0] = idx;
		std::get<BumperPrysmoid>(bumperNode[lowerIndex].bumper).neibSide[0] = idx;
		
		std::get<BumperCapsuloid>(bumperNode[idx].bumper).neibTriangles.push_back(upperIndex);
		std::get<BumperCapsuloid>(bumperNode[idx].bumper).neibTriangles.push_back(lowerIndex);
		
		if (capsuloidAdj[p.indices[0]][p.indices[2]] == -1 || capsuloidAdj[p.indices[2]][p.indices[0]] == -1)
		{
			capsuloidAdj[p.indices[0]][p.indices[2]] = bumperNode.size();
			capsuloidAdj[p.indices[2]][p.indices[0]] = bumperNode.size();
			
			std::get<BumperSphere>(bumperNode[p.indices[0]].bumper).neibCapsules.push_back(bumperNode.size());
			std::get<BumperSphere>(bumperNode[p.indices[2]].bumper).neibCapsules.push_back(bumperNode.size());
			
			bumperNode.push_back(edge1);
		}
		
		idx = capsuloidAdj[p.indices[0]][p.indices[2]];
		
		std::get<BumperPrysmoid>(bumperNode[upperIndex].bumper).neibSide[1] = idx;
		std::get<BumperPrysmoid>(bumperNode[lowerIndex].bumper).neibSide[1] = idx;
		
		std::get<BumperCapsuloid>(bumperNode[idx].bumper).neibTriangles.push_back(upperIndex);
		std::get<BumperCapsuloid>(bumperNode[idx].bumper).neibTriangles.push_back(lowerIndex);
		
		if (capsuloidAdj[p.indices[1]][p.indices[2]] == -1 || capsuloidAdj[p.indices[2]][p.indices[1]] == -1)
		{
			capsuloidAdj[p.indices[1]][p.indices[2]] = bumperNode.size();
			capsuloidAdj[p.indices[2]][p.indices[1]] = bumperNode.size();
			
			std::get<BumperSphere>(bumperNode[p.indices[1]].bumper).neibCapsules.push_back(bumperNode.size());
			std::get<BumperSphere>(bumperNode[p.indices[2]].bumper).neibCapsules.push_back(bumperNode.size());
			
			bumperNode.push_back(edge2);
		}
		
		idx = capsuloidAdj[p.indices[1]][p.indices[2]];
		
		std::get<BumperPrysmoid>(bumperNode[upperIndex].bumper).neibSide[2] = idx;
		std::get<BumperPrysmoid>(bumperNode[lowerIndex].bumper).neibSide[2] = idx;
		
		std::get<BumperCapsuloid>(bumperNode[idx].bumper).neibTriangles.push_back(upperIndex);
		std::get<BumperCapsuloid>(bumperNode[idx].bumper).neibTriangles.push_back(lowerIndex);
	}
}

bool BumperGraph::pushOutsideSphere(glm::vec3& p, glm::vec3& n, int& bumperIndex)
{
	BumperNode& node = bumperNode[bumperIndex];
	BumperSphere bs = std::get<BumperSphere>(node.bumper);
	Sphere& s = sphere[bs.sphereIndex];
	
	bool isInside = glm::length(p - s.c) < s.r;
	
	if (!isInside)
		return false;
	
	n = glm::normalize(p - s.c);
	p = s.c + (p - s.c) * (s.r + PUSH_EPSILON) / glm::length(p - s.c);
	
	return true;
}

bool BumperGraph::pushOutsideCapsuloid (glm::vec3 &p, glm::vec3 &n, int &bumperIndex)
{
	BumperNode& node = bumperNode[bumperIndex];
	BumperCapsuloid bc = std::get<BumperCapsuloid>(node.bumper);
	
	Sphere& a = sphere[bc.sphereIndex[0]];
	Sphere& b = sphere[bc.sphereIndex[1]];
	
	glm::vec3 d = b.c - a.c;
	float t = glm::dot(p - a.c, d) / glm::dot(d, d);
	glm::vec3 q = (1 - t) * a.c + t * b.c;
	
	float pq = glm::length(p - q);
	
	float kd = bc.k * pq;
	float t_prime = t + kd;
	
	if (t_prime < 0)
	{
		int index = bc.sphereIndex[0];
		if (!pushOutsideSphere(p, n, index))
			return false;
		
		bumperIndex = index;
		return true;
	}
	else if (t_prime > 1)
	{
		int index = bc.sphereIndex[1];
		if(!pushOutsideSphere(p, n, index))
			return false;
		
		bumperIndex = index;
		return true;
	}
	
	t_prime = glm::clamp(t_prime, 0.0f, 1.0f);
	
	glm::vec3 q_prime = (1 - t_prime) * a.c + t_prime * b.c;
	float r_prime = (1 - t_prime) * a.r + t_prime * b.r;
	
	float dist = glm::length(p - q_prime);
	
	if (dist >= r_prime)
		return false;
	
	n = glm::normalize(p - q_prime);
	p = q_prime + n * (r_prime + PUSH_EPSILON);
	
	return true;
}

bool BumperGraph::pushOutsidePrysmoid (glm::vec3 &p, glm::vec3 &n, int &bumperIndex)
{
	BumperNode& node = bumperNode[bumperIndex];
	BumperPrysmoid bp = std::get<BumperPrysmoid>(node.bumper);
	
	Sphere& a = sphere[bp.sphereIndex[0]];
	Sphere& b = sphere[bp.sphereIndex[1]];
	Sphere& c = sphere[bp.sphereIndex[2]];
	
	glm::vec3 basePlaneNormal = glm::vec3(bp.basePlane.x, bp.basePlane.y, bp.basePlane.z);
	
	glm::vec3 barycentricCoordWithDist = bp.invM * (p - a.c);
	float thirdCoord = 1.0f - barycentricCoordWithDist.x - barycentricCoordWithDist.y;
	
	glm::vec4 coord = glm::vec4(barycentricCoordWithDist.x, barycentricCoordWithDist.y, thirdCoord,
	                            barycentricCoordWithDist.z);
	
	if (coord.x < 0)
	{
		bumperIndex = bp.neibSide[1];
		return pushOutsideCapsuloid(p, n, bumperIndex);
	}
	else if (coord.y < 0)
	{
		bumperIndex = bp.neibSide[0];
		return pushOutsideCapsuloid(p, n, bumperIndex);
	}
	else if (coord.z < 0)
	{
		bumperIndex = bp.neibSide[2];
		return pushOutsideCapsuloid(p, n, bumperIndex);
	}
	
	float sphereRadius = coord.x * b.r + coord.y * c.r + coord.z * a.r;
	
	if (coord.w > sphereRadius)
		return false;
	
	glm::vec3 sphereCenter = coord.x * b.c + coord.y * c.c + coord.z * a.c;
	
	n = bp.normal;
	p = sphereCenter + n * (sphereRadius + PUSH_EPSILON);
	
	return true;
}

bool BumperGraph::pushOutsidePrysmoidPair (glm::vec3 &p, glm::vec3 &n, int &bumperIndex)
{
	BumperNode& node = bumperNode[bumperIndex];
	BumperPrysmoid bp = std::get<BumperPrysmoid>(node.bumper);
	
	Sphere& a = sphere[bp.sphereIndex[0]];
	Sphere& b = sphere[bp.sphereIndex[1]];
	Sphere& c = sphere[bp.sphereIndex[2]];
	
	glm::vec3 basePlaneNormal = glm::vec3(bp.basePlane.x, bp.basePlane.y, bp.basePlane.z);
	bool sideUp = glm::dot(basePlaneNormal, p - a.c) > 0;
	
	if (!sideUp)
	{
		bumperIndex = bp.niebOpp;
		return pushOutsidePrysmoid(p, n, bumperIndex);
	}
	
	return pushOutsidePrysmoid(p, n, bumperIndex);
}

bool BumperGraph::pushOutsideBruteForce(glm::vec3& p, glm::vec3& n, int& bumperIndex)
{
	bool hasCollided = false;
	
	for (int i = 0; i < bumperNode.size(); i++)
	{
		int index = i;
		
		if (bumperNode[i].shapeType == BumperNode::SPHERE)
		{
			if (std::get<BumperSphere>(bumperNode[i].bumper).neibCapsules.empty())
				continue;
			
			if (pushOutsideSphere(p, n, index))
			{
				bumperIndex = index;
				hasCollided = true;
			}
		}
		else if (bumperNode[i].shapeType == BumperNode::CAPSULOID)
		{
//			if (std::get<BumperCapsuloid>(bumperNode[i].bumper).neibTriangles.empty())
//				continue;
			
			if (pushOutsideCapsuloid(p, n, index))
			{
				bumperIndex = index;
				hasCollided = true;
			}
		}
		else if (bumperNode[i].shapeType == BumperNode::PRYSMOID)
		{
			if (pushOutsidePrysmoidPair(p, n, index))
			{
				bumperIndex = index;
				hasCollided = true;
			}
			
			i++; // Skipping the other plane since I already check that
		}
	}
	
	return hasCollided;
}

bool BumperGraph::pushOutside (glm::vec3 &p, glm::vec3 &n, int &bumperIndex)
{
	if (bumperIndex == -1)
		return pushOutsideBruteForce(p, n, bumperIndex);
	
	if (bumperNode[bumperIndex].shapeType == BumperNode::SPHERE)
	{
		int idx = bumperIndex;
		
		if (std::get<BumperSphere>(bumperNode[idx].bumper).neibCapsules.empty())
			return false;
		
		if (!pushOutsideSphere(p, n, bumperIndex))
			return false;
		
		for (int i : std::get<BumperSphere>(bumperNode[idx].bumper).neibCapsules)
		{
			int index = i;
			if (pushOutsideCapsuloid(p, n, index))
			{
				bumperIndex = index;
				return true;
			}
		}
	}
	else if (bumperNode[bumperIndex].shapeType == BumperNode::CAPSULOID)
	{
		int idx = bumperIndex;
		if (!pushOutsideCapsuloid(p, n, bumperIndex))
			return false;
		
		for (int i = 0; i < std::get<BumperCapsuloid>(bumperNode[idx].bumper).neibTriangles.size(); i += 2)
		{
			int index = std::get<BumperCapsuloid>(bumperNode[idx].bumper).neibTriangles[i];
			if(pushOutsidePrysmoidPair(p, n, index))
			{
				bumperIndex = index;
				return true;
			}
		}
	}
	else if (bumperNode[bumperIndex].shapeType == BumperNode::PRYSMOID)
	{
		if (!pushOutsidePrysmoidPair(p, n, bumperIndex))
			return false;
	}
	
	return true;
}