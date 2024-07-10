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

void removeElement(std::vector<int>& vec, int element) {
	vec.erase(std::remove(vec.begin(), vec.end(), element), vec.end());
}

void BumperGraph::constructFrom (const SphereMesh &sm)
{
	sphere = sm.spheres;
	
	for (int i = 0; i < sphere.size(); i++)
	{
		BumperSphere bs = BumperSphere();
		bs.sphereIndex = i;
		
		BumperNode node {};
		node.shapeType = BumperNode::SPHERE;
		node.bumper.sphere = bs;
		
		bumperNode.push_back(node);
	}
	
	for(const Capsuloid &c : sm.capsuloids)
	{
		BumperCapsuloid bc = BumperCapsuloid(
				c.indices[0],
				c.indices[1],
				sphere[c.indices[0]],
				sphere[c.indices[1]]);
		
		BumperNode node {};
		node.shapeType = BumperNode::CAPSULOID;
		node.bumper.capsuloid = bc;
		
		bumperNode.push_back(node);
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
		nodeUpper.bumper.prysmoid = bpUpper;
		
		BumperNode nodeLower {};
		nodeLower.shapeType = BumperNode::PRYSMOID;
		nodeLower.bumper.prysmoid = bpLower;
		
		BumperNode edge0 {};
		edge0.shapeType = BumperNode::CAPSULOID;
		edge0.bumper.capsuloid = bc0;
		
		BumperNode edge1 {};
		edge1.shapeType = BumperNode::CAPSULOID;
		edge1.bumper.capsuloid = bc1;
		
		BumperNode edge2 {};
		edge2.shapeType = BumperNode::CAPSULOID;
		edge2.bumper.capsuloid = bc2;
		
		bumperNode.push_back(nodeUpper);
		bumperNode.push_back(nodeLower);
		bumperNode.push_back(edge0);
		bumperNode.push_back(edge1);
		bumperNode.push_back(edge2);
	}
	
	for (int i = 0; i < bumperNode.size(); i++)
	{
		BumperNode &node = bumperNode[i];
		
		for (int j = 0; j < bumperNode.size(); j++)
		{
			if (i == j)
				continue;
			
			BumperNode &neighbour = bumperNode[j];
			
			bool check0 = node.shapeType == BumperNode::CAPSULOID && neighbour.shapeType == BumperNode::CAPSULOID
			              && doesShareEdge(node.bumper.capsuloid, neighbour.bumper.capsuloid);
			bool check1 = node.shapeType == BumperNode::CAPSULOID && neighbour.shapeType == BumperNode::PRYSMOID
			              && doesShareEdge(node.bumper.capsuloid, neighbour.bumper.prysmoid);
			bool check2 = node.shapeType == BumperNode::PRYSMOID && neighbour.shapeType == BumperNode::CAPSULOID
			              && doesShareEdge(neighbour.bumper.capsuloid, node.bumper.prysmoid);
			bool check3 = node.shapeType == BumperNode::PRYSMOID && neighbour.shapeType == BumperNode::CAPSULOID
			              && doesShareEdge(neighbour.bumper.prysmoid, node.bumper.prysmoid);
			
			if (check0 || check1 || check2 || check3)
				node.neighbours.push_back(j);
		}
		
		// Removing duplicates and cycles
		std::sort(node.neighbours.begin(), node.neighbours.end());
		auto last = std::unique(node.neighbours.begin(), node.neighbours.end());
		node.neighbours.erase(last, node.neighbours.end());
		
		removeElement(node.neighbours,i);
	}
}

bool BumperGraph::doesShareEdge (BumperCapsuloid &a, BumperCapsuloid &b)
{
	return (a.sphereIndex[0] == b.sphereIndex[0] || a.sphereIndex[0] == b.sphereIndex[1])
	       || (a.sphereIndex[1] == b.sphereIndex[0] || a.sphereIndex[1] == b.sphereIndex[1]);
}

bool BumperGraph::doesShareEdge (BumperPrysmoid &a, BumperPrysmoid &b)
{
	return (a.sphereIndex[0] == b.sphereIndex[0] || a.sphereIndex[0] == b.sphereIndex[1] || a.sphereIndex[0] == b.sphereIndex[2])
	       || (a.sphereIndex[1] == b.sphereIndex[0] || a.sphereIndex[1] == b.sphereIndex[1] || a.sphereIndex[1] == b.sphereIndex[2])
	       || (a.sphereIndex[2] == b.sphereIndex[0] || a.sphereIndex[2] == b.sphereIndex[1] || a.sphereIndex[2] == b.sphereIndex[2]);
}

bool BumperGraph::doesShareEdge (BumperCapsuloid &a, BumperPrysmoid &b)
{
	return (a.sphereIndex[0] == b.sphereIndex[0] || a.sphereIndex[0] == b.sphereIndex[1] || a.sphereIndex[0] == b.sphereIndex[2])
	       || (a.sphereIndex[1] == b.sphereIndex[0] || a.sphereIndex[1] == b.sphereIndex[1] || a.sphereIndex[1] == b.sphereIndex[2]);
}

bool BumperGraph::pushOutsideSphere(glm::vec3& p, glm::vec3& n, int& bumperIndex)
{
	BumperNode& node = bumperNode[bumperIndex];
	Sphere& s = sphere[node.bumper.sphere.sphereIndex];
	
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
	BumperCapsuloid& bc = node.bumper.capsuloid;
	
	Sphere& a = sphere[bc.sphereIndex[0]];
	Sphere& b = sphere[bc.sphereIndex[1]];
	
	glm::vec3 d = b.c - a.c;
	float t = glm::dot(p - a.c, d) / glm::dot(d, d);
	glm::vec3 q = (1 - t) * a.c + t * b.c;
	
	float pq = glm::length(p - q);
	
	float kd = node.bumper.capsuloid.k * pq;
	float t_prime = t + kd;
	
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
	BumperPrysmoid& bp = node.bumper.prysmoid;
	
	Sphere& a = sphere[bp.sphereIndex[0]];
	Sphere& b = sphere[bp.sphereIndex[1]];
	Sphere& c = sphere[bp.sphereIndex[2]];
	
	glm::vec3 basePlaneNormal = glm::vec3(bp.basePlane.x, bp.basePlane.y, bp.basePlane.z);
	bool isOver = glm::dot(basePlaneNormal, p - a.c) > 0;
	
	if (!isOver)
		return false;
	
	glm::vec3 barycentricCoordWithDist = bp.invM * (p - a.c);
	float thirdCoord = 1.0f - barycentricCoordWithDist.x - barycentricCoordWithDist.y;
	
	glm::vec4 coord = glm::vec4(barycentricCoordWithDist.x, barycentricCoordWithDist.y, thirdCoord,
	                            barycentricCoordWithDist.z);
	
	bool isInside = coord.x >= 0
	                && coord.y >= 0
	                && coord.z >= 0;
	
	if (!isInside)
		return false;
	
	float sphereRadius = coord.x * b.r + coord.y * c.r + coord.z * a.r;
	
	if (coord.w > sphereRadius)
		return false;
	
	glm::vec3 sphereCenter = coord.x * b.c + coord.y * c.c + coord.z * a.c;
	
	n = bp.normal;
	p = sphereCenter + n * (sphereRadius + PUSH_EPSILON);
	
	return true;
}

bool BumperGraph::pushOutsideBruteForce(glm::vec3& p, glm::vec3& n, int& bumperIndex)
{
	bool hasCollided = false;
	
	for (int i = 0; i < bumperNode.size(); i++)
	{
		BumperNode& node = bumperNode[i];
		
		if (node.shapeType == BumperNode::SPHERE)
		{
			if (pushOutsideSphere(p, n, i))
			{
				bumperIndex = i;
				hasCollided = true;
			}
		}
		else if (node.shapeType == BumperNode::CAPSULOID)
		{
			if (pushOutsideCapsuloid(p, n, i))
			{
				bumperIndex = i;
				hasCollided = true;
			}
		}
		else if (node.shapeType == BumperNode::PRYSMOID)
		{
			if (pushOutsidePrysmoid(p, n, i))
			{
				bumperIndex = i;
				hasCollided = true;
			}
		}
	}
	
	return hasCollided;
}

bool BumperGraph::pushOutside (glm::vec3 &p, glm::vec3 &n, int &bumperIndex)
{
	if (bumperIndex == -1)
		return pushOutsideBruteForce(p, n, bumperIndex);
	
	BumperNode& node = bumperNode[bumperIndex];
	std::vector<int> neighbours = node.neighbours;
	neighbours.insert(neighbours.begin(), bumperIndex); // testing also against the old
	
	bool hasCollided = false;
	
	for (int& i : neighbours)
	{
		BumperNode& neighbour = bumperNode[i];
		
		if (neighbour.shapeType == BumperNode::SPHERE)
		{
			if (pushOutsideSphere(p, n, i))
			{
				bumperIndex = i;
				hasCollided = true;
			}
		}
		else if (neighbour.shapeType == BumperNode::CAPSULOID)
		{
			if (pushOutsideCapsuloid(p, n, i))
			{
				bumperIndex = i;
				hasCollided = true;
			}
		}
		else if (neighbour.shapeType == BumperNode::PRYSMOID)
		{
			if (pushOutsidePrysmoid(p, n, i))
			{
				bumperIndex = i;
				hasCollided = true;
			}
		}
	}
	
	return hasCollided;
}