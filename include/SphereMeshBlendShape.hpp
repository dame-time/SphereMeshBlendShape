#pragma once

#include <glm/glm.hpp>

#include <vector>
#include <variant>

class Sphere {
public:
	int index;
	
	glm::vec3 c;
	float r;
	
	Sphere() = default;
	Sphere(const glm::vec3& center, float radius) : c(center), r(radius) {}
};

struct Triangle {
	glm::vec3 n;
	glm::vec3 a, b, c;
};

struct Capsuloid {
	int indices[2];
};

struct Prysmoid {
	int indices[3];
};

struct AABB
{
	glm::vec3 minCorner = glm::vec3(FLT_MAX, FLT_MAX, FLT_MAX);
	glm::vec3 maxCorner = glm::vec3(-FLT_MAX, -FLT_MAX, -FLT_MAX);
	
	void addPoint(const glm::vec3& p)
	{
		for (int i = 0; i < 3; i++)
			if (p[i] > maxCorner[i])
				maxCorner[i] = p[i];
		
		for (int i = 0; i < 3; i++)
			if (p[i] < minCorner[i])
				minCorner[i] = p[i];
	}
	
	[[nodiscard]] glm::vec3 getRandomInternalPos() const
	{
		return glm::vec3(
				minCorner.x + static_cast<float>(rand()) / RAND_MAX * (maxCorner.x - minCorner.x),
				minCorner.y + static_cast<float>(rand()) / RAND_MAX * (maxCorner.y - minCorner.y),
				minCorner.z + static_cast<float>(rand()) / RAND_MAX * (maxCorner.z - minCorner.z)
		);
	}
	
	[[nodiscard]] glm::vec3 BDD() const
	{
		return maxCorner - minCorner;
	}
};

class SphereMesh {
public:
	AABB bbox;
	
	std::vector<Sphere> spheres;
	std::vector<Capsuloid> capsuloids;
	std::vector<Prysmoid> prysmoids;
	
	bool loadFromFile(const char* filename);

private:
	void getBBox();
	
	static Sphere extractSphereFromString(const std::string &s);
	Capsuloid extractCapsuloidFromString(const std::string &s);
	Prysmoid extractPrysmoidFromString(const std::string &s);
};

class SphereMeshBlendShape{
public:
	SphereMesh extract(float t); // t is in 0 to n-1
	bool loadFromFile(const char* filename);
	std::vector<SphereMesh> shapes; // n elements
private:
};

class BumperPrysmoid {
public:
	glm::vec4 basePlane;
	glm::vec3 normal;
	glm::mat3x3 invM;
	
	int sphereIndex[3] {-1, -1, -1};
	
	int neibSide[3] {-1, -1, -1};
	int niebOpp;
	
	BumperPrysmoid(
			int indexA,
			int indexB,
			int indexC,
			const Sphere& a,
			const Sphere& b,
			const Sphere& c,
			int direction);

private:
	[[nodiscard]] glm::vec3 getNormal(const Sphere& a, const Sphere& b, const Sphere& c, int direction) const;
	glm::mat3x3 getM(const Sphere& a, const Sphere& b, const Sphere& c);
	glm::vec4 getPlane(const Sphere& a, const Sphere& b, const Sphere& c, int direction);
};

class BumperCapsuloid {
public:
	float k{};
	
	int sphereIndex[2]{ -1, -1 };
	
	int neibExtreme[2] {-1, -1};
	std::vector<int> neibTriangles;
	
	BumperCapsuloid(int indexA, int indexB, const Sphere& a, const Sphere& b);
};

class BumperSphere {
public:
	int sphereIndex {-1};
	
	std::vector<int> neibCapsules;
};

class BumperNode {
public:
	std::variant<BumperPrysmoid, BumperCapsuloid, BumperSphere> bumper;
	
	enum {
		PRYSMOID,
		CAPSULOID,
		SPHERE
	} shapeType;
	
	BumperNode() : bumper(BumperSphere()) {};
};

class BumperGraph {
public:
	std::vector<BumperNode> bumperNode;
	std::vector<Sphere> sphere;
	
	void constructFrom(const SphereMesh& sm);
	
	bool pushOutside(glm::vec3& p, glm::vec3& n, int& bumperIndex);

private:
	static float PUSH_EPSILON;
	
	bool pushOutsideBruteForce(glm::vec3& p, glm::vec3& n, int& bumperIndex);
	
	bool pushOutsideCapsuloid(glm::vec3& p, glm::vec3& n, int& bumperIndex);
	bool pushOutsidePrysmoid(glm::vec3& p, glm::vec3& n, int& bumperIndex);
	bool pushOutsidePrysmoidPair(glm::vec3& p, glm::vec3& n, int& bumperIndex);
	bool pushOutsideSphere(glm::vec3& p, glm::vec3& n, int& bumperIndex);
};