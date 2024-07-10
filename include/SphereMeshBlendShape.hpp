#pragma once

#include <glm/glm.hpp>

#include <vector>

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
	
	BumperPrysmoid() = default;
	BumperPrysmoid(
			int indexA,
			int indexB,
			int indexC,
			const Sphere& a,
			const Sphere& b,
			const Sphere& c,
			int direction);

private:
	glm::vec3 getNormal(const Sphere& a, const Sphere& b, const Sphere& c, int direction) const;
	glm::mat3x3 getM(const Sphere& a, const Sphere& b, const Sphere& c);
	glm::vec4 getPlane(const Sphere& a, const Sphere& b, const Sphere& c, int direction);
};

class BumperCapsuloid {
public:
	float k{};
	
	int sphereIndex[2]{ -1, -1 };
	
	BumperCapsuloid() = default;
	BumperCapsuloid(int indexA, int indexB, const Sphere& a, const Sphere& b);
};

class BumperSphere {
public:
	int sphereIndex {-1};
	
	BumperSphere() = default;
};

class BumperNode {
public:
	union
	{
		BumperPrysmoid prysmoid;
		BumperCapsuloid capsuloid;
		BumperSphere sphere;
	} bumper;
	
	enum {
		PRYSMOID,
		CAPSULOID,
		SPHERE
	} shapeType;
	
	std::vector<int> neighbours;
	
	BumperNode() = delete;
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
	bool pushOutsideSphere(glm::vec3& p, glm::vec3& n, int& bumperIndex);
	
	bool doesShareEdge(BumperCapsuloid& a, BumperCapsuloid& b);
	bool doesShareEdge(BumperPrysmoid& a, BumperPrysmoid& b);
	bool doesShareEdge(BumperCapsuloid& a, BumperPrysmoid& b);
};