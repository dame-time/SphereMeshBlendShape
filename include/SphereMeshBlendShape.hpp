#pragma once

#include <glm/glm.hpp>

#include <vector>

struct Sphere {
	glm::vec3 c;
	float r;
};

struct Capsuloid {
	int a, b;
};

struct Prysmoid {
	int a, b, c;
};

class SphereMesh {
public:
	bool loadFromFile(const char* filename);
	glm::vec3 pushOutside(glm::vec3 p, glm::vec3 &normal);
	glm::vec3 pushOutside(glm::vec3 p, glm::vec3 &normal, int &lastCollided); // if lastCollided == -1, test all

private:
	Sphere extractSphereFromString(const std::string &s);
	Capsuloid extractCapsuloidFromString(const std::string &s);
	Prysmoid extractPrysmoidFromString(const std::string &s);
	
	std::vector<Sphere> spheres;
	std::vector<Capsuloid> capsuloids;
	std::vector<Prysmoid> prysmoids;
};

class SphereMeshBlendShape{
public:
	SphereMesh extract(float t); // t is in 0 to n-1
	bool loadFromFile(const char* filename);
	std::vector<SphereMesh> shapes; // n elements
private:
};