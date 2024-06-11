#pragma once

#include <glm/glm.hpp>

class SphereMesh {
public:
	bool loadFromFile(const char* filename);
	glm::vec3 pushOutside( glm::vec3 p, glm::vec3 &normal);
	glm::vec3 pushOutside( glm::vec3 p, glm::vec3 &normal, int &lastCollided); // if lastCollided == -1, test all

private:
	/// data...
	
};

class SphereMeshBlendShape{
public:
	SphereMesh extract(float t); // t is in 0 to n-1
	bool loadFromFile(const char* filename);
	std::vector<SphereMesh> shapes; // n elements
private:
};