/*
	Licensed under the Apache License, Version 2.0 (the "License");
	you may not use this file except in compliance with the License.
	You may obtain a copy of the License at

		http://www.apache.org/licenses/LICENSE-2.0

	Unless required by applicable law or agreed to in writing, software
	distributed under the License is distributed on an "AS IS" BASIS,
	WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
	See the License for the specific language governing permissions and
	limitations under the License.

	Authors: Davide Paolillo, Marco Tarini
*/

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <algorithm>
#include <unordered_set>
#include <set>

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/norm.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include "sphere_mesh.h"
#include "bumper_grid.h"

using namespace SM;

void SphereMesh::updateBBox ()
{
	bbox.minCorner = glm::vec3(FLT_MAX, FLT_MAX, FLT_MAX);
	bbox.maxCorner = glm::vec3(-FLT_MAX, -FLT_MAX, -FLT_MAX);

	for (const auto& s : spheres)
	{
		bbox.addPoint(s.center + glm::vec3(s.radius, s.radius, s.radius));
		bbox.addPoint(s.center - glm::vec3(s.radius, s.radius, s.radius));
	}
}

Sphere deeperSphere(const Sphere& sphereA, const Sphere& sphereB, const glm::vec3& p)
{
	Sphere s {};

    const Grid::BumperCapsuloid bc (-1, -1, sphereA, sphereB);

	const glm::vec3 pa = p - sphereA.center;
	float t = glm::dot(pa, bc.dNorm);
	const float pq = std::sqrt(glm::dot(pa, pa) - (t * t));

	t +=  bc.k * pq;

	t = glm::clamp(t, 0.0f, bc.dLength);

	if (t == 0.0f) return sphereA;
	if (t == bc.dLength) return sphereB;

	s.center = sphereA.center + t * bc.dNorm;
	s.radius = sphereA.radius + t * bc.dR;

	return s;
}

bool isDegenerateCapsuloid(const Sphere& sphereA, const Sphere& sphereB)
{
	return glm::distance(sphereA.center, sphereB.center) <= std::abs(sphereA.radius - sphereB.radius);
}

void SphereMesh::removeDegeneratePrysmoids()
{
	std::vector<Prysmoid> newPrysmoids;
	for (Prysmoid& prysmoid : prysmoids)
	{
		bool canAdd = true;
		for (int j = 0; j < 3; j++)
		{
			const int indexA = prysmoid.indices[j];
			const int indexB = prysmoid.indices[(j + 1) % 3];
			const int indexC = prysmoid.indices[(j + 2) % 3];

			Sphere s = deeperSphere(spheres[indexA], spheres[indexB], spheres[indexC].center);
			if (isDegenerateCapsuloid(spheres[indexC], s))
			{
				canAdd = false;

				if (spheres[indexC].radius < s.radius)
					capsuloids.push_back({indexA, indexB});
				else
				{
					capsuloids.push_back({indexA, indexC});
					capsuloids.push_back({indexB, indexC});
				}
			}
		}
		if (canAdd)
			newPrysmoids.emplace_back(prysmoid);
	}
	prysmoids = newPrysmoids;
}

void SphereMesh::removeDegenerateCapsuloids()
{
	std::vector<Capsuloid> newCapsuloids;
	for (Capsuloid& capsuloid : capsuloids)
	{
		const int indexA = capsuloid.indices[0];
		const int indexB = capsuloid.indices[1];
		if (isDegenerateCapsuloid(spheres[indexA], spheres[indexB]))
		{
			if (spheres[indexA].radius < spheres[indexB].radius)
				singletons.push_back({indexB});
			else
				singletons.push_back({indexA});
			continue;
		}
		newCapsuloids.push_back(capsuloid);
	}
	capsuloids = newCapsuloids;
}

void SphereMesh::removeDegenerateElements()
{
	removeDegeneratePrysmoids();
	removeDegenerateCapsuloids();
}

void SphereMesh::removeRedundantElements()
{
	std::set<std::string> uniquePrysmoids;
	std::set<std::string> uniqueCapsuloids;
	std::set<std::string> uniqueSingletons;
	std::vector<Prysmoid> newPrysmoids;
	std::vector<Capsuloid> newCapsuloids;
	std::vector<Singleton> newSingletons;

	for (auto& prysmoid : prysmoids)
	{
		std::vector<int> indices = {prysmoid.indices[0], prysmoid.indices[1], prysmoid.indices[2]};
		std::ranges::sort(indices);
		std::string key = std::to_string(indices[0]) + std::to_string(indices[1]) + std::to_string(indices[2]);
		if (uniquePrysmoids.contains(key))
			continue;

		newPrysmoids.push_back(prysmoid);
		uniquePrysmoids.insert(key);
		std::vector<std::string> capsulesOverPrysmoid = {std::to_string(indices[0]) + std::to_string(indices[1]),
														   std::to_string(indices[1]) + std::to_string(indices[2]),
														   std::to_string(indices[0]) + std::to_string(indices[2])};
		uniqueCapsuloids.insert(capsulesOverPrysmoid[0]);
		uniqueCapsuloids.insert(capsulesOverPrysmoid[1]);
		uniqueCapsuloids.insert(capsulesOverPrysmoid[2]);
		uniqueSingletons.insert(std::to_string(indices[0]));
		uniqueSingletons.insert(std::to_string(indices[1]));
		uniqueSingletons.insert(std::to_string(indices[2]));
	}

	for (auto& capsuloid : capsuloids)
	{
		std::vector<int> indices = {capsuloid.indices[0], capsuloid.indices[1]};
		std::ranges::sort(indices);
		std::string key = std::to_string(indices[0]) + std::to_string(indices[1]);
		if (uniqueCapsuloids.contains(key))
			continue;

		newCapsuloids.push_back(capsuloid);
		uniqueCapsuloids.insert(key);
		uniqueSingletons.insert(std::to_string(indices[0]));
		uniqueSingletons.insert(std::to_string(indices[1]));
	}

	for (auto &[index] : singletons)
	{
		std::string key = std::to_string(index);
		if (uniqueSingletons.contains(key))
			continue;

		newSingletons.push_back({index});
		uniqueSingletons.insert(key);
	}

	prysmoids = newPrysmoids;
	capsuloids = newCapsuloids;
	singletons = newSingletons;
}

std::vector<int> setIntersection(const std::unordered_set<int>& set1, const std::unordered_set<int>& set2) {
	std::vector<int> intersection;
	for (const int& elem : set1)
		if (set2.contains(elem))
			intersection.push_back(elem);

	return intersection;
}

std::vector<int> setDifference(const std::unordered_set<int>& set1, const std::unordered_set<int>& set2) {
	std::vector<int> difference;
	std::unordered_set<int> intersectionSet;

	for (const int& elem : set1)
		if (set2.contains(elem))
			intersectionSet.insert(elem);

	for (const int& elem : set1)
		if (!intersectionSet.contains(elem))
			difference.push_back(elem);

	for (const int& elem : set2)
		if (!intersectionSet.contains(elem))
			difference.push_back(elem);

	return difference;
}

JoinedQuadrilateral SphereMesh::joinPrysmoids(const Prysmoid &p1, const Prysmoid &p2) const
{
	std::unordered_set<int> uniqueIndices;
	std::unordered_set<int> p1Indices;
	std::unordered_set<int> p2Indices;
	for (int i = 0; i < 3; i++)
	{
		uniqueIndices.insert(p1.indices[i]);
		uniqueIndices.insert(p2.indices[i]);
		p1Indices.insert(p1.indices[i]);
		p2Indices.insert(p2.indices[i]);
	}

	const std::vector<int> ab = setIntersection(p1Indices, p2Indices);
	const std::vector<int> cd = setDifference(p1Indices, p2Indices);

	if (ab.size() != 2 || cd.size() != 2)
	{
		std::cerr << "Error joining prysmoids: invalid indices" << std::endl;
		return {};
	}

	glm::vec3 diagonalAB = spheres[ab[0]].center - spheres[ab[1]].center;
	glm::vec3 diagonalCD = spheres[cd[0]].center - spheres[cd[1]].center;
	glm::vec3 newPlaneNormal = glm::normalize(glm::cross(diagonalAB, diagonalCD));

	float k0 = glm::dot(-newPlaneNormal, spheres[ab[0]].center);
	float k1 = glm::dot(-newPlaneNormal, spheres[ab[1]].center);
	float k2 = glm::dot(-newPlaneNormal, spheres[cd[0]].center);
	float k3 = glm::dot(-newPlaneNormal, spheres[cd[1]].center);
	float k = (k0 + k1 + k2 + k3) / 4.0f;

	Plane p {};
	p.n = newPlaneNormal;
	p.k = k;

	Sphere s0 = {p.project(spheres[ab[0]].center), 0};
	Sphere s1 = {p.project(spheres[ab[1]].center), 0};
	Sphere s2 = {p.project(spheres[cd[0]].center), 0};
	Sphere s3 = {p.project(spheres[cd[1]].center), 0};

	return {p, s0, s1, s2, s3, ab[0], ab[1], cd[0], cd[1], 1};
}

void SphereMesh::translateScale(glm::vec3 t , float s){
	for (Sphere& sp : spheres) {
		sp.center = (sp.center + t) * s;
		sp.radius *= s;
	}
}

void SphereMesh::scaleTranslate(glm::vec3 t, float s)
{
	for (Sphere& sp : spheres) {
		sp.center = (sp.center * s) + t;
		sp.radius *= s;
	}

	updateBBox();
}

void SphereMesh::rotoTranslate (const glm::mat3 &rot, const glm::vec3& trasl)
{
	assert(std::abs(glm::determinant(rot) - 1.0f) < 0.0001f);
	for (Sphere& sp:spheres)
		sp.center = rot * sp.center + trasl;

	// Recompute AABB
	updateBBox();
}

void SphereMesh::translate(const glm::vec3 &t)
{
	for (Sphere &s : spheres)
		s.center += t;

	updateBBox();
}

void SphereMesh::scale(const float s)
{
	for (Sphere &sph : spheres)
	{
		sph.center *= s;
		sph.radius *= s;
	}

	updateBBox();
}

void SphereMesh::rotateY(const int angle)
{
	const float radians = glm::radians(static_cast<float>(angle));
	const auto rotationMatrix = glm::mat3(glm::rotate(glm::mat4(1.0f), radians, glm::vec3(0, 1, 0)));

	for (Sphere &s : spheres)
		s.center = s.center * rotationMatrix;

	updateBBox();
}

void SphereMesh::printf() const{
	std::cout
			<< "Sphere mesh: "
			<< spheres.size() << " spheres, "
			<< capsuloids.size() << " caps, "
			<< prysmoids.size() << " prys"
			<< std::endl;
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

	if (content.size() < 3)
	{
		std::cout << "Invalid file format!" << std::endl;
		return false;
	}

	std::string structCounter = content[1];
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
		Sphere s = extractSphereFromString(content[2 + i]);
		spheres.push_back(s);
	}

	for (int i = 0; i < nPrysmoids; i++)
	{
		Prysmoid p = extractPrysmoidFromString(content[2 + nSpheres + i]);
		prysmoids.push_back(p);
	}

	for (int i = 0; i < nCapsuloids; i++)
	{
		Capsuloid c = extractCapsuloidFromString(content[2 + nSpheres + nPrysmoids + i]);
		capsuloids.push_back(c);
	}

	updateBBox();

	removeDegenerateElements();
	removeRedundantElements();

	return true;
}

int Plane::cellTest(const std::vector<glm::vec3> &cell) const
{
	int behind = 0;
	for (auto& pos : cell)
		if(isBehind(pos))
			++behind;

	if (behind == 8)
		return -1; // totally behind SHOULD NOT HAPPEN

	if (behind > 0)
		return 0; // mixed

	return 1; // totally in front
}

std::string Plane::serialize() const
{
	std::string res;

	res += std::to_string(n.x) + " " + std::to_string(n.y) + " " + std::to_string(n.z) + " ";
	res += std::to_string(k);

	return res;
}

JoinedQuadrilateral::JoinedQuadrilateral(const Plane& p, const Sphere &a, const Sphere &b, const Sphere &c, const Sphere &d,
	int i, int j, int k, int l, int dir)
{
	midPlane = p;

	spheres[0] = a;
	spheres[1] = b;
	spheres[2] = c;
	spheres[3] = d;

	indices[0] = i;
	indices[1] = j;
	indices[2] = k;
	indices[3] = l;

	sortIndices();
	upperPlane = computeUpperPlane(dir);
	// TODO: Check this function
	computeSphereRadii();
}

void JoinedQuadrilateral::sortIndices()
{
	if (indices[0] < indices[2])
	{
		std::swap(indices[0], indices[2]);
		std::swap(spheres[0], spheres[2]);
	}

	if (indices[1] < indices[3])
	{
		std::swap(indices[1], indices[3]);
		std::swap(spheres[1], spheres[3]);
	}

	if (indices[0] < indices[1])
	{
		std::swap(indices[0], indices[1]);
		std::swap(indices[2], indices[3]);
		std::swap(spheres[0], spheres[1]);
		std::swap(spheres[2], spheres[3]);
	}
}

Plane JoinedQuadrilateral::computeUpperPlane(int direction) const
{
	glm::vec3 a = spheres[0].center;
	glm::vec3 b = spheres[1].center;
	glm::vec3 c = spheres[2].center;
	glm::vec3 d = spheres[3].center;

	auto sign = static_cast<float>(direction);

	glm::vec3 n = sign * midPlane.n;
	glm::vec3 startN = n;

	for (int i = 0; i < 1000; i++){
		a = spheres[0].center + n * spheres[0].radius;
		b = spheres[1].center + n * spheres[1].radius;
		c = spheres[2].center + n * spheres[2].radius;
		d = spheres[3].center + n * spheres[3].radius;

		glm::vec3 new_n = glm::normalize(sign * glm::cross(a - c, b - d));
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

// TODO: Check this function
void JoinedQuadrilateral::computeSphereRadii()
{
	spheres[0].radius = glm::distance(spheres[0].center, midPlane.project(spheres[0].center));
	spheres[1].radius = glm::distance(spheres[1].center, midPlane.project(spheres[1].center));
	spheres[2].radius = glm::distance(spheres[2].center, midPlane.project(spheres[2].center));
	spheres[3].radius = glm::distance(spheres[3].center, midPlane.project(spheres[3].center));
}

bool SphereMesh::loadFromText(const char *text)
{
	std::stringstream inputFile(text);

	std::string line;
	std::vector<std::string> content;

	while (std::getline(inputFile, line))
		content.push_back(line);

	if (content.size() < 3)
	{
		std::cout << "Invalid file format!" << std::endl;
		return false;
	}

	std::string structCounter = content[1];
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
		Sphere s = extractSphereFromString(content[2 + i]);
		spheres.push_back(s);
	}

	for (int i = 0; i < nPrysmoids; i++)
	{
		Prysmoid p = extractPrysmoidFromString(content[2 + nSpheres + i]);
		prysmoids.push_back(p);
	}

	for (int i = 0; i < nCapsuloids; i++)
	{
		Capsuloid c = extractCapsuloidFromString(content[2 + nSpheres + nPrysmoids + i]);
		capsuloids.push_back(c);
	}

	updateBBox();

	removeDegenerateElements();
	removeRedundantElements();

	return true;
}

Sphere SphereMesh::extractSphereFromString(const std::string &sphereString)
{
	std::string s = sphereString;

	Sphere sphere{};

	const std::string delimiter = " ";
	std::string token = s.substr(0, s.find(delimiter));
	sphere.center.x = std::stof(token);
	s.erase(0, s.find(delimiter) + delimiter.length());

	token = s.substr(0, s.find(delimiter));
	sphere.center.y = std::stof(token);
	s.erase(0, s.find(delimiter) + delimiter.length());

	token = s.substr(0, s.find(delimiter));
	sphere.center.z = std::stof(token);
	s.erase(0, s.find(delimiter) + delimiter.length());

	token = s.substr(0, s.find(delimiter));
	sphere.radius = std::stof(token);

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

	Capsuloid capsuloid{};
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

	Prysmoid prysmoid{};

	prysmoid.indices[0] = prysmoidA;
	prysmoid.indices[1] = prysmoidB;
	prysmoid.indices[2] = prysmoidC;

	return prysmoid;
}

bool SphereMesh::saveToFile (const char *path, const char *name) const
{
	std::ostringstream fileContent;

	fileContent << "Sphere Mesh 2.0" << std::endl;
	fileContent << spheres.size();
	fileContent << " " << prysmoids.size();
	fileContent << " " << capsuloids.size() << std::endl;

	for (auto& s : spheres)
		fileContent << s.center.x << " " << s.center.y << " " << s.center.z << " " << s.radius << std::endl;

	for (auto& p : prysmoids)
		fileContent << p.indices[0] << " " << p.indices[1] << " " << p.indices[2] << std::endl;

	for (auto& c : capsuloids)
		fileContent << c.indices[0] << " " << c.indices[1] << std::endl;

	namespace fs = std::__fs::filesystem;
	fs::path cwd = fs::current_path();

	std::string separator = std::string(1, fs::path::preferred_separator);

	const std::string& filePath = name;
	std::string folderPath = "." + separator;
	if (path != std::string("."))
		folderPath = path;

	std::ofstream fout(folderPath + filePath);
	std::cout << "File location: " << folderPath + filePath << std::endl;
	fout << fileContent.str();
	fout.close();

	return false;
}

std::vector<std::pair<Prysmoid, Prysmoid>> SphereMesh::findPrysmoidPairs()
{

}

int SphereMesh::intersectedSphereAlongRay(const glm::vec3& rayOrigin, const glm::vec3& rayDir, glm::vec3& hitPos) {
    int closest = -1;
    float closestDist = FLT_MAX;

    for (int i = 0; i < spheres.size(); ++i) {
        const Sphere& sphere = spheres[i];
        float t;
        glm::vec3 currHitPos;
        if (raySphereIntersection(rayOrigin, rayDir, sphere.center, sphere.radius, t, currHitPos))
        {
            float currDist = glm::distance2(sphere.center, rayOrigin);
            if (currDist < closestDist)
            {
                hitPos = currHitPos;
                closestDist = currDist;
                closest = i;
            }
        }
    }

    return closest;
}

bool SphereMesh::raySphereIntersection(const glm::vec3& rayOrigin, const glm::vec3& rayDir,
                                       const glm::vec3& sphereCenter, float sphereRadius, float& t, glm::vec3& hitPos) {
    glm::vec3 L = sphereCenter - rayOrigin;
    float tca = glm::dot(L, rayDir);

    float d2 = glm::dot(L, L) - tca * tca;
    float radius2 = sphereRadius * sphereRadius;
    if (d2 > radius2)
        return false;

    float thc = sqrt(radius2 - d2);
    float t0 = tca - thc;
    float t1 = tca + thc;

    if (t1 < 0)
        return false;

    t = (t0 < 0) ? t1 : t0;
    hitPos = rayOrigin + rayDir * t;

    return true;
}

void SphereMesh::duplicateSphere(int i)
{
    if (i >= spheres.size() || i < 0)
        return;

    Sphere newS = spheres[i];
    newS.center += 1.0f;
    spheres.push_back(newS);
}

void SphereMesh::removeSphere(int i)
{
    if (i >= spheres.size() || i < 0)
        return;

    for (int idx = 0; idx < capsuloids.size();)
    {
        int x = capsuloids[idx].indices[0];
        int y = capsuloids[idx].indices[1];
        if (x == i || y == i)
            capsuloids.erase(capsuloids.begin() + idx);
        else
        {
            if (x > i) x--;
            if (y > i) y--;
            capsuloids[idx] = {x, y};
            idx++;
        }
    }

    for (int idx = 0; idx < prysmoids.size();)
    {
        int x = prysmoids[idx].indices[0];
        int y = prysmoids[idx].indices[1];
        int z = prysmoids[idx].indices[2];

        if (x == i || y == i || z == i)
        {
            std::vector<int> remaining;
            if (x != i) remaining.push_back(x);
            if (y != i) remaining.push_back(y);
            if (z != i) remaining.push_back(z);

            if (remaining.size() == 2)
                addCapsuloid(remaining[0], remaining[1]);

            prysmoids.erase(prysmoids.begin() + idx);
        }
        else
        {
            if (x > i) x--;
            if (y > i) y--;
            if (z > i) z--;
            prysmoids[idx] = {x, y, z};
            idx++;
        }
    }

    spheres.erase(spheres.begin() + i);
}

void SphereMesh::addCapsuloid(int i, int j)
{
    if (i >= spheres.size() || i < 0 || j >= spheres.size() || j < 0 || i == j)
        return;

    if (i > j)
        std::swap(i, j);

    for (int idx = 0; idx < capsuloids.size(); idx++)
    {
        int x = capsuloids[idx].indices[0];
        int y = capsuloids[idx].indices[1];
        if (x == i && y == j)
            return;
    }

    capsuloids.push_back({i, j});
}

void SphereMesh::addPrysmoid(int i, int j, int k)
{
    if (i >= spheres.size() || i < 0 || j >= spheres.size() || j < 0
        || k >= spheres.size() || k < 0 || i == k || i == j || j == k)
        return;

    std::vector<int> indices;
    indices.push_back(i);
    indices.push_back(j);
    indices.push_back(k);
    std::sort(indices.begin(), indices.end());

    for (int idx = 0; idx < prysmoids.size(); idx++)
    {
        int x = prysmoids[idx].indices[0];
        int y = prysmoids[idx].indices[1];
        int z = prysmoids[idx].indices[2];
        if (x == indices[0] && y == indices[1] && z == indices[2])
            return;
    }

    prysmoids.push_back({indices[0], indices[1], indices[2]});
}

void SphereMesh::removeLink(int i, int j)
{
    if (i >= spheres.size() || i < 0 || j >= spheres.size() || j < 0 || i == j)
        return;

    if (i > j)
        std::swap(i, j);

    for (int idx = 0; idx < capsuloids.size(); idx++)
    {
        int x = capsuloids[idx].indices[0];
        int y = capsuloids[idx].indices[1];
        if ((x == i && y == j) || (x == j && y == i))
        {
            capsuloids.erase(capsuloids.begin() + idx);
            return;
        }
    }
}

void SphereMesh::removeLink(int i, int j, int k)
{
    if (i >= spheres.size() || i < 0 || j >= spheres.size() || j < 0 || k >= spheres.size() || k < 0)
        return;

    if (i == j || i == k || j == k)
        return;

    std::vector<int> indices = {i, j, k};
    std::sort(indices.begin(), indices.end());

    for (int idx = 0; idx < prysmoids.size(); idx++)
    {
        int x = prysmoids[idx].indices[0];
        int y = prysmoids[idx].indices[1];
        int z = prysmoids[idx].indices[2];

        if (x == indices[0] && y == indices[1] && z == indices[2])
        {
            prysmoids.erase(prysmoids.begin() + idx);
            addCapsuloid(x, y);
            addCapsuloid(x, z);
            addCapsuloid(z, y);
            return;
        }
    }
}

// TODO: Implement the SphereMeshBlendShape class
bool SphereMeshBlendShape::loadFromFile (const char *filename)
{
	SphereMesh sm;

	/*if (!sm.loadFromFile(filename))
		return false;*/

	shapes.push_back(sm);

	return true;
}
