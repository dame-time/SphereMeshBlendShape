#include "SphereMeshBlendShape.hpp"

#include <fstream>
#include <iostream>
#include <string>

bool SphereMesh::loadFromFile (const char *filename)
{
	std::ifstream inputFile(filename);
	
	if (!inputFile.is_open()) {
		std::cerr << "Error opening the file!" << std::endl;
		return false;
	}
	
	std::string line;
	std::vector<std::string> content;
	
	while (std::getline(inputFile, line))
		content.push_back(line);
	
	inputFile.close();
	
	if (content.size() < 5)
	{
		std::cerr << "Invalid file format!" << std::endl;
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
	
	return true;
}

Sphere SphereMesh::extractSphereFromString(const std::string &sphereString)
{
	std::string s = sphereString;
	
	Sphere sphere;
	
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
	
	Capsuloid capsuloid;
	
	std::string delimiter = " ";
	std::string token = s.substr(0, s.find(delimiter));
	capsuloid.a = std::stoi(token);
	s.erase(0, s.find(delimiter) + delimiter.length());
	
	token = s.substr(0, s.find(delimiter));
	capsuloid.b = std::stoi(token);
	
	return capsuloid;
}

Prysmoid SphereMesh::extractPrysmoidFromString(const std::string &prysmoidString)
{
	std::string s = prysmoidString;
	
	Prysmoid prysmoid;
	
	std::string delimiter = " ";
	std::string token = s.substr(0, s.find(delimiter));
	prysmoid.a = std::stoi(token);
	s.erase(0, s.find(delimiter) + delimiter.length());
	
	token = s.substr(0, s.find(delimiter));
	prysmoid.b = std::stoi(token);
	s.erase(0, s.find(delimiter) + delimiter.length());
	
	token = s.substr(0, s.find(delimiter));
	prysmoid.c = std::stoi(token);
	
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