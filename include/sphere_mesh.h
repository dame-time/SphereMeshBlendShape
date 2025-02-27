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

#pragma once

#include <glm/glm.hpp>

#include <vector>
#include <string>
#include <istream>
#include <iostream>
#include <unordered_map>

#define MAX_GRID_BUMPERS 5

namespace SM {
	struct Plane
	{
        bool valid { true };
		glm::vec3 n {}; // unitary vector
		float k {};

		Plane() = default;
        Plane(const glm::vec3& norm, const float _k) : n(norm), k(_k) {}
		Plane(const glm::vec3& dir, const glm::vec3& p) : n(glm::normalize(dir)) { setPoint(p); }

		[[nodiscard]] bool isBehind(const glm::vec3& p) const { return glm::dot(n, p) < k; }
		[[nodiscard]] int cellTest(const std::vector<glm::vec3>& cell) const;
		void setPoint(const glm::vec3& p) { k = glm::dot(n, p); }

		void translate (const glm::vec3& t) { k += glm::dot(n, t); };
		void scale (const float s) { k *= s; }
		void rotate (const glm::mat3& rot) { n = rot * n; }

		[[nodiscard]] float distance(const glm::vec3& p) const { return glm::dot(n, p) - k; }
		[[nodiscard]] glm::vec3 project(const glm::vec3& p) const { return p - distance(p) * n; }

		[[nodiscard]] std::string serialize() const;

		bool operator == (const Plane& other) const;

		void flip () { n = -n; k = -k; }
	};

	class Sphere {
	public:
		glm::vec3 center;
		float radius;

		Sphere() = default;
		Sphere(const glm::vec3& center, float radius) : center(center), radius(radius) {}
	};

	struct Triangle {
		glm::vec3 n;
		glm::vec3 a, b, c;
	};

	struct Singleton {
		int index;
	};

	struct Capsuloid {
		int indices[2];
	};

	struct Prysmoid {
		int indices[3];
	};

	struct Quadrilateral {
		int indices[4];
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
			return {
					minCorner.x + static_cast<float>(rand()) / RAND_MAX * (maxCorner.x - minCorner.x),
					minCorner.y + static_cast<float>(rand()) / RAND_MAX * (maxCorner.y - minCorner.y),
					minCorner.z + static_cast<float>(rand()) / RAND_MAX * (maxCorner.z - minCorner.z)
			};
		}

		[[nodiscard]] glm::vec3 BDD() const
		{
			return maxCorner - minCorner;
		}
	};

    class FourSpheres
	{
		public:
			int indices[4] {};
            int prysmoidIndices[2] {};
			Sphere spheres[4] {};
            bool failed {};

            bool degenerate {};
            float error {};

			Plane midPlane {};
			Plane upperPlane {};

            FourSpheres() = default;
            explicit FourSpheres(bool isFailed);
            FourSpheres(const Sphere& a, const Sphere& b, const Sphere& c, const Sphere& d);

    		[[nodiscard]] Plane computeCommonTangentPlane(int direction) const;
    		[[nodiscard]] Plane computeCommonMidPlane(int direction) const;

            void computeError(const FourSpheres& other);

		private:
			void sortIndices();
	};

    struct CompareByError {
        bool operator()(const FourSpheres& a, const FourSpheres& b) const {
            return a.error > b.error;
        }
    };

	class SphereMesh {
	public:
		AABB bbox;

		std::vector<Sphere> spheres;

		std::vector<Singleton> singletons;
		std::vector<Capsuloid> capsuloids;
		std::vector<Prysmoid> prysmoids;
		std::vector<Quadrilateral> quadrilaterals;

		bool loadFromText(const char* text);
		bool loadFromFile(const char* text);

		bool saveToFile(const char* path, const char* name = "outSM.sm") const;

        [[nodiscard]] std::priority_queue<FourSpheres, std::vector<FourSpheres>, CompareByError> findJoinableQuadrilaterals() const;
        [[nodiscard]] FourSpheres joinPrysmoids(int p1, int p2) const;
        [[nodiscard]] FourSpheres joinPrysmoids(int p1, int p2, int a, int b, int c, int d) const;
        bool generateQuadrilateral(float threshold);
        [[nodiscard]] bool isSphereInQuadrilateral(int sphereIdx) const;
        std::vector<int> getSymmetricalSpheres(int sphereIdx);

		int intersectedSphereAlongRay(const glm::vec3& rayOrigin, const glm::vec3& rayDir, glm::vec3& hitPos) const;

		static bool raySphereIntersection(const glm::vec3& rayOrigin, const glm::vec3& rayDir,
		                                  const glm::vec3& sphereCenter, float sphereRadius, float& t, glm::vec3& hitPos);

		void duplicateSphere(int i);
		void removeSphere(int i);
		void addCapsuloid(int i, int j);
		void addPrysmoid(int i, int j, int k);
		void removeLink(int i, int j);
		void removeLink(int i, int j, int k);

		void translateScale( glm::vec3 t , float s);
		void scaleTranslate( glm::vec3 t , float s);
		void rotoTranslate(const glm::mat3 &rot, const glm::vec3& trasl);

		void translate (const glm::vec3& t);
		void scale (float s);
		void rotateY (int angle);

        void printf() const;

	private:
        std::unordered_map<int, std::vector<int>> quadsSpheres;

		void updateBBox();

		void removeDegeneratePrysmoids();
		void removeDegenerateCapsuloids();

		void removeDegenerateElements();
		void removeRedundantElements();

        bool checkPrysmEdgesAgainstQuadDiagonal(int s1, int s2) const;

		static Sphere extractSphereFromString(const std::string &sphereString);
		static Capsuloid extractCapsuloidFromString(const std::string &s);
		static Prysmoid extractPrysmoidFromString(const std::string &s);
		Quadrilateral extractQuadFromString(const std::string &s);

		void computeFourSpheresFrom(const Quadrilateral& quad);
	};

	// TODO: implement this class
	class SphereMeshBlendShape{
	public:
		SphereMesh extract(float t); // t is in 0 to n-1
		bool loadFromFile(const char* filename);
		std::vector<SphereMesh> shapes; // n elements
	private:
	};

	inline std::istream& operator>>(std::istream& is, Plane& plane)
	{
		is >> plane.n.x >> plane.n.y >> plane.n.z >> plane.k;
		return is;
	}

	inline Plane commonTangentPlane(const Sphere& s0, const Sphere& s1, const Sphere& s2, const Sphere& s3, int direction)
	{
		glm::vec3 a = s0.center;
		glm::vec3 b = s1.center;
		glm::vec3 c = s2.center;
		glm::vec3 d = s3.center;

		auto sign = static_cast<float>(direction);

		auto n = glm::vec3(0);
		glm::vec3 startN = n;

		Plane p;

		for (int i = 0; i < 1000; i++){
			a = s0.center + n * s0.radius;
			b = s1.center + n * s1.radius;
			c = s2.center + n * s2.radius;
			d = s3.center + n * s3.radius;

			glm::vec3 new_n = glm::normalize(sign * glm::cross(a - c, b - d));
			if (glm::dot(n, new_n) >= 0.999f)
			{
				if (glm::dot(startN, new_n) < 0)
				{
					std::cerr << "Degenerate Prysmoid detected (Normal flipped)" << std::endl;
					p.valid = false;
					return p;
				}

				p = {new_n, (a + b + c + d) / 4.0f};
				return p;
			}
			n = new_n;
		}

		std::cerr << "Degenerate Prysmoid detected (Normal diverged)" << std::endl;
		p.valid = false;
		return p;
	}
}
