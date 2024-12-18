#pragma once

#include <variant>

#include "sphere_mesh.h"

namespace SM::Graph
{
	class Cone
	{

	};

	class BumperPrysmoid {
	public:
		Plane upperPlane;

		Plane m1{};
		Plane m2{};
		Plane m0{};

		Plane midPlane {};

		int sphereIndex[3] {-1, -1, -1};

		int neibSide[3] {-1, -1, -1};
		int neibOpp{};

		bool operator == (const BumperPrysmoid& bp) const;

		void translate(const glm::vec3& t);
		void scale(float r);
		void rotate(const glm::mat3& rot);

		BumperPrysmoid() = default;
		BumperPrysmoid(
				int indexA,
				int indexB,
				int indexC,
				const Sphere& v0,
				const Sphere& v1,
				const Sphere& v2,
				int direction);

	private:
		[[nodiscard]] static Plane computeUpperPlane(const Sphere &a, const Sphere &b, const Sphere &c, int direction);
	};

	class BumperCapsuloid {
	public:
		float k{};

		std::vector<Plane> neib;
		int sphereIndex[2]{ -1, -1 };

		glm::vec3 dNorm{};
		float dLength{};
		float dR{};

		int neibExtreme[2] {-1, -1};

		bool hasFather = true;

		bool operator == (const BumperCapsuloid& bc) const;

		BumperCapsuloid() = default;
		BumperCapsuloid(int indexA, int indexB, const Sphere& a, const Sphere& b);
	};

	class BumperSphere {
	public:
		int sphereIndex {-1};

		std::vector<Cone> cones;

		BumperSphere() = default;
		explicit BumperSphere(const int sphereIndex) : sphereIndex(sphereIndex) {};
	};

	class Bumper {
	public:
		std::variant<BumperPrysmoid, BumperCapsuloid, BumperSphere> bumper;

		enum {
			SPHERE,
			CAPSULOID,
			PRYSMOID,
		} shapeType {};

		[[nodiscard]] bool hasAParent() const;
		[[nodiscard]] std::string serialize() const;

		bool operator == (const Bumper &) const;

		Bumper() : bumper(BumperSphere()) {};
	};

	class BumperGraph {
	public:
		std::vector<Bumper> bumper;
		std::vector<Sphere> sphere;

		void constructFrom (const SphereMesh& sm);

		void translate (const glm::vec3& t);
		void scale (float s);
		void rotateY (int angle);

		bool pushOutside(glm::vec3& p, glm::vec3& n, int& bumperIndex)  const;
		bool pushOutsideBruteForce(glm::vec3& p, glm::vec3& n, int& bumperIndex)  const;
		bool projectOn(glm::vec3& p, glm::vec3& n, int& bumperIndex) const;

	private:
		static float PUSH_EPSILON;

		float projectOnBruteForce(glm::vec3& p, glm::vec3& n, int& bumperIndex)  const;

		bool pushOutsideSphere(glm::vec3& p, glm::vec3& n, int& bumperIndex)  const;
		bool pushOutsideCapsuloid(glm::vec3& p, glm::vec3& n, int& bumperIndex)  const;
		bool pushOutsidePrysmoid(glm::vec3& p, glm::vec3& n, int& bumperIndex)  const;

		[[nodiscard]] Sphere getInterpolatedSphere(const BumperCapsuloid& bc, float t) const;
		[[nodiscard]] float closestSphereOn(const glm::vec3& p, const BumperCapsuloid &bc) const;

		[[nodiscard]] bool isPointOverPrysmoid(int bumperIndex, const glm::vec3& p) const;

		float signedDistanceFromSphere(int bumperIndex, const glm::vec3& p, glm::vec3& closestPos, glm::vec3&
		closestNorm) const;
		float signedDistanceFromCapsuloid(int bumperIndex, const glm::vec3& p, glm::vec3& closestPos, glm::vec3&
		closestNorm) const;
		float signedDistanceFromPrysmoid(int bumperIndex, const glm::vec3& p, glm::vec3& closestPos, glm::vec3&
		closestNorm) const;

		[[nodiscard]] std::pair<int, float> sampleSignedDistanceFromSphere(int bumperIndex, const glm::vec3& p) const;
		[[nodiscard]] std::pair<int, float> sampleSignedDistanceFromCapsuloid(int bumperIndex, const glm::vec3& p) const;
		[[nodiscard]] std::pair<int, float> sampleSignedDistanceFromPrysmoid(int bumperIndex, const glm::vec3& p) const;

		[[nodiscard]] std::pair<int, float> signedDistanceFromBumper(int i, const glm::vec3& p) const;
		void sampleDistanceFromBumperWithPosition(const glm::vec<3, float> &p, int proposedIndex, glm::vec3 &closestPos) const;

		void initializeBumperSpheres();
		void initializeBumperCapsuloids(const SphereMesh &sm, std::vector<std::vector<int>> &capsuloidAdj);
		void initializeBumperPrysmoids(const SphereMesh &sm, std::vector<std::vector<int>> &capsuloidAdj);
		void initializeBumperNodes(const SphereMesh &sm);

		void sortByType();
	};
}
