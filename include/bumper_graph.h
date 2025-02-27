#pragma once

#include <unordered_set>
#include <variant>

#include "sphere_mesh.h"

namespace SM::Graph
{
	struct Cone {
		int capsuloidIndex {-1};

		glm::vec3 apex {};
		glm::vec3 axis {};
		float halfAngle {};

		Cone() = default;
		Cone(const glm::vec3& apex, const glm::vec3& axis, const float halfAngle)
			: apex(apex), axis(normalize(axis)), halfAngle(halfAngle)
		{}

		[[nodiscard]] bool contains(const glm::vec3& point) const {
			const glm::vec3 v = point - apex;

			const float dotWithAxis = dot(v, axis);

			if (dotWithAxis < 0.0f)
				return false;

			const float dist = length(v);
			if (dist < 1e-6f)
				return true;

			const float cosTheta = dotWithAxis / dist;

			return cosTheta >= std::cos(halfAngle);
		}

		void translate(const glm::vec3& t) {
			apex += t;
		}

		void rotate(const glm::mat3& rot) {
			apex = rot * apex;
			axis = normalize(rot * axis);
		}

		void scale(float s) {
			apex *= s;
		}
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

		void rotateAround(const glm::vec3& axis, float angle);

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
		[[nodiscard]] static Plane computeUpperPlane(const Sphere &sa, const Sphere &sb, const Sphere &sc, int direction);
	};

	class BumperCapsuloid {
	public:
		float k{};

		std::vector<Plane> neib;
		int sphereIndex[2]{ -1, -1 };

		glm::vec3 dNorm{};
		float dLength{};
		float dR{};

		std::vector<int> neibBumpers;
		Cone cones[2] {};

		void translate(const glm::vec3& t);
		void scale(float r);
		void rotate(const glm::mat3& rot);

		void rotateAround(const glm::vec3& axis, float angle);

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

		void translate(const glm::vec3& t);
		void scale(float r);
		void rotate(const glm::mat3& rot);

		void rotateAround(const glm::vec3& axis, float angle);

		BumperSphere() = default;
		explicit BumperSphere(const int sphereIndex) : sphereIndex(sphereIndex) {}
	};

	class BumperQuad
	{
	public:
		Plane upperPlane {};
		Plane midPlane {};
		Plane sidePlanes[4] {};

		int sphereIndex[4] {-1, -1, -1, -1};
		int neibSide[4] {-1, -1, -1, -1};
		int neibOpp{};

		void translate(const glm::vec3& t);
		void scale(float r);
		void rotate(const glm::mat3& rot);

		void rotateAround(const glm::vec3& axis, float angle);

		BumperQuad() = default;
		explicit BumperQuad(const FourSpheres& fs);

	private:
		void getSidePlanes(const BumperPrysmoid& bp, const BumperPrysmoid& bp1);
	};

	class Bumper {
	public:
		std::variant<BumperPrysmoid, BumperCapsuloid, BumperSphere, BumperQuad> bumper;

		enum {
			SPHERE,
			CAPSULOID,
			PRYSMOID,
			QUAD
		} shapeType {};

		[[nodiscard]] bool hasAParent() const;
		[[nodiscard]] std::string serialize() const;

		void translate(const glm::vec3& t);
		void scale(float r);
		void rotate(const glm::mat3& rot);

		void rotateAround(const glm::vec3& axis, float angle);

		bool operator == (const Bumper &) const;

		Bumper() : bumper(BumperSphere()) {}
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
		bool pushOutsideBruteForce(glm::vec3& p, glm::vec3& n, const int& bumperIndex)  const;
		bool projectOn(glm::vec3& p, glm::vec3& n, int& bumperIndex) const;

		void animateRig(int iterations);

	private:
		static float PUSH_EPSILON;

		struct Transform
		{
			glm::vec3 axis;
			float angle;
			glm::vec3 root;
		};

		std::unordered_map<int, int> rig;
		int safeLeft = 5, safeRight = 3, cancer = 6;
		std::unordered_map<int, Transform> transformMapper;

		[[nodiscard]] std::vector<std::vector<bool>> getBumperAdjMatrix() const;
		void fillBranchWithTransform(int value, const std::vector<std::vector<bool>> &adjMatrix, std::queue<int>& q,
		                             std::unordered_set<int>& visited);
		void floodfill();

		float projectOnBruteForce(glm::vec3& p, glm::vec3& n, int& bumperIndex)  const;

		bool pushOutsideSphere(glm::vec3& p, glm::vec3& n, const int& bumperIndex)  const;
		bool pushOutsideCapsuloid(glm::vec3& p, glm::vec3& n, int& bumperIndex)  const;
		bool pushOutsidePrysmoid(glm::vec3& p, glm::vec3& n, int& bumperIndex)  const;
		bool pushOutsideQuad(glm::vec3& p, glm::vec3& n, int& bumperIndex)  const;

		[[nodiscard]] Sphere getInterpolatedSphere(const BumperCapsuloid& bc, float t) const;
		[[nodiscard]] float closestSphereOn(const glm::vec3& p, const BumperCapsuloid &bc) const;

		[[nodiscard]] bool isPointOverPrysmoid(int bumperIndex, const glm::vec3& p) const;
		[[nodiscard]] bool isPointOverQuad(int bumperIndex, const glm::vec3& p) const;

		float signedDistanceFromSphere(int bumperIndex, const glm::vec3& p, glm::vec3& closestPos, glm::vec3&
		closestNorm) const;
		float signedDistanceFromCapsuloid(int bumperIndex, const glm::vec3& p, glm::vec3& closestPos, glm::vec3&
		closestNorm) const;
		float signedDistanceFromPrysmoid(int bumperIndex, const glm::vec3& p, glm::vec3& closestPos, glm::vec3&
		closestNorm) const;
		float signedDistanceFromQuad(int bumperIndex, const glm::vec3& p, glm::vec3& closestPos, glm::vec3&
		closestNorm) const;

		[[nodiscard]] std::pair<int, float> sampleSignedDistanceFromSphere(int bumperIndex, const glm::vec3& p) const;
		[[nodiscard]] std::pair<int, float> sampleSignedDistanceFromCapsuloid(int bumperIndex, const glm::vec3& p) const;
		[[nodiscard]] std::pair<int, float> sampleSignedDistanceFromPrysmoid(int bumperIndex, const glm::vec3& p) const;
		[[nodiscard]] std::pair<int, float> sampleSignedDistanceFromQuad(int bumperIndex, const glm::vec3& p) const;

		[[nodiscard]] std::pair<int, float> signedDistanceFromBumper(int i, const glm::vec3& p) const;
		void sampleDistanceFromBumperWithPosition(const glm::vec<3, float> &p, int proposedIndex, glm::vec3 &closestPos) const;

		void initializeBumperSpheres();
		void initializeBumperCapsuloids(const SphereMesh &sm, std::vector<std::vector<int>> &capsuloidAdj);
		void initializeBumperPrysmoids(const SphereMesh &sm, std::vector<std::vector<int>> &capsuloidAdj);
		void initializeBumperQuads(const SphereMesh &sm, std::vector<std::vector<int>> &capsuloidAdj);
		void initializeBumperNodes(const SphereMesh &sm);

		void sortByType();
		[[nodiscard]] FourSpheres computeFourSpheresFrom(const Quadrilateral& quad, int direction) const;
	};

	struct SubMesh {
		std::string name;
		std::vector<glm::vec3> vertices;
		std::vector<std::array<int, 3>> faces;

		void addQuad(int a, int b, int c, int d);
		void addTriangle(int a, int b, int c);
	};

	struct hash_pair {
		template <typename T1, typename T2>
		std::size_t operator()(const std::pair<T1,T2>& p) const {
			return std::hash<T1>{}(p.first) ^ std::hash<T2>{}(p.second) << 1;
		}
	};

	struct VertexColor {
		glm::vec3 pos;
		unsigned char r, g, b;
	};

	class BumperGraphMesh {
	public:
		std::vector<SubMesh> submeshes;

		explicit BumperGraphMesh(const BumperGraph &bg);
		void exportObj(const std::string &filename) const;

		[[nodiscard]] static SubMesh createSphereSubMesh(const Sphere &s, int resolution);
		[[nodiscard]] static SubMesh createCapsuleSubMesh(const Sphere &s0, const Sphere &s1, int resolution);
		[[nodiscard]] static SubMesh createClosedCapsuleSubMesh(const Sphere &s0, const Sphere &s1, int resolution);
		[[nodiscard]] static SubMesh createPlaneSubMesh(const Plane &p, float width, float height);
		[[nodiscard]] static SubMesh createSubMeshFromBumper(const Bumper &b, const std::vector<Sphere>& spheres);
		[[nodiscard]] static SubMesh createPrysmoidSubMesh(const BumperPrysmoid &prys, const std::vector<Sphere> &spheres);
		[[nodiscard]] static SubMesh createExtraLateralCutPlaneForPrysmoidSideWithElongation(
																		const BumperPrysmoid &prys,
																		const std::vector<Sphere> &spheres,
																		int side,
																		float topOffset,
																		float sideElongation);

		void exportPly(const std::string &filename) const;

	private:
		static void appendSubMesh(SubMesh &dest, const SubMesh &src);
	};
}
