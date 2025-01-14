#pragma once

#include <variant>
#include <set>

#include "sphere_mesh.h"

namespace SM::Grid
{
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

		BumperSphere() = default;
		explicit BumperSphere(const int sphereIndex) : sphereIndex(sphereIndex) {};
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

			BumperQuad() = default;
			BumperQuad(const FourSpheres& fs, int direction);

		private:
			std::vector<Plane> getSidePlanes(const BumperPrysmoid& bp, const BumperPrysmoid& bp1);
	};

	struct CompositeBumper
	{
		int bumperIndices[MAX_GRID_BUMPERS] {};
		char flags[MAX_GRID_BUMPERS] {};

		CompositeBumper()
		{
			for (char& flag : flags)
				flag = 0;

			for (int& bumperIdx : bumperIndices)
				bumperIdx = -1;
		}

		bool operator == (const CompositeBumper& bn) const;

		[[nodiscard]] std::string serialize() const;
	};

	class Bumper {
	public:
		std::variant<BumperPrysmoid, BumperCapsuloid, BumperSphere, BumperQuad, CompositeBumper> bumper;

		enum {
			SPHERE,
			CAPSULOID,
			PRYSMOID,
			QUAD,
			COMPOSITE
		} shapeType {};

		[[nodiscard]] bool hasAParent() const;
		[[nodiscard]] std::string serialize() const;

		bool operator == (const Bumper &) const;

		Bumper() : bumper(BumperSphere()) {};
	};

	struct GridSample
	{
		glm::vec3 position;
		glm::vec3 dest;
		int index;
	};

	class BumperGrid;
	class SpatialGrid;

	class WIPGrid
	{
	public:
		WIPGrid() = default;

		void format(const AABB& boundingBox, int requiredNumCells);
		void populate(const std::vector<GridSample>& samples);

		void bake(BumperGrid& bg) const;

		[[nodiscard]] int getIndex(const glm::vec3& position) const;
		[[nodiscard]] glm::ivec3 getCellCoords(const glm::vec3& position) const;

	private:
		std::vector<std::vector<int>> cells;
		std::vector<std::set<int>> cellsSet;

		glm::vec3 offset{};
		float scale{};
		glm::ivec3 resolution{};

		static int howManyCellsForScale(const AABB& bbox, float scale);

		void setSize(const AABB& boundingBox, float requiredScale);
		[[nodiscard]] int numCells() const;
	};

	class SpatialGrid {
	friend class WIPGrid;
	friend class BumperGrid;
	public:
		SpatialGrid() = default;

		[[nodiscard]] std::string serialize() const;
		void deserialize(const std::string& data);

		[[nodiscard]] int indicesAt(const glm::vec3& position) const;
		[[nodiscard]] int getIndex(const glm::vec3& position) const;
		[[nodiscard]] glm::ivec3 getCellCoords(const glm::vec3& position) const;
		[[nodiscard]] std::vector<glm::vec3> getCellCorners(int i) const;

		void translateGrid(const glm::vec3& t);
		void scaleGrid(float s);
		void rotateYGrid(int angle);

#ifdef DEBUG
		void drawGrid() const;
		void drawGridWithSpheres(const BumperGrid &bg) const;
		void drawGridWithCapsules(const BumperGrid &bg) const;
		void drawGridWithPrysmoids(const BumperGrid &bg) const;
		void drawGridWithDensityColors(const BumperGrid &bg, bool fillCubes) const;
		void drawGridWithBumpers(const BumperGrid &bg, int n) const;
		void printHistogram(const BumperGrid &bg) const;
		[[nodiscard]] std::vector<std::vector<int>> getCellsWithMaxElements(const BumperGrid &bg) const;
#endif

	private:
		std::vector<int> cells;

		glm::vec3 offset{};
		float scale{};
		glm::ivec3 resolution{};

		[[nodiscard]] bool isIndexValid(const glm::ivec3& i) const;

		[[nodiscard]] int getFlatIndex(int i, int j, int k) const;
	};

	class BumperGrid {
	friend class WIPGrid;
	friend class SpatialGrid;
	public:
		const char IGNORE_OPP_SIDE_0 = 1 << 0;
		const char IGNORE_OPP_SIDE_1 = 1 << 1;
		const char IGNORE_OPP_SIDE_2 = 1 << 2;
		const char IGNORE_OPP_PLANE = 1 << 3;
		const char DONT_FOLLOW_SIDE_0 = 1 << 4;
		const char DONT_FOLLOW_SIDE_1 = 1 << 5;
		const char DONT_FOLLOW_SIDE_2 = 1 << 6;

		std::vector<Bumper> bumper;
		std::vector<Sphere> sphere;
		std::vector<GridSample> gridSamples;

		SpatialGrid grid;

		void constructFrom (const SphereMesh& sm);
		void constructGrid (int numberOfCells);

		bool pushOutsideBruteForce(glm::vec3& p, glm::vec3& n, int& bumperIndex)  const;
		bool pushOutsideGrid(glm::vec3& p, glm::vec3& n) const;
		bool projectOn(glm::vec3& p, glm::vec3& n, int& bumperIndex) const;

		void createGridSamplesFrom(const std::vector<glm::vec3>& pos, float skinThickness);

		void translate (const glm::vec3& t);
		void scale (float s);
		void rotateY (int angle);

		[[nodiscard]] bool serialize(const std::string& filepath) const;

		void deserializeCapsuloids(std::ifstream &file, int numSpheres, int numCapsuloids);
		void deserializePrysmoids(std::ifstream &file, int numSpheres, int numCapsuloids, int numPrysmoids);
		void deserializeQuads(std::ifstream &file, int numSpheres, int numCapsuloids, int numPrysmoids, int numQuads);
		void deserializeComposite(std::ifstream &file, int numSpheres, int numCapsuloids, int numPrysmoids, int numQuads,
		                          int numCompositeBumpers);
		void deserializeBumpers(std::ifstream& file, int numSpheres, int numCapsuloids, int numPrysmoids, int numQuads,
		                        int numCompositeBumpers);

		void deserialize(const std::string& filepath);

		void clearAllFlagsForTesting(); // for testing purposes only

#ifdef DEBUG
		void printOperationHistogram();
		void drawGrid () const;
		void exportGridSamplesToPLY(const std::string& filename) const;
		void printGridMaxCells () const;
		void printElement (int i) const;
		void printFlagsStat() const;
		[[nodiscard]] int parentOf(int i) const;
		void compareGrids(const BumperGrid &bg1) const;
#endif

	private:
		static float PUSH_EPSILON;

		AABB bbox;

		std::vector<std::set<int>> gridWIP;

		float projectOnBruteForce(glm::vec3& p, glm::vec3& n, int& bumperIndex)  const;

		bool constructionPushOutsideCapsuloid(glm::vec3& p, glm::vec3& n, int& bumperIndex)  const;
		bool constructionPushOutsidePrysmoid(glm::vec3& p, glm::vec3& n, int& bumperIndex)  const;
		bool constructionPushOutsideQuad(glm::vec3& p, glm::vec3& n, int& bumperIndex)  const;
		bool pushOutsideSphere(glm::vec3& p, glm::vec3& n, const int& bumperIndex)  const;

		bool pushOutsideCapsuloid(glm::vec3& p, glm::vec3& n, int bumperIndex)  const;
		bool pushOutsidePrysmoid(glm::vec3& p, glm::vec3& n, int bumperIndex, char flags = 0)  const;
		bool pushOutsideComposite(glm::vec3& p, glm::vec3 n, int bumperIndex) const;

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

		[[nodiscard]] std::pair<int, float> signedDistanceFromBumper(int i, const glm::vec3& p) const;
		void sampleDistanceFromBumperWithPosition(const glm::vec<3, float> &p, int proposedIndex, glm::vec3 &closestPos) const;

		void populateGrid(WIPGrid& wipGrid) const;

		std::vector<std::pair<Bumper, int>> sortCompositeIndices(CompositeBumper &cb);

		void addDontFollowFlags(CompositeBumper &cb, const std::vector<std::pair<Bumper, int>> &nodes, int i) const;

		Bumper makeCompositeBumper(const std::set<int> &cell, const std::set<int> &clean, int cellIdx);
		int makeOrFindCompositeBumper(const std::set<int>& orig, const std::set<int> &clean, int cellIdx);
		void processCompositeBumperFlags(CompositeBumper& cb, const std::set<int>& orig, int cellIdx);

		[[nodiscard]] std::set<int> compressChildren (const std::set<int>& cell) const;
		[[nodiscard]] std::set<int> compressSiblings (const std::set<int>& cell) const;
		[[nodiscard]] char computePrysmoidFlagsFromSet(int i, const std::set<int> &set) const;
		[[nodiscard]] char computePrysmoidFlagsFromCell(int element, int cellIdx) const;

		[[nodiscard]] int gridIndexOf(const std::set<int>& orig, int cellIdx);
		void computeFinalGrid();

		void initializeBumperSpheres();
		void initializeBumperCapsuloids(const SphereMesh &sm, std::vector<std::vector<int>> &capsuloidAdj);
		void initializeBumperPrysmoids(const SphereMesh &sm, std::vector<std::vector<int>> &capsuloidAdj);
		void initializeBumperQuads(const SphereMesh &sm, std::vector<std::vector<int>> &capsuloidAdj);
		void initializeBumperNodes(const SphereMesh &sm);

		[[nodiscard]] std::set<int> descendents(int i) const;

		void sortByType();
	};
}
