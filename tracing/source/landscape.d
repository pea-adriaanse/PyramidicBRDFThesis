module landscape;

import distribution;
import ply;
import pyramid_shape;
import ray;
import tracing;
import vdmath;
import vdmath.misc : degreesToRadians;

import std.algorithm.comparison : clamp;
import std.algorithm.sorting : sort;
import std.exception : enforce;
import std.math.algebraic : sqrt;
import std.math.exponential : log;
import std.math.rounding : floor;
import std.math.trigonometry : tan;
import std.parallelism : parallel;
import std.path : isValidFilename;
import std.random : uniform;
import std.typecons : Nullable;

struct Landscape {
	float width;
	float density;
	float approxHeight;
	Vec!3 minBound, maxBound;
	PyramidShape shape;

	Vec!3[] peaks;

	bool useBins = false;
	uint binCount;
	Vec!3[][] bins;
	Vec!2 binWidth;

	this(float width, float density = 0.6, float slope = degreesToRadians(54.7), Distribution heightDistribution = null) {
		setShape(width, density, slope);
		this.peaks = createPeaks(heightDistribution);
		this.density = peaks.length / width ^^ 2; // Update
		this.approxHeight = calcHeight(this.shape.slope, this.density);
	}

	this(float width, Vec!3[] peaks, float slope = degreesToRadians(54.7)) {
		setShape(width, peaks.length / (width ^^ 2), slope);
		this.peaks = peaks;
		this.approxHeight = calcHeight(this.shape.slope, this.density);
	}

	private void setShape(float width, float density, float slope) {
		this.width = width;
		this.density = density;
		this.shape = PyramidShape(slope);
	}

	static float calcHeight(float slope, float density, float heightError = 1e-5) {
		return tan(slope) * sqrt(-log(heightError) / (4.0 * density));
	}

	@disable this(Vec!3[] peaks, float slope) {
		this.peaks = peaks;
		this.shape = PyramidShape(slope);
	}

	@disable this(Vec!3[] peaks, const PyramidShape* shape) {
		this.peaks = peaks;
		this.shape = *shape;
	}

	Nullable!(Vec!(2, uint)) getBin(Vec!3 pos) const {
		Nullable!(Vec!(2, uint)) ret;
		static foreach (i; 0 .. 3)
			if (pos[i] < minBound[i] || pos[i] > maxBound[i])
				return ret; // null
		Vec!(3, double) local = (cast(Vec!(3, double)) pos) - minBound; // float rounding error causes issues
		uint binX = clamp(cast(uint) floor(local.x / binWidth.x), 0, binCount - 1);
		uint binY = clamp(cast(uint) floor(local.y / binWidth.y), 0, binCount - 1);
		ret = Vec!(2, uint)(binX, binY);
		return ret;
	}

	uint binIndex(uint x, uint y) const {
		return y * binCount + x;
	}

	const(Vec!3[]) getBin(uint x, uint y) const {
		assert(x < binCount && y < binCount);
		return bins[binIndex(x, y)];
	}

	/// Creates Default Bins
	/// Roof at 0, with simple width/2 bounds.
	/// Params:
	///   binCount = Number of bins along x & y axes.
	void createBins(uint binCount) {
		Vec!3 minBound = Vec!3(-width / 2, -width / 2, -approxHeight);
		Vec!3 maxBound = Vec!3(width / 2, width / 2, 0);
		createBins(minBound, maxBound, binCount);
	}

	void createBins(Vec!3 minBound, Vec!3 maxBound, uint binCount) {
		this.minBound = minBound;
		this.maxBound = maxBound;
		this.useBins = true;
		this.binCount = binCount;
		this.bins = new Vec!3[][binCount ^^ 2];

		float binXWidth = (maxBound.x - minBound.x) / binCount;
		float binYWidth = (maxBound.y - minBound.y) / binCount;
		this.binWidth = Vec!2(binXWidth, binYWidth);
		float cot_slope = shape.cot_slope;

		foreach (Vec!3 peak; peaks) {
			Vec!3 local = peak - minBound;
			float baseWidth = local.z * cot_slope;
			uint minXBin = cast(uint) clamp(floor(cast(int)(local.x - baseWidth) / binWidth.x), 0, binCount - 1);
			uint maxXBin = cast(uint) clamp(floor(cast(int)(local.x + baseWidth) / binWidth.x), 0, binCount - 1);
			uint minYBin = cast(uint) clamp(floor(cast(int)(local.y - baseWidth) / binWidth.y), 0, binCount - 1);
			uint maxYBin = cast(uint) clamp(floor(cast(int)(local.y + baseWidth) / binWidth.y), 0, binCount - 1);
			for (uint x = minXBin; x <= maxXBin; x++)
				for (uint y = minYBin; y <= maxYBin; y++)
					this.bins[binIndex(x, y)] ~= peak;
		}

		foreach (i, bin; bins)
			assert(bin.length > 0);
	}

	static Vec!3[] createGrid(float width, float density, float height, bool center = true) {
		uint root = cast(uint) sqrt(width * density);
		return createGrid(width, root, height, center);
	}

	static Vec!3[] createGrid(float width, uint cells, float height, bool center = true) {
		uint count = cells ^^ 2;
		float step = width / cells;
		Vec!3[] peaks = new Vec!3[count];
		foreach (x; 0 .. cells) {
			foreach (y; 0 .. cells) {
				uint index = x * cells + y;
				peaks[index] = Vec!3(x * step - width / 2.0, y * step - width / 2.0, height);
				if (center) {
					peaks[index].x += step / 2;
					peaks[index].y += step / 2;
				}
			}
		}
		return peaks;
	}

	static Vec!3[] createPoints(float width, uint count, float height) {
		Vec!3[] points = new Vec!3[count];
		foreach (i; 0 .. count) {
			float x = uniform!"()"(-width / 2, width / 2);
			float y = uniform!"()"(-width / 2, width / 2);
			points[i] = Vec!3(x, y, height);
		}
		return points;
	}

	/// Generate Pyramid Peaks
	/// Params:
	///   width = square width of sample (micrometers)
	///   density = density of pyramids (#/micrometer)
	///   testBurried = ensure peaks are not burried (ignored for constant height)
	///   heightDist = Distribution of pyramid peak heights, null = constant 0.
	Vec!3[] createPeaks(Distribution heightDist = null, bool testBurried = true) {
		uint count = cast(uint)(width * width * density);

		if (heightDist is null)
			return createPoints(width, count, 0);

		Vec!3[] peaks = new Vec!3[count];

		// Decide on heights first.
		float[] heights = new float[count];
		foreach (i; 0 .. count)
			heights[i] = heightDist.sample();

		// Sort to place largest first, better at preventing burrying.
		if (testBurried)
			sort!"a>b"(heights);

		uint i = 0;
		while (i < count) { // TODO: Prevent infinite loops
			Vec!3 newPeak;
			newPeak[0] = uniform!"()"(-width / 2, width / 2);
			newPeak[1] = uniform!"()"(-width / 2, width / 2);
			newPeak[2] = heights[i];
			if (testBurried) { // Cull if burried
				Hit hit = tracePeaks(shape, peaks, Ray(newPeak, Vec!3(0, 0, 1)));
				if (hit.hit)
					continue;
				peaks[i] = newPeak;
			} else
				peaks[i] = newPeak;
			i += 1;
		}
		return peaks;
	}

	void save(bool splitTriangles = true, bool globalBottom = false)(string fileName) {
		if (globalBottom)
			approxHeight = 1.0;
		enforce(isValidFilename(fileName), "File name invalid: " ~ fileName);
		Vec!3[5] shapeVertices = Vec!3(0, 0, 0) ~ shape.legs.dup;
		shapeVertices[] = shapeVertices[] * (-approxHeight / shape.legs[0].z);
		uint[3][4] shapeFaces;
		if (splitTriangles)
			shapeFaces = [[0, 1, 2], [3, 4, 5], [6, 7, 8], [9, 10, 11]];
		else
			shapeFaces = [[0, 1, 2], [0, 2, 3], [0, 3, 4], [0, 4, 1]];

		enum vertCount = splitTriangles ? 12 : 5;
		__gshared Vec!3[] vertices;
		vertices = new Vec!3[this.peaks.length * vertCount];
		__gshared uint[3][] faces;
		faces = new uint[3][this.peaks.length * 4];
		static if (splitTriangles) {
			__gshared Vec!3[] normals;
			normals = new Vec!3[vertices.length];
		}

		foreach (index, Vec!3 peak; parallel(this.peaks)) {
			Vec!3[5] pyramidVertices = shapeVertices;
			// pyramidVertices[] = pyramidVertices[] + peak; // Compiler bug
			foreach (i, ref v; pyramidVertices) {
				static if (globalBottom)
					v = (v * peak.z) + peak;
				else
					v = v + peak;
			}

			static if (splitTriangles) {
				vertices[index * vertCount .. (index + 1) * vertCount] = [
					pyramidVertices[0], pyramidVertices[1], pyramidVertices[2],
					pyramidVertices[0],
					pyramidVertices[2], pyramidVertices[3], pyramidVertices[0],
					pyramidVertices[3],
					pyramidVertices[4], pyramidVertices[0], pyramidVertices[4],
					pyramidVertices[1],
				];
				normals[index * vertCount .. (index + 1) * vertCount] = [
					shape.normalN, shape.normalN, shape.normalN, shape.normalW,
					shape.normalW,
					shape.normalW, shape.normalS, shape.normalS, shape.normalS,
					shape.normalE,
					shape.normalE, shape.normalE
				];
			} else
				vertices[index * vertCount .. (index + 1) * vertCount] = pyramidVertices;

			foreach (fi, f; shapeFaces) {
				uint[3] newFace;
				newFace[] = f[] + cast(uint)(index * vertCount);
				faces[index * shapeFaces.length + fi] = newFace;
			}
		}

		static if (splitTriangles)
			Ply.saveToFile(cast(immutable) vertices, cast(immutable) faces, cast(immutable) normals, fileName);
		else
			Ply.saveToFile(cast(immutable) vertices, cast(immutable) faces, fileName);
	}
}
