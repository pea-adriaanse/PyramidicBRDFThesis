module landscape;

import ply;
import vdmath;
import pyramid_shape;

import std.algorithm.comparison : clamp;
import std.math.rounding : floor;
import std.parallelism : parallel;
import std.path : isValidFilename;
import std.typecons : Nullable;
import std.exception : enforce;
import distribution;

struct Landscape {
	Vec!3[] peaks;
	float width;
	Vec!3 minBound, maxBound;
	float roof;
	bool useBins = false;
	PyramidShape shape;

	uint binCount;
	Vec!3[][] bins;
	Vec!2 binWidth;

	this(float width, float density = 0.6, float slope = 54.7) {
		this.width = width;
		this.peaks = createPeaks(width, density);
		this.shape = PyramidShape(slope);
	}

	this(float width, Vec!3[] peaks, float slope = 54.7) {
		this.width = width;
		this.peaks = peaks.dup;
		this.shape = PyramidShape(slope);
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
		Vec!3 local = pos - minBound;
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
			uint minXBin = clamp(cast(uint) floor((local.x - baseWidth) / binWidth.x), 0, binCount - 1);
			uint maxXBin = clamp(cast(uint) floor((local.x + baseWidth) / binWidth.x), 0, binCount - 1);
			uint minYBin = clamp(cast(uint) floor((local.y - baseWidth) / binWidth.y), 0, binCount - 1);
			uint maxYBin = clamp(cast(uint) floor((local.y + baseWidth) / binWidth.y), 0, binCount - 1);
			for (uint x = minXBin; x <= maxXBin; x++)
				for (uint y = minYBin; y <= maxYBin; y++)
					this.bins[binIndex(x, y)] ~= peak;
		}
	}

	static Vec!3[] createGrid(float width, uint cells) {
		uint count = cells ^^ 2;
		float step = width / cells;
		Vec!3[] peaks = new Vec!3[count];
		foreach (x; 0 .. cells) {
			foreach (y; 0 .. cells) {
				uint index = x * cells + y;
				peaks[index] = Vec!3(x * step - width / 2.0, y * step - width / 2.0, 0);
			}
		}
		return peaks;
	}

	/// Generate Pyramid Peaks
	/// Params:
	///   width = square width of sample (micrometers)
	///   density = density of pyramids (#/micrometer)
	///   testBurried = ensure peaks are not burried (ignored for constant height)
	///   heightDist = Distribution of pyramid peak heights, null = constant 0.
	static Vec!3[] createPeaks(float width, float density, bool testBurried = true, Distribution heightDist) {
		uint count = cast(uint)(width * width * density);

		if (heightDist is null) {
			Vec!3[] peaks = new Vec!3[count];
			foreach (i; 0 .. peaks) {
				float x = uniform!"()"(-width / 2, width / 2);
				float y = uniform!"()"(-width / 2, width / 2);
				peaks[i] = Vec!3(x, y, 0);
			}
			return peaks;
		}

		Vec!3[] peaks = new Vec!3[count];

		// Decide on heights first.
		float[] heights = new float[count];
		foreach (i; 0 .. count)
			heights[i] = heightDistribution.sample();

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

	void save(bool splitTriangles = true)(string fileName, float height) {
		enforce(isValidFilename(fileName), "File name invalid: " ~ fileName);
		// if (height == 0) {
		// 	float error = 1e-5;
		// 	height = ceil(tan(shape.slope) * sqrt(-log(error) / (4.0 * _density))); // ceil optional
		// }
		Vec!3[5] shapeVertices = Vec!3(0, 0, 0) ~ shape.legs.dup;
		shapeVertices[] = shapeVertices[] * (-height / shape.legs[0].z);
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
			// pyramidVertices[] = pyramidVertices[] + peak; // Compiler bug'
			foreach (ref v; pyramidVertices)
				v = v + peak;

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
