import math;
import misc;
import ply;
import sphere_samplers;
import std.algorithm.comparison : clamp;
import std.algorithm.searching : canFind;
import std.algorithm.sorting : sort;
import std.algorithm.sorting;
import std.conv : to;
import std.conv : to;
import std.exception : enforce;
import std.math;
import std.parallelism;
import std.random;
import std.stdio;
import std.path;
import std.string : strip;

struct Face {
	Vec!3[2] legs;
	Vec!3 normal;
}

struct PyramidShape {
	float slope;
	float cot_slope;

	union {
		struct {
			Vec!3 NE;
			Vec!3 NW;
			Vec!3 SW;
			Vec!3 SE;
		}

		Vec!3[4] legs;
	}

	union {
		struct {
			Vec!3 normalE;
			Vec!3 normalN;
			Vec!3 normalW;
			Vec!3 normalS;
		}

		Vec!3[4] normals;
	}

	union {
		struct {
			Face E;
			Face N;
			Face W;
			Face S;
		}

		Face[4] faces;
	}

	this(double slope) {
		this.slope = cast(float) slope;
		double cot = 1.0 / tan(slope);
		this.cot_slope = cast(float) cot;

		Vec!(3, double) leg = (Vec!(3, double)(cot, cot, -1)).normalize!double();
		this.NE = cast(Vec!3) leg;
		this.NW = [-leg[0], leg[0], leg[2]];
		this.SW = [-leg[0], -leg[0], leg[2]];
		this.SE = [leg[0], -leg[0], leg[2]];

		// this.legAngleCos = cos(2 * atan(cos(slope))); // == 2/(cos²(slope)+1)-1

		double normalZ = cos(slope);
		double normalX = sin(slope);
		this.normalE = Vec!3(normalX, 0, normalZ);
		this.normalN = Vec!3(0, normalX, normalZ);
		this.normalW = Vec!3(-normalX, 0, normalZ);
		this.normalS = Vec!3(0, -normalX, normalZ);

		this.E = Face([SE, NE], normalE);
		this.N = Face([NE, NW], normalN);
		this.W = Face([NW, SW], normalW);
		this.S = Face([SW, SE], normalS);
	}
}

alias Landscape = Vec!3[];

string toStr(float f) {
	return to!string(f);
}

interface Distribution {
	float sample();
}

class ConstantDistribution : Distribution {
	float height;
	this(float height) {
		this.height = height;
	}

	float sample() {
		return this.height;
	}
}

/// Models distribution using histogram
class HistogramDistribution : Distribution {
	float[] cdf = [];
	float rangeSize; // relates #buckets to range.
	float bucketSize;

	this(string file, float rangeSize) {
		this.rangeSize = rangeSize;
		foreach (char[] l; File(file).byLine()) {
			this.cdf ~= l.to!float;
		}
		this.bucketSize = rangeSize / cdf.length;
	}

	float sample() {
		float prob = uniform01();
		foreach_reverse (bucket, cumulativeProb; cdf) {
			if (cumulativeProb <= prob) {
				return (bucket + 1) * bucketSize;
			}
		}
		assert(0);
	}
}

PyramidShape shape = PyramidShape(degreesToRadians(_slope));
enum _slope = 54.7;
enum _width = 20.0; //? in micron
enum float width = 40; // used
enum _density = 0.6; //? per micron²

uint simpleSphereLatCount = 16;
uint simpleSphereLongCount = 10;
uint fibonacciSphereCount = 400;

uint sphereSampleCount; // determined given sampler used ^

uint sampleCount = 4000;
float maxVarHeight = 7.0;
float maxConstHeight = 1.0;
uint heightBins = 10;
uint landWidthSampleCount = cast(uint)(5 * width);

void main(string[] args) {
	if (args.length <= 1)
		return writeln("Choose:\n\t- experiment\n\t- generate\n\t- test");
	switch (args[1]) {
		case "experiment":
			return experiment();
		case "generate":
			write("landscape size: ");
			float size = readln().strip().to!float;
			enforce(size > 0, "Size cannot be negative.");
			write("output file: ");
			string fileName = readln().strip();
			enforce(isValidFilename(fileName), "File name invalid.");
			return generate(size, fileName);
		case "test":
			return test_shadowing();
		default:
			return writeln("Choose:\n\t- experiment\n\t- generate\n\t- test");
	}
}

void generate(bool splitTriangles = true)(float size, string fileName) {
	import std.datetime.stopwatch;

	StopWatch watch = StopWatch(AutoStart.yes);

	float error = 0.001;
	float L = ceil(tan(shape.slope) * sqrt(-log(error) / (4.0 * _density))); // ceil optional

	Distribution heightDistribution = new ConstantDistribution(L);
	Landscape land = createLandscape(size, _density, heightDistribution, false);

	Vec!3[5] shapeVertices = Vec!3(0, 0, 0) ~ shape.legs.dup;
	shapeVertices[] = shapeVertices[] * (-L / shape.legs[0].z);
	uint[3][4] shapeFaces;
	if (splitTriangles)
		shapeFaces = [[0, 1, 2], [3, 4, 5], [6, 7, 8], [9, 10, 11]];
	else
		shapeFaces = [[0, 1, 2], [0, 2, 3], [0, 3, 4], [0, 4, 1]];

	enum vertCount = splitTriangles ? 12 : 5;
	__gshared Vec!3[] vertices;
	vertices = new Vec!3[land.length * vertCount];
	__gshared uint[3][] faces;
	faces = new uint[3][land.length * 4];
	static if (splitTriangles) {
		__gshared Vec!3[] normals;
		normals = new Vec!3[vertices.length];
	}

	foreach (index, Vec!3 peak; parallel(land)) {
		Vec!3[5] pyramidVertices = shapeVertices;
		// pyramidVertices[] = pyramidVertices[] + peak; // Compiler bug'
		foreach (ref v; pyramidVertices)
			v = v + peak;

		static if (splitTriangles) {
			vertices[index * vertCount .. (index + 1) * vertCount] = [
				pyramidVertices[0], pyramidVertices[1], pyramidVertices[2], pyramidVertices[0],
				pyramidVertices[2], pyramidVertices[3], pyramidVertices[0], pyramidVertices[3],
				pyramidVertices[4], pyramidVertices[0], pyramidVertices[4], pyramidVertices[1],
			];
			normals[index * vertCount .. (index + 1) * vertCount] = [
				shape.normalN, shape.normalN, shape.normalN, shape.normalW, shape.normalW,
				shape.normalW, shape.normalS, shape.normalS, shape.normalS, shape.normalE,
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
	writeln("Mesh created");

	static if (splitTriangles)
		Ply.saveToFile(cast(immutable) vertices, cast(immutable) faces, cast(immutable) normals, fileName);
	else
		Ply.saveToFile(cast(immutable) vertices, cast(immutable) faces, fileName);
	writeln("PLY File created");
	writeln("Milliseconds: ", watch.peek().total!"msecs");
}

void test_shadowing() {
	Vec!3 zeroAzimuth(Vec!3 r) {
		float x = sqrt(1.0f - r.z * r.z);
		return Vec!3(x, 0, r.z);
	}

	float shadow_lyanne(Vec!3 o, Vec!3 v) {
		if (o.dot(v) < 0.0)
			return 0.0;
		Vec!3 r = zeroAzimuth(o);

		float rg = r.z;
		if (rg < 0.0)
			return 0.0;
		float D = 1.0 / (4.0 * cos(degreesToRadians(54.7)));
		float rfe = (r.z * shape.normalE.z) + (r.x * shape.normalE.x);
		float rfo = r.z * shape.normalN.z;

		float numerator = rg;
		float denominator = D * (rfe + 2 * rfo);
		float shadowing = numerator / denominator;
		if (shadowing > 1)
			return 1.0;
		return shadowing;
	}

	File file = File("results/shadow_lyanne.csv", "w");
	file.writeln("theta,shadow");
	uint steps = 200; // >= 2
	foreach (i; 0 .. steps) {
		float theta = i * PI_2 / (steps - 1.0);
		Vec!3 o = Vec!3(sin(theta), 0, cos(theta));
		float shadow = shadow_lyanne(o, shape.normalE);
		file.writeln(theta, ",", shadow);
	}
}

void experiment() {
	string heightsFile = "heightCumulative.csv"; // Source: "Opto-electrical modelling and optimization study of a novel IBC c-Si solar cellOpto-electrical modelling and optimization study of a novel IBC c-Si solar cell"
	Distribution heightDistribution = new HistogramDistribution(heightsFile, maxVarHeight);
	Distribution heightDistributionConstant = new ConstantDistribution(maxConstHeight);
	Landscape land = createLandscape(width, _density, heightDistribution); // 20 micron, 0.6 per micron²
	Landscape landConstant = createLandscape(width, _density, heightDistributionConstant);

	// SphereSampler sphereSampler = new SphereFibonacci(sphereSampleCount);
	SphereSampler sphereSampler = new SphereSimple(simpleSphereLongCount, simpleSphereLatCount, true);
	sphereSampleCount = cast(uint) sphereSampler.data.length;
	sphereSampler.save("results/sphereSamples.csv");

	File settingsFile = File("results/settings.csv", "w");
	settingsFile.writeln("sphereSampler,", sphereSampler.name);
	settingsFile.writeln("sphereSampleCount,", sphereSampleCount.to!string);
	settingsFile.writeln("sampleCount,", sampleCount.to!string);
	settingsFile.writeln("landWidthSampleCount,", landWidthSampleCount.to!string);
	settingsFile.writeln("width,", width);
	settingsFile.writeln("maxConstHeight,", maxConstHeight);
	settingsFile.writeln("heightBins,", heightBins);
	settingsFile.writeln("simpleSphereLongCount,", simpleSphereLongCount);
	settingsFile.writeln("simpleSphereLatCount,", simpleSphereLatCount);

	measureLand(landConstant, "results/land.csv");
	// measure(land, sphereSamples, maxVarHeight, "hits.csv");
	measure(landConstant, sphereSampler, maxConstHeight, heightBins, "results/hitsConstant.csv",
		"results/hitsHeightConstant.csv", "results/heightDistributionConstant.csv",
		"results/hitsConstantByHeight.csv");

	// measureHighPoint(landConstant, sphereSamples, "results/highPointSampling.csv");
}

void measureHighPoint(const Landscape land, const Vec!3[] sphereSamples, string fileName) {
	Vec!3 sample;
	uint samplePeakID;
	bool goodSample = false;
	while (!goodSample) {
		Vec!3 org = Vec!3(uniform(-width / 2, width / 2), uniform(-width / 2, width / 2), maxConstHeight + 1);
		Ray ray = Ray(org, Vec!3(0, 0, -1));
		Hit hit = trace(land, ray);
		sample = hit.pos;
		samplePeakID = hit.peakID;
		goodSample = sample.z > 0.85 * maxConstHeight && hit.face == 0;
	}
	measureTracing(land, sphereSamples, [sample], [samplePeakID], fileName);
}

void measureTracing(const Landscape land, const Vec!3[] sphereSamples, const Vec!3[] samples,
	const uint[] samplePeakIDs, string fileName) {
	File file = File(fileName, "w");
	file.writeln("orgx,orgy,orgz,dirPointx,dirPointy,dirPointz,hit");
	foreach (i, Vec!3 sample; samples) {
		foreach (Vec!3 sphereSample; sphereSamples) {
			Ray ray = Ray(sample, sphereSample, samplePeakIDs[i]);
			Hit hit = trace(land, ray);

			bool hitPeak = hit.hit == false; // Not obscured
			Vec!3 dirPoint = sample + sphereSample;
			file.writeln(sample.x, ",", sample.y, ",", sample.z, ",", dirPoint.x, ",",
				dirPoint.y, ",", dirPoint.z, ",", hitPeak);
		}
	}
}

void measureLand(const Landscape land, string fileName) {
	File file = File(fileName, "w");
	file.writeln("x,y,z");
	foreach (x; 0 .. landWidthSampleCount) {
		foreach (y; 0 .. landWidthSampleCount) {
			Ray ray = Ray(Vec!3(x * width / landWidthSampleCount - width / 2,
					y * width / landWidthSampleCount - width / 2, maxConstHeight + 1), Vec!3(0, 0, -1));
			Hit hit = trace(land, ray);
			assert(hit.hit);
			file.writeln(hit.pos.x, ",", hit.pos.y, ",", hit.pos.z);
		}
	}
}

void measure(const Landscape land, SphereSampler sphereSampler, float maxHeight, uint heightBins,
	string outFile, string outHeightFile, string outHeightDistributionFile, string outHitsByHeightFile) {
	uint[] hits = new uint[sphereSampleCount];
	uint[2][] hitsByHeight = new uint[2][sphereSampleCount * heightBins];
	uint[2][] heightHits = new uint[2][heightBins];
	Vec!3[] samplePoss;
	uint[] samplePeakIDs;

	float minSampleHeight = float.infinity;
	float maxSampleHeight = -float.infinity;

	// Create Samples
	samplePoss.reserve(sampleCount);
	samplePeakIDs.reserve(sampleCount);
	foreach (i; 0 .. sampleCount) {
		bool east = false;
		Vec!3 pos;
		Hit hit;
		while (!east) { // Filter for east faces
			Vec!3 org = Vec!3(uniform(-width / 2, width / 2), uniform(-width / 2, width / 2), maxHeight + 1);
			Ray ray = Ray(org, Vec!3(0, 0, -1));
			hit = trace(land, ray);
			assert(hit.hit);
			org.assertAlmostEq(Vec!3(hit.pos.x, hit.pos.y, maxHeight + 1));
			east = hit.face == 0;
			pos = hit.pos;
		}
		samplePoss ~= pos;
		samplePeakIDs ~= hit.peakID;

		if (pos.z < minSampleHeight)
			minSampleHeight = pos.z;
		if (pos.z > maxSampleHeight)
			maxSampleHeight = pos.z;
	}

	// Bin sample heights
	uint[] bins = new uint[heightBins];
	foreach (p; samplePoss) {
		uint heightBin = cast(uint) floor(heightBins * (p.z - minSampleHeight) / (maxSampleHeight - minSampleHeight));
		heightBin = clamp(heightBin, 0, heightBins - 1); // Floating error precaution
		bins[heightBin] += 1;
	}

	// Sample for visibility
	foreach (sampleID, samplePos; samplePoss) {
		// Loading bar
		write("[", sampleID.to!string, "/", sampleCount.to!string, "]\t\t\r");
		stdout.flush();

		foreach (sphereID, sphereSample; sphereSampler.data) {
			Ray ray = Ray(samplePos, sphereSample, samplePeakIDs[sampleID]); // Assuming East Face
			Hit hit = trace(land, ray);

			uint heightBin = cast(uint) floor(
				heightBins * (samplePos.z - minSampleHeight) / (maxSampleHeight - minSampleHeight));
			heightBin = clamp(heightBin, 0, heightBins - 1); // Floating error precaution
			heightHits[heightBin][0] += 1;
			hitsByHeight[heightBin * sphereSampleCount + sphereID][0] += 1;
			if (!hit.hit) { // Not obscured
				heightHits[heightBin][1] += 1;
				hits[sphereID] += 1;
				hitsByHeight[heightBin * sphereSampleCount + sphereID][1] += 1;
			}
		}
	}

	// Save results
	File hitsFile = File(outFile, "w");
	hitsFile.writeln("org_id,hits");
	foreach (i, h; hits)
		hitsFile.writeln(i.to!string, ",", h.to!string);

	File hitsByHeightFile = File(outHitsByHeightFile, "w");
	hitsByHeightFile.writeln("lat,hits");
	// Processing hits by height
	foreach (heightBin; 0 .. heightBins) {
		uint[2][float] binHits;
		foreach (sphereID; 0 .. sphereSampleCount) {
			uint index = heightBin * sphereSampleCount + sphereID;
			float lat = sphereSampler.angles[sphereID][1];
			if (lat !in binHits)
				binHits[lat] = [0, 0];
			binHits[lat][0] += hitsByHeight[index][0];
			binHits[lat][1] += hitsByHeight[index][1];
		}
		foreach (lat; binHits.keys.sort()) {
			if (binHits[lat][0] == 0)
				hitsByHeightFile.writeln(lat, ",", 0);
			else
				hitsByHeightFile.writeln(lat, ",", binHits[lat][1] / cast(float) binHits[lat][0]);
		}
	}

	File heightHitsFile = File(outHeightFile, "w");
	heightHitsFile.writeln("occurrance,hit");
	foreach (h; heightHits)
		heightHitsFile.writeln(h[0].to!string, ",", h[1].to!string);

	File sampleHeightsFile = File(outHeightDistributionFile, "w");
	sampleHeightsFile.writeln(
		"count,minHeight=" ~ minSampleHeight.to!string ~ ",maxHeight=" ~ maxSampleHeight.to!string);
	foreach (i, b; bins) {
		// if (i == 0) {
		// 	sampleHeightsFile.writeln(p.z.to!string, ",", minSampleHeight.to!string, ",", maxSampleHeight.to!string);
		// }
		sampleHeightsFile.writeln(b.to!string);
	}
}

/// Generate Pyramid Landscape
/// Params:
///   width = square width of sample (micrometers)
///   density = density of pyramids (#/micrometer)
///   heightDistribution = Distribution of pyramid peak heights
///   testBurried = test if peaks are burried under surface & cull. (Not necessary at constant peak heights)
/// Returns: Pyramid Landscape
Landscape createLandscape(float width, float density, Distribution heightDistribution, bool testBurried = true) {
	uint count = cast(uint)(width * width * density);
	Vec!3[] peaks;
	peaks.reserve(count);

	// Decide on heights first.
	float[] heights = new float[count];
	foreach (i; 0 .. count)
		heights[i] = heightDistribution.sample();

	// Sort to place largest first.
	sort!"a>b"(heights);

	while (peaks.length < count) {
		uint i = cast(uint) peaks.length;
		Vec!3 newPeak;
		newPeak[0] = uniform!"()"(-width / 2, width / 2);
		newPeak[1] = uniform!"()"(-width / 2, width / 2);
		newPeak[2] = heights[i];

		if (testBurried) { // Cull if burried
			Hit hit = trace(peaks, Ray(newPeak, Vec!3(0, 0, 1)));
			if (!hit.hit)
				peaks ~= newPeak;
		} else {
			peaks ~= newPeak;
		}
	}
	return peaks;
}

struct Ray {
	Vec!3 org;
	Vec!3 dir;
	uint exculdePeak = uint.max;
}

struct Hit {
	bool hit = false;
	Vec!3 pos;
	float t = float.infinity;
	uint peakID;
	ubyte face;
}

Hit trace(const Landscape land, Ray ray) {
	shared Hit minHit;
	foreach (id, peak; parallel(land)) {
		if (ray.exculdePeak == id)
			continue;
		Hit hit = trace(peak, ray);
		if (hit.hit) {
			hit.peakID = cast(uint) id;
			synchronized if (hit.t < minHit.t)
				minHit = cast(shared Hit) hit;
		}
	}
	return minHit;
}

Hit trace(Vec!3 peak, Ray ray) {
	Hit minHit;
	foreach (ubyte i, Face f; shape.faces) {
		Hit hit = trace(peak, f, ray);
		if (hit.hit && hit.t < minHit.t) {
			minHit = hit;
			minHit.face = i;
		}
	}
	return minHit;
}

unittest {
	shape = PyramidShape(degreesToRadians(45.0));
	Hit hit = trace(Vec!3(0, 0.25, 0), Ray(Vec!3(0, -1, 1), Vec!3(0, 0, -1)));
	assert(hit.hit);
	writeln(hit.pos);
	assert(hit.pos == Vec!3(0, -1, -1.25));
}

unittest {
	shape = PyramidShape(degreesToRadians(54.7));
	Hit hit = trace(Vec!3(0, 0, 1), Ray(Vec!3(0), Vec!3(0, 0, 1)));
	assert(hit.hit);
}

Hit trace(Vec!3 peak, Face face, Ray ray) {
	Hit hit = tracePlane(peak, face.normal, ray);
	if (!hit.hit)
		return hit;
	Vec!3 P = hit.pos - peak;

	float delta = 1e-5;
	float sliceWidth = -shape.cot_slope * P.z + delta; // equals abs x or y.

	hit.hit = abs(P.x) <= sliceWidth && abs(P.y) <= sliceWidth;
	// float Pu = P.dot(face.legs[0]);
	// float Pv = (P - face.legs[0] * Pu).dot(face.legs[1]);
	// bool inside = Pu >= 0 && Pv >= 0;
	return hit;
}

Hit tracePlane(Vec!3 point, Vec!3 normal, Ray ray) {
	float distToPlane = (point - ray.org).dot(normal);
	float slope = ray.dir.dot(normal);

	// version (backfaceCulling) { //WARNING version breaks createLandscape culling
	// 	if (slope > 0)
	// 		return Hit(false);
	// }
	if (slope == 0 || (signbit(slope) != signbit(distToPlane)))
		return Hit(false);

	float t = distToPlane / slope;
	Vec!3 pos = ray.org + ray.dir * t;
	return Hit(true, pos, t);
}
