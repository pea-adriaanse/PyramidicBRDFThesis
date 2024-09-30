import vdmath;
import vdmath.misc;
import ply;
import csv;

import sphere_samplers;
import std.algorithm.comparison : clamp;
import std.algorithm.searching : canFind;
import std.algorithm.sorting : sort;
import std.algorithm.sorting;
import std.conv : to;
import std.conv : to;
import std.exception : enforce, ErrnoException;
import std.file : exists;
import std.math;
import std.parallelism;
import std.path;
import std.random;
import std.stdio;
import std.string : strip;
import std.typecons : Nullable;
import std.datetime.stopwatch;

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

struct Landscape {
	Vec!3[] peaks;
	Vec!3 minBound, maxBound;
	float roof;
	bool useBins = false;

	uint binCount;
	Vec!3[][] bins;
	Vec!2 binWidth;

	this(Vec!3[] peaks) {
		this.peaks = peaks;
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
}

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
	StopWatch stopwatch = StopWatch(AutoStart.yes);
	scope (exit)
		writeln("Finished in ", stopwatch.peek().total!"msecs"(), "milliseconds");
	if (args.length <= 1)
		return writeln("Choose:\n\t- reflectDist\n\t- experiment\n\t- generate\n\t- test");
	switch (args[1]) {
	case "reflectDist":
		return measureReflectPathDist(args[2 .. $]);
	case "experiment":
		return experiment();
	case "generate":
		return generate();
	case "test":
		return test_shadowing();
	case "measure":
		return measure_shadowing();
	default:
		return writeln("Choose:\n\t- reflectDist\n\t- experiment\n\t- generate\n\t- test");
	}
}

uint pow4(uint exponent) {
	return 1u << (2u * exponent);
}

uint pow4sum(uint exponent) {
	return (4 * (pow4(exponent) - 1)) / 3; // Modified geometric series
}

unittest {
	assert(pow4sum(0u) == 0);
	assert(pow4sum(1u) == 4);
	assert(pow4sum(2u) == 20);
	assert(pow4sum(3u) == 84);
	assert(pow4sum(4u) == 340);
}

void measureReflectPathDist(string[] args) {
	rndGen().seed(0);

	float width = 1000;
	float error = 0.001;
	const uint binCount = 25;
	float L = ceil(tan(shape.slope) * sqrt(-log(error) / (4.0 * _density))); // ceil optional
	Distribution heightDistribution = new ConstantDistribution(L);
	Landscape land = createLandscape(width, _density, heightDistribution, false, true);
	land.createBins(Vec!3(-width - L, -width - L, -L - 1), Vec!3(width + L, width + L, L + 1), binCount);

	Vec!3 wo;
	uint sampleCount, reflectCount;
	string identifier;
	const argsError = "Expect 0 or 6 arguments: (wo.x wo.y wo.z sampleCount reflectCount identifier) but got: " ~ args
		.to!string;
	if (args.length == 6) {
		wo.x = args[0].to!float;
		wo.y = args[1].to!float;
		wo.z = args[2].to!float;
		sampleCount = args[3].to!uint;
		reflectCount = args[4].to!uint;
		identifier = args[5];
	} else {
		enforce(args.length == 0, argsError);
		write("wo.x: ");
		wo.x = readln().strip().to!float;
		write("wo.y: ");
		wo.y = readln().strip().to!float;
		write("wo.z: ");
		wo.z = readln().strip().to!float;
		write("sample count: ");
		sampleCount = readln().strip().to!uint;
		write("reflect count: ");
		reflectCount = readln().strip().to!uint;
		write("output file identifier: ");
		identifier = readln().strip();
	}
	writeln("\nwo: ", wo.toString());
	Vec!3 dir = -wo;

	uint optionCount = pow4sum(reflectCount);
	uint[] pathCounts = new uint[optionCount];
	uint captureCount = 0;

	Vec!2[] offsets = new Vec!2[sampleCount];
	foreach (i; 0 .. sampleCount) {
		float x = uniform(-width / 3.0, width / 3.0);
		float y = uniform(-width / 3.0, width / 3.0);
		offsets[i] = Vec!2(x, y);
	}

	CSV csv = CSV(',', false, "index", "indexStr", "count", "prob");

	foreach (i; 0 .. sampleCount) {
		debug writef("\rRaytracing: %u/%u", i, sampleCount);
		float minHeight = L / abs(dir.z);

		Vec!3 target = Vec!3(offsets[i].x, offsets[i].y, 0);
		Vec!3 cam = (target - dir * 2 * minHeight);

		Ray ray = Ray(cam, dir);
		ReflectData reflectData = reflectRecurse(land, ray, reflectCount);

		if (!reflectData.exits) {
			captureCount += 1;
			continue;
		} else {
			pathCounts[reflectData.reflectID] += 1;
		}
	}
	debug writefln("\rRaytracing: %u/%u", sampleCount, sampleCount);

	foreach (uint i; 0 .. optionCount) {
		csv.addEntry(i, reflectIDToString(i), pathCounts[i], pathCounts[i].to!float / sampleCount);
	}
	csv.save("../scripts/temp/dist_" ~ identifier ~ "_S.csv");
}

void measure_shadowing() {
	Landscape land = generate_land();

	// TODO: replace region with parallel input direction
	writeln("Region:");
	write("x_min: ");
	float xMin = readln().strip().to!float;
	write("x_max: ");
	float xMax = readln().strip().to!float;

	write("y_min: ");
	float yMin = readln().strip().to!float;
	write("y_max: ");
	float yMax = readln().strip().to!float;

	write("cam_x: ");
	float camX = readln().strip().to!float;
	write("cam_y: ");
	float camY = readln().strip().to!float;
	write("cam_z: ");
	float camZ = readln().strip().to!float;
	Vec!3 cam = Vec!3(camX, camY, camZ);

	write("sample count: ");
	uint sampleCount = readln().strip().to!uint;

	write("reflect count: ");
	uint reflectCount = readln().strip().to!uint;

	Vec!2[] samples = new Vec!2[sampleCount];
	foreach (i; 0 .. sampleCount) {
		float x = uniform(xMin, xMax);
		float y = uniform(yMin, yMax);
		samples[i] = Vec!2(x, y);
	}

	// CSV csv = CSV(',', true, "exits", "reflectCount", "reflectID", "outDirX", "outDirY", "outDirZ");
	CSV csv = CSV(',', true, "exitID", "count");
	uint[uint] data;

	foreach (i; 0 .. sampleCount) {
		Vec!3 target = Vec!3(samples[i].x, samples[i].y, 0);
		Vec!3 dir = (target - cam).normalize();
		Ray ray = Ray(cam, dir);

		ReflectData reflectData = reflectRecurse(land, ray, reflectCount);
		data[reflectData.reflectID] += reflectData.exits;
		// csv.addEntry(reflectData.exits, reflectData.reflectCount, reflectData.reflectID,
		// 	reflectData.outDir.x, reflectData.outDir.y, reflectData.outDir.z);
	}

	foreach (d, v; data) {
		csv.addEntry(d, v);
	}

	csv.save("results/measure_shadowing.csv");
}

Landscape generate_land() {
	write("landscape size (um): ");
	float size = readln().strip().to!float;
	enforce(size > 0, "Size cannot be negative.");
	write("landscape scale (unit): ");
	float scale = readln().strip().to!float;
	enforce(scale > 0, "Scale cannot be negative.");

	float scaling = scale / size; // units/micron
	float error = 0.001;
	float L = scaling * ceil(tan(shape.slope) * sqrt(-log(error) / (4.0 * _density))); // ceil optional
	writeln("Height: ", L);
	Distribution heightDistribution = new ConstantDistribution(L);
	Landscape land = createLandscape(scaling * size, _density / (scaling ^^ 2), heightDistribution, false);
	return land;
}

void generate(bool splitTriangles = true)() {
	write("landscape size (um): ");
	float size = readln().strip().to!float;
	enforce(size > 0, "Size cannot be negative.");
	write("landscape scale (unit): ");
	float scale = readln().strip().to!float;
	enforce(scale > 0, "Scale cannot be negative.");
	write("output file: ");
	string fileName = readln().strip();
	enforce(isValidFilename(fileName), "File name invalid.");

	import std.datetime.stopwatch;

	float scaling = scale / size; // units/micron

	StopWatch watch = StopWatch(AutoStart.yes);

	float error = 0.001;
	float L = scaling * ceil(tan(shape.slope) * sqrt(-log(error) / (4.0 * _density))); // ceil optional

	Distribution heightDistribution = new ConstantDistribution(L);
	Landscape land = createLandscape(scaling * size, _density / (scaling ^^ 2), heightDistribution, false, true);

	Vec!3[5] shapeVertices = Vec!3(0, 0, 0) ~ shape.legs.dup;
	shapeVertices[] = shapeVertices[] * (-L / shape.legs[0].z);
	uint[3][4] shapeFaces;
	if (splitTriangles)
		shapeFaces = [[0, 1, 2], [3, 4, 5], [6, 7, 8], [9, 10, 11]];
	else
		shapeFaces = [[0, 1, 2], [0, 2, 3], [0, 3, 4], [0, 4, 1]];

	enum vertCount = splitTriangles ? 12 : 5;
	__gshared Vec!3[] vertices;
	vertices = new Vec!3[land.peaks.length * vertCount];
	__gshared uint[3][] faces;
	faces = new uint[3][land.peaks.length * 4];
	static if (splitTriangles) {
		__gshared Vec!3[] normals;
		normals = new Vec!3[vertices.length];
	}

	foreach (index, Vec!3 peak; parallel(land.peaks)) {
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
			org.assertAlmostEquals(Vec!3(hit.pos.x, hit.pos.y, maxHeight + 1));
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
		uint heightBin = cast(uint) floor(
			heightBins * (
				p.z - minSampleHeight) / (maxSampleHeight - minSampleHeight));
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
				heightBins * (
					samplePos.z - minSampleHeight) / (maxSampleHeight - minSampleHeight));
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
		"count,minHeight=" ~ minSampleHeight.to!string ~ ",maxHeight=" ~ maxSampleHeight
			.to!string);
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
Landscape createLandscape(float width, float density, Distribution heightDistribution, bool testBurried = true, bool grid = false) {
	uint count = cast(uint)(width * width * density);
	if (grid)
		count = (cast(uint) sqrt(width * width * density)) ^^ 2;

	Vec!3[] peaks;
	peaks.reserve(count);

	// Decide on heights first.
	float[] heights = new float[count];
	foreach (i; 0 .. count)
		heights[i] = heightDistribution.sample();

	// Sort to place largest first.
	sort!"a>b"(heights);

	if (grid) {
		uint root = cast(uint) sqrt(cast(float) count);
		float step = width / root;
		foreach (x; 0 .. root) {
			foreach (y; 0 .. root) {
				uint index = x * root + y;
				peaks ~= Vec!3(x * step, y * step, heights[index]);
			}
		}
	} else {
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
	}
	return Landscape(peaks);
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

Hit trace(const Vec!3[] peaks, Ray ray) {
	shared Hit minHit;
	foreach (id, peak; parallel(peaks)) {
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

bool findNextBin(const Landscape land, ref Ray ray, ref Nullable!(Vec!(2, uint)) currentBin) {
	assert(land.useBins);
	if (currentBin.isNull) {
		if (ray.org.z > land.maxBound.z) { // Outside z bounds, trace to roof
			if (ray.dir.z > 0)
				return false;
			Hit hit = tracePlane(land.maxBound, Vec!3(0, 0, 1), ray);
			assert(hit.hit);
			ray.org = hit.pos;
			currentBin = land.getBin(ray.org);
			assert(!currentBin.isNull);
			return true;
		} else if (ray.org.z < land.minBound.z) { // Outside z bounds, trace to floor
			if (ray.dir.z < 0)
				return false;
			Hit hit = tracePlane(land.minBound, Vec!3(0, 0, -1), ray);
			assert(hit.hit);
			ray.org = hit.pos;
			currentBin = land.getBin(ray.org);
			assert(!currentBin.isNull);
			return true;
		} else { // Initialize current bin
			currentBin = land.getBin(ray.org);
			assert(!currentBin.isNull);
			return true;
		}
	} else {
		Hit[2] hits;
		foreach (ubyte axis; 0 .. 2) { // Calculate distances to potential next bins walls
			if (ray.dir[axis] == 0) {
				hits[axis] = Hit();
			} else {
				float offset = land.binWidth[axis] * (currentBin.get()[axis] + (ray.dir[axis] > 0 ? 1
						: 0));
				Vec!3 point = land.minBound;
				point[axis] += offset;
				Vec!3 normal = Vec!3(0);
				normal[axis] = (ray.dir[axis] > 0) ? -1 : 1;
				hits[axis] = tracePlane(point, normal, ray);
				assert(hits[axis].hit);
				assert(hits[axis].t > 0);
			}
		}
		// Choose closest bin
		if (hits[0].t == hits[1].t) { // perfect slanted progression, should (almost) never happen
			if (hits[0].t == float.infinity)
				return false;
			currentBin.get()[0] += ((ray.dir[0] > 0) ? 1 : -1);
			currentBin.get()[1] += ((ray.dir[1] > 0) ? 1 : -1);
			ray.org = hits[0].pos;
		} else {
			ubyte minAxis = (hits[1].t > hits[0].t) ? 0 : 1;
			currentBin.get()[minAxis] += (ray.dir[minAxis] > 0) ? 1 : -1; // increment/decrement bin
			ray.org = hits[minAxis].pos;
		}
		// Check if next bin exits bounds
		// Note unsigned int underflow (< 0) wraps around to > binCount.
		if (currentBin.get()[0] >= land.binCount || currentBin.get()[1] >= land.binCount)
			return false;
		return true;
	}
}

unittest {
	Landscape land = Landscape([]);
	land.createBins(Vec!3(-1), Vec!3(1), 2);
	Nullable!(Vec!(2, uint)) bin;

	// From outside
	Ray ray = Ray(Vec!3(0.5, 0.5, 2), Vec!3(0, 0, -1));
	bool found = findNextBin(land, ray, bin);
	assert(found);
	assert(ray.org == Vec!3(0.5, 0.5, 1));
	assert(bin.get() == [1, 1]);
	// Exit
	found = findNextBin(land, ray, bin);
	assert(!found);

	// From inside
	bin = bin.init;
	ray = Ray(Vec!3(0.5, 0.5, 0.5), Vec!3(0, 0, -1));
	found = findNextBin(land, ray, bin);
	assert(found);
	assert(ray.org == Vec!3(0.5, 0.5, 0.5));
	assert(bin.get() == [1, 1]);

	// Go through
	bin = bin.init;
	Vec!3 start = Vec!3(-0.5, -0.25, -1);
	Vec!3 stop = Vec!3(0.5, 1.0, 0.5);
	Vec!3 dir = (stop - start).normalize();
	ray = Ray(start - dir, dir);
	found = findNextBin(land, ray, bin);
	assert(found);
	assert(ray.org.almostEquals(start));
	assert(bin.get() == [0, 0]);
	found = findNextBin(land, ray, bin);
	assert(found);
	assert(bin.get() == [0, 1]);
	found = findNextBin(land, ray, bin);
	assert(found);
	assert(bin.get() == [1, 1]);
	found = findNextBin(land, ray, bin);
	assert(!found);
	assert(ray.org.almostEquals(stop));
}

Hit trace(const Landscape land, Ray ray) {
	if (!land.useBins) {
		return trace(land.peaks, ray);
	} else {
		Nullable!(Vec!(2, uint)) currentBin;
		while (findNextBin(land, ray, currentBin)) {
			const Vec!3[] peaks = land.getBin(currentBin.get().x, currentBin.get().y);
			Hit hit = trace(peaks, ray);
			if (hit.hit)
				return hit;
		}
		return Hit();
	}
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

/// Reflects dir using the normal.
/// Params:
///   dir = normalized in direction (pointing into the surface)
///   normal = normalized normal
/// Returns: normalized reflected direction
Vec!3 reflect(Vec!3 dir, Vec!3 normal) {
	return dir - (normal * (2 * dir.dot(normal)));
}

unittest {
	Vec!3 dir = Vec!3(0.5, 0, -0.5).normalize();
	Vec!3 normal = Vec!3(0, 0, 1);
	Vec!3 reflected = reflect(dir, normal);
	reflected.assertAlmostEquals(Vec!3(0.5, 0, 0.5).normalize());
}

struct ReflectData {
	bool exits;
	Vec!3 outDir;
	uint reflectCount;
	uint reflectID; // gives what normals were hit. (E,N,W,S)=(0,1,2,3) tree id.
}

ReflectData reflectRecurse(Landscape land, Ray ray, uint reflectCount) {
	uint reflectID = 0;
	foreach (r; 0 .. reflectCount + 1) { //TODO: prefer not to use +1
		Hit hit = trace(land, ray);
		if (!hit.hit) {
			assert(r > 0);
			return ReflectData(true, ray.dir, r, reflectID);
		}
		if (r > 0)
			reflectID = 4 * reflectID + 4 + hit.face;
		else
			reflectID = hit.face;
		ray.dir = reflect(ray.dir, shape.normals[hit.face]);
		ray.org = hit.pos;
		ray.exculdePeak = hit.peakID;
	}
	return ReflectData(false, ray.dir, reflectCount);
}

unittest {
	Landscape land = [Vec!3(-2, 0, 1), Vec!3(2, 0, 1)];
	shape = PyramidShape(degreesToRadians(45.0));

	Vec!3 org = Vec!3(-1, 0, 2);
	Vec!3 dir = Vec!3(0, 0, -1);
	Ray ray = Ray(org, dir);
	ReflectData reflectData = reflectRecurse(land, ray, 2);
	assert(reflectData.exits);
	assert(reflectData.reflectCount == 2);
	assert(reflectData.reflectID == 4 + 2);
	reflectData.outDir.assertAlmostEquals(Vec!3(0, 0, 1));
}

string reflectIDToString(uint reflectID) {
	string reflectPath = "";
	immutable char[4] faceNames = ['E', 'N', 'W', 'S'];

	while (true) {
		reflectPath = faceNames[reflectID % 4] ~ reflectPath;
		if (reflectID < 4)
			break;
		reflectID = (reflectID / 4) - 1;
	}
	return reflectPath;
}

unittest {
	string str1 = reflectIDToString(0);
	assert(str1 == "E");
	string str2 = reflectIDToString(70);
	assert(str2 == "SEW");
}
