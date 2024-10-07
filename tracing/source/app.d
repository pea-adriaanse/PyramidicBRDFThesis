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

import landscape;
import pyramid_shape;
import ply;
import distribution;

immutable PyramidShape shape = PyramidShape(degreesToRadians(_slope));
enum _slope = 45.0; //54.7;
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

	float width = 20;
	float error = 0.001;
	const uint binCount = 25;
	float L = ceil(tan(shape.slope) * sqrt(-log(error) / (4.0 * _density))); // ceil optional
	Distribution heightDistribution = new ConstantDistribution(L);
	Vec!3[] peaks = createPeaks(width, _density, heightDistribution);
	Landscape land = Landscape(peaks, &shape);
	land.createBins(Vec!3(-width - L, -width - L, -L - 1), Vec!3(width + L, width + L, L + 1), binCount);

	land.save("reflectDist.ply", L);

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

	Vec!2[] offsets;

	if (sampleGrid) {
		uint root = sampleCount;
		sampleCount = sampleCount ^^ 2;

		offsets = new Vec!2[sampleCount];

		float step = width / root;
		foreach (x; 0 .. root) {
			foreach (y; 0 .. root) {
				uint index = x * root + y;
				offsets[index] = Vec!2(x * step - width / 2.0, y * step - width / 2.0);
			}
		}
	} else {
		offsets = new Vec!2[sampleCount];
		foreach (i; 0 .. sampleCount) {
			float x = uniform(-width / 3.0, width / 3.0);
			float y = uniform(-width / 3.0, width / 3.0);
			offsets[i] = Vec!2(x, y);
		}
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
	Landscape land = Landscape(createPeaks(scaling * size, _density / (scaling ^^ 2), heightDistribution), &shape);
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
	Landscape land = Landscape(createPeaks(scaling * size, _density / (scaling ^^ 2), heightDistribution), &shape);

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
	Landscape land = Landscape(createPeaks(width, _density, heightDistribution), &shape); // 20 micron, 0.6 per micron²
	Landscape landConstant = Landscape(createPeaks(width, _density, heightDistributionConstant), &shape);

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
		Hit hit = traceLand(land, ray);
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
			Hit hit = traceLand(land, ray);

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
			Hit hit = traceLand(land, ray);
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
			hit = traceLand(land, ray);
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
			Hit hit = traceLand(land, ray);

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

struct Ray {
	Vec!3 org;
	Vec!3 dir;
	uint excludePeak = uint.max;
}

struct Hit {
	bool hit = false;
	Vec!3 pos;
	float t = float.infinity;
	uint peakID;
	ubyte face;
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
				float offset = land.binWidth[axis] * (
					currentBin.get()[axis] + (ray.dir[axis] > 0 ? 1 : 0));
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
	Landscape land = Landscape([], &shape);
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

Hit traceLand(const Landscape land, Ray ray) {
	if (!land.useBins) {
		return tracePeaks(land.shape, land.peaks, ray);
	} else {
		Nullable!(Vec!(2, uint)) currentBin;
		while (findNextBin(land, ray, currentBin)) {
			const Vec!3[] peaks = land.getBin(currentBin.get().x, currentBin.get().y);
			Hit hit = tracePeaks(land.shape, peaks, ray);
			if (hit.hit)
				return hit;
		}
		return Hit();
	}
}

Hit tracePeaks(const ref PyramidShape shape, const Vec!3[] peaks, Ray ray) {
	shared Hit minHit;
	foreach (id, peak; parallel(peaks)) {
		// foreach (id, peak; peaks) {
		if (ray.excludePeak == id)
			continue;
		Hit hit = tracePeak(shape, peak, ray);
		if (hit.hit) {
			hit.peakID = cast(uint) id;
			synchronized if (hit.t < minHit.t)
				minHit = cast(shared Hit) hit;
		}
	}
	return minHit;
}

Hit tracePeak(const ref PyramidShape shape, Vec!3 peak, Ray ray) {
	Hit minHit;
	foreach (ubyte i, Face f; shape.faces) {
		Hit hit = traceFace(shape, peak, f, ray);
		if (hit.hit && hit.t < minHit.t) {
			minHit = hit;
			minHit.face = i;
		}
	}
	return minHit;
}

unittest {
	PyramidShape shape = PyramidShape(degreesToRadians(45.0));
	Hit hit = tracePeak(shape, Vec!3(0, 0.25, 0), Ray(Vec!3(0, -1, 1), Vec!3(0, 0, -1)));
	assert(hit.hit);
	assert(hit.pos == Vec!3(0, -1, -1.25));
}

unittest {
	PyramidShape shape = PyramidShape(degreesToRadians(54.7));
	Hit hit = tracePeak(shape, Vec!3(0, 0, 1), Ray(Vec!3(0), Vec!3(0, 0, 1)));
	assert(hit.hit);
}

Hit traceFace(const ref PyramidShape shape, Vec!3 peak, Face face, Ray ray) {
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
	Ray outRay;
	uint reflectCount;
	uint reflectID; // gives what normals were hit. (E,N,W,S)=(0,1,2,3) tree id.
}

ReflectData reflectRecurse(Landscape land, Ray ray, uint reflectCount) {
	uint reflectID = 0;
	foreach (r; 0 .. reflectCount + 1) { //TODO: prefer not to use +1
		Hit hit = traceLand(land, ray);
		if (!hit.hit) {
			assert(r > 0, "Initial ray missed landscape! " ~ ray.to!string);
			return ReflectData(true, ray, r, reflectID);
		}
		if (r > 0)
			reflectID = 4 * reflectID + 4 + hit.face;
		else
			reflectID = hit.face;
		ray.dir = reflect(ray.dir, shape.normals[hit.face]);
		ray.org = hit.pos;
		ray.excludePeak = hit.peakID;
	}
	return ReflectData(false, ray, reflectCount);
}

unittest {
	Landscape land = Landscape([Vec!3(-2, 0, 1), Vec!3(2, 0, 1)], PI_4);
	Vec!3 org = Vec!3(-1, 0, 2);
	Vec!3 dir = Vec!3(0, 0, -1);
	Ray ray = Ray(org, dir);
	ReflectData reflectData = reflectRecurse(land, ray, 2);
	assert(reflectData.exits);
	assert(reflectData.reflectCount == 2);
	assert(reflectData.reflectID == 4 + 2);
	reflectData.outRay.dir.assertAlmostEquals(Vec!3(0, 0, 1));
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

uint reflectStringToID(string path) {
	uint id = 0;

	ubyte faceID(char c) {
		final switch (c) {
		case 'E':
			return 0;
		case 'N':
			return 1;
		case 'W':
			return 2;
		case 'S':
			return 3;
		}
	}

	foreach (i, char c; path) {
		if (i > 0)
			id = (id + 1) * 4;
		id += faceID(c);
	}
	return id;
}

unittest {
	uint id1 = reflectStringToID("E");
	assert(id1 == 0);
	uint id2 = reflectStringToID("SEW");
	assert(id2 == 70);
}

// 45 degrees triangle exit directions
unittest {
	PyramidShape shape = PyramidShape(PI_4);
	Landscape land = Landscape([Vec!3(0, 0, 0), Vec!3(10, 0, 0)], &shape);
	Ray ray = Ray(Vec!3(1, 0, 1), Vec!3(0, 0, -1));
	ReflectData reflect = reflectRecurse(land, ray, 2);
	assert(reflect.exits);
	assert(reflect.reflectCount == 2);
	assert(reflect.outRay.dir.almostEquals(Vec!3(0, 0, 1)));
	assert(reflect.outRay.org.almostEquals(Vec!3(9, 0, -1)));
	assert(reflect.reflectID == reflectStringToID("EW"));

	Ray rayBack = Ray(Vec!3(9, 0, 1), Vec!3(0, 0, -1));
	ReflectData reflectBack = reflectRecurse(land, rayBack, 2);
	assert(reflectBack.exits);
	assert(reflectBack.reflectCount == 2);
	assert(reflectBack.outRay.dir.almostEquals(Vec!3(0, 0, 1)));
	assert(reflectBack.outRay.org.almostEquals(Vec!3(1, 0, -1)));
	assert(reflectBack.reflectID == reflectStringToID("WE"));
}
