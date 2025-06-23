import csv;
import distribution;
import landscape;
import ply;
import pyramid_shape;
import ray;
import sphere_samplers;
import tracing;
import vdmath;
import vdmath.misc;

import std.algorithm.comparison : clamp;
import std.algorithm.searching : canFind;
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

// immutable PyramidShape shape = PyramidShape(degreesToRadians(_slope));
enum _slope = degreesToRadians(54.7);
enum _width = 20.0; //? in micron
enum float width = 40; // used
enum _density = 0.6; //? per micron²

uint simpleSphereLatCount = 16;
uint simpleSphereLongCount = 3; //10;
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
	case "measureBackBounceProb":
		return measureBackBounceProb(args[2 .. 5]);
	case "reflectDist":
		return measureReflectPathDist(args[2 .. $]);
	case "experiment":
		return experiment();
	case "generate":
		return generate();
	case "generateVar":
		return generateVar();
	case "test":
		return test_shadowing();
	case "measure":
		return measure_shadowing();
	default:
		return writeln("Choose:\n\t- reflectDist\n\t- experiment\n\t- generate\n\t- test");
	}
}

uint pow4(uint exponent) pure {
	return 1u << (2u * exponent);
}

///See_Also: reflectIDCount (identical)
uint pow4sum(uint exponent) pure {
	return (4 * (pow4(exponent) - 1)) / 3; // Modified geometric series
}

unittest {
	assert(pow4sum(0u) == 0);
	assert(pow4sum(1u) == 4);
	assert(pow4sum(2u) == 20);
	assert(pow4sum(3u) == 84);
	assert(pow4sum(4u) == 340);
}

// unittest {
// 	Vec!3[] peaks = [
// 		Vec!3(-1, -1, 0),
// 		Vec!3(1, -1, 0),
// 		Vec!3(1, 1, 0),
// 		Vec!3(-1, 1, 0)
// 	];
// 	Landscape land = Landscape(2, peaks, PI_4);
// 	Vec!3[] samples = Landscape.createPoints()
// }

import plt = matplotlibd.pyplot;

void measureBackBounceProb(string[] args) {
	float slope = degreesToRadians(54.7);
	float density = 0.6;

	float z = args[0].to!float; // 0 to -L
	float beta = args[1].to!float; // -1 to 1
	float x = -z / tan(slope);
	float y = beta * x;

	// float theta = degreesToRadians(args[2].to!float); x-axis
	float phi = degreesToRadians(args[2].to!float); // 0 to 90
	// Vec!3 c = Vec!3(cos(theta) * cos(phi), cos(theta) * sin(phi), sin(theta));

	Vec!3 p = Vec!3(x, y, z);
	float width = 2000;
	Landscape land = Landscape(width, density, slope);
	land.createBins(25, 0, -land.approxHeight, true);
	uint samples = 500;

	Vec!3[] sampleStart;
	uint[] sampleIndex;
	sampleStart.reserve(samples);
	foreach (i; 0 .. land.peaks.length) {
		Vec!3 peak = land.peaks[i];
		if (peak.x > width * 2 / 3.0 || peak.y > width * 2 / 3.0)
			continue;

		Vec!3 start = peak + p;
		Ray testRay = Ray(start, Vec!3(0, 0, 1));
		Hit hit = traceLand(land, testRay);
		if (!hit.hit) {
			sampleStart ~= start;
			sampleIndex ~= cast(uint) i;
			if (sampleStart.length == samples)
				break;
		}
	}
	assert(sampleStart.length == samples, "Could not find enough peaks");

	enum thetaSamples = 50;
	float[] thetas = new float[thetaSamples];
	float[] rebounceCounts = new float[thetaSamples];
	foreach (i; 0 .. thetaSamples)
		thetas[i] = i * PI_2 / thetaSamples;
	foreach (index, theta; thetas) {
		write("\r", index, " / ", thetaSamples, "\t\t");
		Vec!3 c = Vec!3(cos(theta) * cos(phi), cos(theta) * sin(phi), sin(theta));
		uint rebounces = 0; // bounces back and hits previous pyramid
		foreach (startIndex, start; sampleStart) {
			Ray ray = Ray(start, start - c);
			ReflectData!true reflectData = reflectRecurse!true(land, ray, 3);
			if (reflectData.history.length <= 2 || reflectData.history[0] != sampleIndex[startIndex]) // no bounce back
				continue;
			if (reflectData.history[0] == reflectData.history[2])
				rebounces += 1;
		}
		rebounceCounts[index] = rebounces.to!float;
	}
	writeln();

	plt.plot(thetas, rebounceCounts, "r-");
	plt.savefig("rebounceProb.png");
	plt.show();
}

void measureReflectPathDist(string[] args) {
	rndGen().seed(1);

	bool gridSample = false;
	bool gridLand = false;

	float width = 500;
	const uint binCount = 25;
	Landscape land;
	if (gridLand)
		land = Landscape(width, Landscape.createGrid(width, _density, 0));
	else
		land = Landscape(width, _density);
	// land = Landscape(width, [Vec!3(0)], PI_4);
	// land.save("reflectDist.ply");
	land.createBins(binCount, 0, -land.approxHeight, true);

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
	uint[] reBouncedCounts = new uint[optionCount];
	uint captureCount = 0;

	Vec!3[] offsets;
	if (gridSample) {
		uint root = cast(uint) sqrt(cast(float) sampleCount);
		enforce(sqrt(cast(float) sampleCount) == root); // whole root
		offsets = Landscape.createGrid(width * 2 / 3.0, root, 0);
		// Mat!4 = rotZ = []
	} else
		offsets = Landscape.createPoints(width * 2 / 3.0, sampleCount, 0);

	CSV csv = CSV(',', false, "index", "indexStr", "count", "prob", "rebounceProb");

	foreach (i; 0 .. sampleCount) {
		debug writef("\rRaytracing: %u/%u", i, sampleCount);
		float minHeight = land.approxHeight / abs(dir.z);

		Vec!3 target = offsets[i];
		Vec!3 cam = (target - dir * 2 * minHeight);

		Ray ray = Ray(cam, dir);
		ReflectData!true reflectData = reflectRecurse!true(land, ray, reflectCount);

		if (!reflectData.hitSurface()) // ignore in release mode
			continue; // TODO: could invalidate & remove from sampleCount in final results
		else if (!reflectData.exits) {
			captureCount += 1;
			continue;
		} else {
			pathCounts[reflectData.reflectID] += 1;
			bool reBounced = (reflectData.history.length >= 2 && reflectData.history[$ - 2] == reflectData
					.history[$ - 1]);
			if (reBounced)
				reBouncedCounts[reflectData.reflectID] += 1;
		}
	}
	debug writefln("\rRaytracing: %u/%u", sampleCount, sampleCount);

	foreach (uint i; 0 .. optionCount) {
		csv.addEntry(i, reflectIDToString(i), pathCounts[i], pathCounts[i].to!float / sampleCount, reBouncedCounts[i]
				.to!float / sampleCount);
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

		ReflectData!false reflectData = reflectRecurse!false(land, ray, reflectCount);
		data[reflectData.reflectID] += reflectData.exits;
		//TODO: no clue what to do when !reflectData.hitSurface() in release mode
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
	Landscape land = Landscape(scaling * size, _density / (scaling ^^ 2));
	return land;
}

void generateVar() {
	write("landscape size (um): ");
	float size = readln().strip().to!float;
	enforce(size > 0, "Size cannot be negative.");

	write("output file: ");
	string fileName = readln().strip();
	enforce(isValidFilename(fileName), "File name invalid.");

	string heightsFile = "heightCumulative.csv"; // Source: "Opto-electrical modelling and optimization study of a novel IBC c-Si solar cellOpto-electrical modelling and optimization study of a novel IBC c-Si solar cell"
	Distribution heightDistribution = new HistogramDistribution(heightsFile, maxVarHeight);
	Landscape land = Landscape(size, _density, _slope, heightDistribution); // 20 micron, 0.6 per micron²
	land.approxHeight = maxVarHeight;
	land.save(fileName);
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
	PyramidShape shape = PyramidShape(_slope);
	float L = scaling * ceil(tan(shape.slope) * sqrt(-log(error) / (4.0 * _density))); // ceil optional

	Landscape land = Landscape(size, _density);

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

	PyramidShape shape = PyramidShape(_slope);

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
	Landscape land = Landscape(width, _density, _slope, heightDistribution); // 20 micron, 0.6 per micron²
	Landscape landConstant = Landscape(width, _density, _slope, heightDistributionConstant);

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

	// measureLand(landConstant, "results/land.csv");
	measureLand(landConstant, "results/landVar.csv");
	// measure(land, sphereSamples, maxVarHeight, "hits.csv");
	measure(land, sphereSampler, maxVarHeight, heightBins, "results/hitsVar.csv",
		"results/hitsHeightVar.csv", "results/heightDistributionVar.csv",
		"results/hitsVarByHeight.csv");
	measure(landConstant, sphereSampler, maxConstHeight, heightBins, "results/hitsConstant.csv",
		"results/hitsHeightConstant.csv", "results/heightDistributionConstant.csv",
		"results/hitsConstantByHeight.csv");

	// measureHighPoint(landConstant, sphereSampler.data, "results/highPointSampling.csv");
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
