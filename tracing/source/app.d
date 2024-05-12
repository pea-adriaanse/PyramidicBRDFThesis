import std.stdio;
import std.conv : to;
import std.random;
import std.math;
import math;
import misc;
import fibonacci;
import std.algorithm.searching : canFind;
import std.algorithm.sorting;
import std.conv : to;

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

// version = backfaceCulling;
enum _slope = 54.7;
enum _width = 20.0; //?
enum _density = 0.6; //?

PyramidShape shape = PyramidShape(degreesToRadians(_slope));

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

void main() {
	string heightsFile = "heightCumulative.csv"; // Source: "Opto-electrical modelling and optimization study of a novel IBC c-Si solar cellOpto-electrical modelling and optimization study of a novel IBC c-Si solar cell"
	// Distribution heightDistribution = new HistogramDistribution(heightsFile, 7.0);
	Distribution heightDistribution = new ConstantDistribution(1);
	float width = 100;
	Landscape land = createLandscape(width, _density, heightDistribution); // 20 micron, 0.6 per micron²

	// File landFile = File("land.csv", "w");
	// landFile.writeln("sep=,");
	// foreach (Vec!3 peak; land) {
	// 	landFile.writeln(peak.x.toStr, ",", peak.y.toStr, ",", peak.z.toStr);
	// }

	uint sphereSampleCount = 1000;
	SpherePoints spherePoints = SpherePoints(sphereSampleCount);
	uint peakSampleCount = 300;
	assert(peakSampleCount >= land.length);

	ulong[] peakSamples;
	peakSamples.reserve(peakSampleCount);
	while (peakSamples.length < peakSampleCount) {
		ulong sample = cast(ulong) uniform(0, land.length);
		if (!peakSamples.canFind(sample))
			peakSamples ~= sample;
	}

	// File results = File("results.csv", "w");
	// results.writeln("sep=,");
	// results.writeln("org_id,org_x,org_y,org_z,choice_id,peak_id,hit_x,hit_y,hit_z,peakSampleCount:,",
	// 	peakSampleCount.to!string, ",", "sphereSampleCount:", ",", sphereSampleCount.to!string);

	// foreach (i; 0 .. N * N) {
	// Vec!3 org = Vec!3((_width / N) * (i % N) - _width / 2, (_width / N) * (floor((cast(float) i) / N)) - _width / 2,
	// 	8.0);
	// Ray ray = Ray(org, Vec!3(0, 0, -1));
	// if (!hit.hit)
	// 	hit.pos = Vec!3(0, 0, 1); // debugging

	ulong[] hits = new ulong[sphereSampleCount];
	Vec!3[] sphereSamples = new Vec!3[sphereSampleCount];

	foreach (i; 0 .. sphereSampleCount)
		sphereSamples[i] = spherePoints.getPoint(i) * cast(float) _width;

	foreach (j, p; peakSamples) {
		write("[", j.to!string, "/", peakSampleCount.to!string, "]\t\t\r");
		stdout.flush();
		Vec!3 peakSample = land[p];
		foreach (i, sphereSample; sphereSamples) {
			Vec!3 org = Vec!3(peakSample.x, peakSample.y, 0) + sphereSample;
			// Vec!3 target = Vec!3(peakSample.x, peakSample.y, 0); // TODO sample across one face!!
			Ray ray = Ray(org, -sphereSample);

			Hit hit = trace(land, ray);
			if (hit.peakID == p)
				hits[i] += 1;

			// results.writeln(i.to!string, ",", org.x.to!string, ",", org.y.to!string, ",",
			// 	org.z.to!string, ",", p.to!string, ",", hit.peakID.to!string, ",",
			// 	hit.pos.x.to!string, ",", hit.pos.y.to!string, ",", hit.pos.z.to!string);
		}
	}

	File settingsFile = File("settings.csv", "w");
	settingsFile.writeln("sphereSampleCount,",sphereSampleCount.to!string);
	settingsFile.writeln("peakSampleCount,",peakSampleCount.to!string);

	File hitsFile = File("hits.csv", "w");
	hitsFile.writeln("org_id,hits");
	foreach (i, h; hits) {
		hitsFile.writeln(i.to!string, ",", h.to!string);
	}
}

Landscape createLandscape(float width, float density, Distribution heightDistribution) {
	ulong count = cast(ulong)(width * width * density);
	Vec!3[] peaks;
	peaks.reserve(count);

	// Decide on heights first.
	float[] heights = new float[count];
	foreach (i; 0 .. count)
		heights[i] = heightDistribution.sample();

	// Sort to place largest first.
	sort!"a>b"(heights);

	while (peaks.length < count) {
		ulong i = peaks.length;
		Vec!3 newPeak;
		newPeak[0] = uniform!"()"(-width / 2, width / 2);
		newPeak[1] = uniform!"()"(-width / 2, width / 2);
		newPeak[2] = heights[i];

		// Cull if burried
		Hit hit = trace(peaks, Ray(newPeak, Vec!3(0, 0, 1)));
		if (!hit.hit) {
			peaks ~= newPeak;
		}
	}
	return peaks;
}

struct Ray {
	Vec!3 org;
	Vec!3 dir;
}

struct Hit {
	bool hit = false;
	Vec!3 pos;
	float t = float.infinity;
	ulong peakID;
}

Hit trace(Landscape land, Ray ray) {
	Hit minHit;
	foreach (id, Vec!3 peak; land) {
		Hit hit = trace(peak, ray);
		if (hit.hit && (hit.t < minHit.t)) {
			minHit = hit;
			minHit.peakID = id;
		}
	}
	return minHit;
}

Hit trace(Vec!3 peak, Ray ray) {
	foreach (Face f; shape.faces) {
		Hit hit = trace(peak, f, ray);
		if (hit.hit)
			return hit;
	}
	return Hit(false);
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
	if (slope == 0 || signbit(slope) != signbit(distToPlane))
		return Hit(false);

	float t = distToPlane / slope;
	Vec!3 pos = ray.org + ray.dir * t;
	return Hit(true, pos, t);
}
