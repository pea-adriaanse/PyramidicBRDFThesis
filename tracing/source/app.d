import std.stdio;
import std.conv : to;
import std.random;
import std.math;
import math;
import misc;
import fibonacci;
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
enum _width = 20.0;
enum _density = 0.6;

PyramidShape shape = PyramidShape(degreesToRadians(_slope));

string toStr(float f) {
	return to!string(f);
}

void main() {
	Landscape land = createLandscape(_width, _density); // 20 micron, 0.6 per micron²

	File landFile = File("land.csv", "w");
	landFile.writeln("sep=,");
	foreach (Vec!3 peak; land) {
		landFile.writeln(peak.x.toStr, ",", peak.y.toStr, ",", peak.z.toStr);
	}

	uint N = 100;
	SpherePoints points = SpherePoints(N);

	File results = File("results.csv", "w");
	results.writeln("sep=,");
	results.writeln("org_x,org_y,org_z,hit,hit_x,hit_y,hit_z");

	foreach (i; 0 .. N * N) {
		Vec!3 org = Vec!3((_width / N) * (i % N) - _width / 2, (_width / N) * (floor((cast(float) i) / N)) - _width / 2,
			1); // debugging
		Ray ray = Ray(org, Vec!3(0, 0, -1));

		// Vec!3 org = points.getPoint(i);
		// Ray ray = Ray(org, Vec!3(0));
		Hit hit = trace(land, ray);

		if (!hit.hit)
			hit.pos = Vec!3(0, 0, 1); // debugging
		results.writeln(org.x.to!string, ",", org.y.to!string, ",", org.z.to!string, ",",
			hit.hit ? "1" : "0", ",", hit.pos.x.to!string, ",", hit.pos.y.to!string, ",", hit.pos.z.to!string);
	}

}

Landscape createLandscape(float width, float density) {
	ulong count = cast(ulong)(width * width * density);
	Vec!3[] peaks = new Vec!3[count];
	foreach (i; 0 .. count) {
		peaks[i][0] = uniform!"()"(-width / 2, width / 2);
		peaks[i][1] = uniform!"()"(-width / 2, width / 2);
		peaks[i][2] = 0;
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
}

Hit trace(Landscape land, Ray ray) {
	Hit minHit;
	foreach (Vec!3 peak; land) {
		Hit hit = trace(peak, ray);
		if (hit.hit && (hit.t < minHit.t))
			minHit = hit;
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

Hit trace(Vec!3 peak, Ray ray) {
	foreach (Face f; shape.faces) {
		Hit hit = trace(peak, f, ray);
		if (hit.hit)
			return hit;
	}
	return Hit(false);
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

	version (backfaceCulling) {
		if (slope > 0)
			return Hit(false);
	}
	if (slope == 0 || signbit(slope) != signbit(distToPlane))
		return Hit(false);

	float t = distToPlane / slope;
	Vec!3 pos = ray.org + ray.dir * t;
	return Hit(true, pos, t);
}
