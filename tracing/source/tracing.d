module tracing;

import landscape;
import pyramid_shape;
import ray;
import vdmath;
import vdmath.misc;

import std.conv : to;
import std.math.algebraic : abs;
import std.math.constants : PI_4;
import std.math.traits : signbit;
import std.parallelism : parallel;
import std.typecons : Nullable;

bool findNextBin(const Landscape land, ref Ray ray, ref Nullable!(Vec!(2, uint)) currentBin) {
	assert(land.useBins);
	if (currentBin.isNull) {
		foreach (ubyte axis; 0 .. 3) {
			bool tooMuch = ray.org[axis] > land.maxBound[axis];
			bool tooLittle = ray.org[axis] < land.minBound[axis];
			Vec!3 pointOnPlane;
			Vec!3 dir = Vec!3(0);
			if (tooMuch) {
				if (ray.dir[axis] > 0)
					return false;
				pointOnPlane = land.maxBound;
				dir[axis] = -1;
			} else if (tooLittle) {
				if (ray.dir[axis] < 0)
					return false;
				pointOnPlane = land.minBound;
				dir[axis] = 1;
			} else
				continue;

			Hit hit = tracePlane(pointOnPlane, dir, ray);
			assert(hit.hit);
			ray.org = hit.pos;
		}
		currentBin = land.getBin(ray.org);
		return !currentBin.isNull;
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
	Landscape land = Landscape(2, [Vec!3(0)]);
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
		Hit hit = traceFace(peak, f, i, ray);
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

Hit traceFace(Vec!3 peak, Face face, uint faceID, Ray ray) {
	Hit hit = tracePlane(peak, face.normal, ray);
	if (!hit.hit)
		return hit;
	Vec!(3, double) P = (cast(Vec!(3, double)) hit.pos) - peak; // float error mitigation

	final switch (faceID) {
	case 0:
		hit.hit = P.x >= 0 && abs(P.y) <= abs(P.x);
		break;
	case 1:
		hit.hit = P.y >= 0 && abs(P.x) <= abs(P.y);
		break;
	case 2:
		hit.hit = P.x <= 0 && abs(P.y) <= abs(P.x);
		break;
	case 3:
		hit.hit = P.y <= 0 && abs(P.x) <= abs(P.y);
		break;
	}

	return hit;
}

unittest {
	PyramidShape shape = PyramidShape(PI_4);
	Vec!3 peak = Vec!3(0);
	Face face = shape.E;
	Ray ray = Ray(Vec!3(0.1, 0.09999, 0), Vec!3(0, 0, -1));
	Hit hit = traceFace(peak, face, 0, ray);
	assert(hit.hit);
	assert(hit.face == 0);
	assert(hit.pos.almostEquals(Vec!3(0.1, 0.09999, -0.1)));

	ray = Ray(Vec!3(0.1, 0.1001, 0), Vec!3(0, 0, -1));
	hit = traceFace(peak, face, 0, ray);
	assert(!hit.hit);
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
/// Returns: reflected direction
Vec!3 reflect(Vec!3 dir, Vec!3 normal) {
	return dir - (normal * (2 * dir.dot(normal)));
}

unittest {
	Vec!3 dir = Vec!3(0.5, 0, -0.5).normalize();
	Vec!3 normal = Vec!3(0, 0, 1);
	Vec!3 reflected = reflect(dir, normal);
	reflected.assertAlmostEquals(Vec!3(0.5, 0, 0.5).normalize());
}

auto reflectRecurse(bool trackHistory = false)(const Landscape land, Ray ray, uint reflectCount) {
	uint reflectID = 0;
	static if (trackHistory)
		uint[] history;
	foreach (r; 0 .. reflectCount + 1) { //TODO: prefer not to use +1
		Hit hit = traceLand(land, ray);
		if (!hit.hit) {
			assert(r > 0, "Initial ray missed landscape! " ~ ray.toString());
			static if (trackHistory)
				return ReflectData!true(true, ray, r, reflectID, history);
			else
				return ReflectData!false(true, ray, r, reflectID);
		}
		if (r > 0) {
			reflectID = 4 * reflectID + 4 + hit.face;
		} else
			reflectID = hit.face;
		static if (trackHistory)
			history ~= hit.peakID;
		ray.dir = reflect(ray.dir, land.shape.normals[hit.face]);
		ray.org = hit.pos;
		ray.excludePeak = hit.peakID;
	}
	static if (trackHistory)
		return ReflectData!true(false, ray, reflectCount, reflectID, history);
	else
		return ReflectData!false(false, ray, reflectCount, reflectID);
}

unittest {
	Landscape land = Landscape(5, [Vec!3(-2, 0, 1), Vec!3(2, 0, 1)], PI_4);
	Vec!3 org = Vec!3(-1, 0, 2);
	Vec!3 dir = Vec!3(0, 0, -1);
	Ray ray = Ray(org, dir);
	ReflectData reflectData = reflectRecurse(land, ray, 2);
	assert(reflectData.exits);
	assert(reflectData.reflectCount == 2);
	assert(reflectData.reflectID == 4 + 2);
	reflectData.outRay.dir.assertAlmostEquals(Vec!3(0, 0, 1));
}

///See_Also: pow4sum (identical)
uint reflectIDCount(uint pathLength) pure {
	uint count = 0;
	foreach (i; 1 .. pathLength + 1)
		count += 4 ^^ i;
	return count;
}

unittest {
	assert(reflectIDCount(0) == 0);
	assert(reflectIDCount(1) == 4);
	assert(reflectIDCount(2) == 20);
	assert(reflectIDCount(3) == 84);
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
	Landscape land = Landscape(20, [Vec!3(0, 0, 0), Vec!3(10, 0, 0)], PI_4);
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

// Correct Faces
unittest {
	Landscape land = Landscape(4, [Vec!3(0)], PI_4);
	Vec!3[] offsets = [
		Vec!3(1, 0, 0),
		Vec!3(0, 1, 0),
		Vec!3(-1, 0, 0),
		Vec!3(0, -1, 0)
	];
	ubyte[] faces = [0, 1, 2, 3];
	foreach (i, offset; offsets) {
		Ray ray = Ray(offset, Vec!3(0, 0, -1));
		Hit hit = traceLand(land, ray);
		assert(hit.hit);
		assert(hit.t == 1.0);
		assert(hit.pos == offset + Vec!3(0, 0, -1));
		assert(hit.face == faces[i]);
		assert(hit.peakID == 0);
	}
}

// Correct Faces
unittest {
	Landscape land = Landscape(4, [Vec!3(0)], PI_4);
	Vec!3[] offsets = Landscape.createGrid(4, 100u, 0, true);

	byte getFace(Vec!3 offset) {
		if (abs(offset.x) == abs(offset.y))
			return -1;
		if (offset.x > 0 && offset.x > abs(offset.y))
			return 0;
		if (offset.x < 0 && -offset.x > abs(offset.y))
			return 2;
		if (offset.y > 0 && offset.y > abs(offset.x))
			return 1;
		if (offset.y < 0 && -offset.y > abs(offset.x))
			return 3;
		assert(0);
	}

	import std.algorithm.comparison : max;

	foreach (i, offset; offsets) {
		byte expectedFace = getFace(offset);
		if (expectedFace == -1)
			continue;
		float distance = max(abs(offset.x), abs(offset.y));
	rerun:
		Ray ray = Ray(offset, Vec!3(0, 0, -1));
		Hit hit = traceLand(land, ray);
		assert(hit.hit);
		assert(abs(hit.t - distance) < 1e-5);
		assert(hit.pos.almostEquals(offset + Vec!3(0, 0, -distance)));
		bool fail = hit.face != expectedFace;
		if (fail)
			goto rerun;
		assertEqual(hit.face, expectedFace);
		assert(hit.peakID == 0);
	}
}
