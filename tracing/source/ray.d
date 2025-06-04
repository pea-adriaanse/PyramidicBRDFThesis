module ray;
import vdmath;

struct Ray {
	Vec!3 org;
	Vec!3 dir;
	uint excludePeak = uint.max;

	string toString() {
		import std.conv : to;

		return "Ray(org: " ~ org.toString() ~ ", dir: " ~ dir.toString() ~ ", excludePeak: " ~ excludePeak.to!string ~ ")";
	}
}

struct Hit {
	bool hit = false;
	Vec!3 pos;
	float t = float.infinity;
	uint peakID;
	ubyte face;
}

struct ReflectData(bool trackHistory = false) {
	bool exits;
	Ray outRay;
	uint reflectCount;
	uint reflectID; // gives what normals were hit. (E,N,W,S)=(0,1,2,3) tree id.
	static if (trackHistory)
		uint[] history;
}
