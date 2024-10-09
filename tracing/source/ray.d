module ray;
import vdmath;

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

struct ReflectData {
	bool exits;
	Ray outRay;
	uint reflectCount;
	uint reflectID; // gives what normals were hit. (E,N,W,S)=(0,1,2,3) tree id.
}
