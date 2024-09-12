module sphere_samplers;

import vdmath;
import std.math : sqrt, PI, asin, sin, cos;
import std.conv : to;
import std.stdio : File;

abstract class SphereSampler {
	Vec!3[] data;
	Vec!2[] angles;

	string name();
	void save(string fileName) {
		File file = File(fileName, "w");
		file.writeln("x,y,z,phi,lat");
		foreach (i; 0 .. data.length) {
			file.writeln(data[i].x.to!string, ",", data[i].y.to!string, ",", data[i].z.to!string,
				",", angles[i][0].to!string, ",", angles[i][1].to!string);
		}
	}
}

Vec!3 angularToXYZ(float longitude, float latitude) {
	float x = cos(longitude) * cos(latitude);
	float y = sin(longitude) * cos(latitude);
	float z = sin(latitude);
	return Vec!3(x, y, z);
}

class SphereFibonacci : SphereSampler {
	static float two_pi_phi_inv = 2.0 * PI * 2.0 / (1.0 + sqrt(5.0));
	static float lat_factor;

	this(uint N) {
		this.lat_factor = 2.0 / (2.0 * N + 1);
		data.reserve(N);
		angles.reserve(N);
		foreach (i; 0 .. N) {
			float longitude = i * two_pi_phi_inv;
			float latitude = asin(i * lat_factor);
			data ~= angularToXYZ(longitude, latitude);
			angles ~= Vec!2(longitude, latitude);
		}
	}

	override string name() {
		return "fibonacci";
	}
}

class SphereSimple : SphereSampler {
	bool quadrant;

	this(uint Nphi, uint Nlat, bool quadrant) {
		assert(Nphi > 0 && Nlat > 0);
		this.quadrant = quadrant;
		data.reserve(Nphi * Nlat);
		angles.reserve(Nphi * Nlat);

		foreach (i; 0 .. Nphi) {
			float longitude = i * 2.0 * PI / (Nphi - 1);
			if (Nphi == 1)
				longitude = 0;
			if (quadrant)
				longitude /= 4;
			foreach (j; 0 .. Nlat) {
				float lat = j * PI / (2.0 * (Nlat - 1));
				if (Nlat == 1)
					lat = PI / 2.0;
				data ~= angularToXYZ(longitude, lat);
				angles ~= Vec!2(longitude, lat);
			}
		}
	}

	override string name() {
		if (quadrant)
			return "simple_quad";
		else
			return "simple";
	}
}
