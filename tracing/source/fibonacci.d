module fibonacci;

import math.mat;
import std.math : sqrt, PI, asin, sin, cos;

float two_pi_phi_inv = 2.0 * PI * 2.0 / (1.0 + sqrt(5.0));

struct SpherePoints {
	uint N;
	float lat_factor;

	this(uint N) {
		this.N = N;
		this.lat_factor = 2.0 / (2.0 * N + 1);
	}

	Vec!3 getPoint(uint i) {
		float longitude = i * two_pi_phi_inv;
		float latitude = asin(i * lat_factor);
		float x = cos(longitude) * cos(latitude);
		float y = sin(longitude) * cos(latitude);
		float z = sin(latitude);
		return Vec!3(x, y, z);
	}
}
