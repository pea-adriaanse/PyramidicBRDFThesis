import pyd.pyd;
import pyd.embedded;
import pyd.extra;
import std.math.constants : PI, PI_2;
import std.range : iota;
import vdmath;
import std.math;

immutable script = import("plot.py");
immutable double alpha = 54.7 * PI / 180.0;

shared static this() {
	py_init();
}

struct ReflectContext {
	Vec!3 c1;
	Vec!3 c2;
	Vec!3 c3;

	Vec!3 b1;
	Vec!3 b2;

	this(double theta, double phi) {
		this.c1 = Vec!3(-sin(theta) * cos(phi), -sin(theta) * sin(phi), -cos(theta));
		this.c2 = Vec!3(cos(2 * alpha) * c1.x - sin(2 * alpha) * c1.z, c1.y, -sin(2 * alpha) * c1.x - cos(
				2 * alpha) * c1.z);
		this.c3 = Vec!3(cos(2 * alpha) * c2.x + sin(2 * alpha) * c2.z, c2.y, sin(2 * alpha) * c2.x - cos(
				2 * alpha) * c2.z);

		this.b1 = Vec!3(-c3.z + tan(alpha) * c3.y, -c3.z - tan(alpha) * c3.x, c3.y + c3.x);
		this.b2 = Vec!3(c3.z + tan(alpha) * c3.y, -c3.z - tan(alpha) * c3.x, c3.y - c3.x);

	}
}

double d1(double theta, double phi, ReflectContext context) {

}

double d2(double theta, double phi, ReflectContext context) {

}

double E(double theta, double phi) {
	ReflectContext context = ReflectContext(theta, phi);

}

void main() {
	auto context = new InterpContext();
	double xmin = 0;
	double xmax = PI_2;
	double ymin = -PI;
	double ymax = PI;
	context.xmin = xmin;
	context.xmax = xmax;
	context.ymin = ymin;
	context.ymax = ymax;

	double stepsize = 0.1;
	double[][] x;
	double[][] y;
	double[][] z;

	foreach (i; iota(xmin, xmax, stepsize)) {
		foreach (j; iota(ymin, ymax, stepsize)) {
			x[j][i] = i;
			y[j][i] = j;
			z[j][i] = E(i, j);
		}
	}

	context.x = d_to_python_numpy_ndarray(x);
	context.y = d_to_python_numpy_ndarray(y);
	context.z = d_to_python_numpy_ndarray(z);

	context.py_stmts(script);
}
