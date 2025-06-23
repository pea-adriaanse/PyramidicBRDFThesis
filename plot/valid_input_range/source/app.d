import backscattering;
import integration;
import plotting;

import std.array : array;
import std.conv : to;
import std.stdio;
import std.math.constants : PI, PI_2;
import pyd.pyd;

immutable script = import("plot.py");

void main(string[] args) {
	Bounds!double thetaBounds = Bounds!double(0, PI_2, PI / 128);
	Bounds!double phiBounds = Bounds!double(0, PI, PI / 128);
	ZLimits zLimits = ZLimits(0, 1, false);

	py_init();
	plotFunction(thetaBounds, phiBounds, zLimits, &backBounce, "backBounce.bin", "");

	// double function(double, double) xFunc = (double d, double e) => d+e/4;
	// plotFunction(context, thetaBounds, phiBounds, zLimits, xFunc, "backBounce.bin", "");

	// import table;
	// BRDF brdf = new BRDF("backBounce.bin");
	// brdf.readTable(degreesToRadians(15.0), 0).writeln;

	// Bounds!double thetaBounds = Bounds!double(-PI_2, PI_2, PI / 32);
	// Bounds!double phiBounds = Bounds!double(-PI, PI, PI / 32);
	// plotFunction(context, thetaBounds, phiBounds, zLimits, &isValidAngle, "", "");

	// plotFunction(context, thetaBounds, phiBounds, zLimits, &backBounce, "", "");
	// plotFunction(context, thetaBounds, phiBounds, zLimits, &backBounceNormal, "", "");
}

// void main(string[] args) {
// 	double pyramidAngle = degreesToRadians!double(54.7);
// 	if (args.length >= 2)
// 		pyramidAngle = degreesToRadians!double(args[1].to!double);

// 	Vec!3 eastNormal = Vec!3(sin(pyramidAngle), 0, cos(pyramidAngle)).normalize();
// 	Vec!3 westNormal = Vec!3(-sin(pyramidAngle), 0, cos(pyramidAngle)).normalize();
// bool isValidBackbounce(double theta, double phi) {
// 	double x = sin(theta) * cos(phi);
// 	double y = sin(theta) * sin(phi);
// 	double z = cos(theta);
// 	Vec!3 c0 = (-Vec!3(x, y, z));

// 	if (c0.dot(eastNormal) >= 0) // does not hit east
// 		return false;

// 	Vec!3 c1 = reflect(c0, eastNormal);

// 	if (c1.dot(westNormal) >= 0) // does not hit west
// 		return false;

// 	Vec!3 c2 = reflect(c1, westNormal);

// 	if (c2.dot(eastNormal) >= 0) // does not hit east
// 		return false;

// 	if (c2.z <= 0) // does not escape on 2 bounce: is a precondition
// 		return false;
// 	return true;
// }

// 	auto context = new InterpContext();

// 	double[2] xBounds = [0, PI_2];
// 	double[2] yBounds = [-PI, PI];

// 	context.xBounds = xBounds;
// 	context.yBounds = yBounds;

// 	double stepsize = 0.001;

// 	double[] xCoords = iota(xBounds[0], xBounds[1], stepsize).array();
// 	double[] yCoords = iota(yBounds[0], yBounds[1], stepsize).array();

// 	double[][] xs = new double[][yCoords.length];
// 	double[][] ys = new double[][yCoords.length];
// 	double[][] zs = new double[][yCoords.length];

// 	foreach (yIndex, y; yCoords) {
// 		xs[yIndex] = new double[xCoords.length];
// 		ys[yIndex] = new double[xCoords.length];
// 		zs[yIndex] = new double[xCoords.length];
// 		foreach (xIndex, x; xCoords) {
// 			xs[yIndex][xIndex] = x;
// 			ys[yIndex][xIndex] = y;
// 			zs[yIndex][xIndex] = isValidBackbounce(x, y) ? 0.0f : double.nan;
// 		}
// 	}

// 	context.xs = d_to_python_numpy_ndarray(xs);
// 	context.ys = d_to_python_numpy_ndarray(ys);
// 	context.zs = d_to_python_numpy_ndarray(zs);

// 	context.py_stmts(script);
// }
