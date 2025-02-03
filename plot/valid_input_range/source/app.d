import pyd.pyd;
import pyd.embedded;
import pyd.extra;
import std.math.constants : PI, PI_2;
import std.range : iota;
import vdmath;
import vdmath.misc : degreesToRadians;
import std.math;
import std.conv;
import std.array : array;
import std.algorithm : min, max;
import std.stdio;

immutable script = import("plot.py");

shared static this() {
	py_init();
}

enum double alpha = degreesToRadians(54.7);
enum double rho = 0.6;
Vec!(3, double) eastNormal = Vec!(3, double)(sin(alpha), 0, cos(alpha));
Vec!(3, double) westNormal = Vec!(3, double)(-sin(alpha), 0, cos(alpha));

Vec!(3, double) reflect(Vec!(3, double) inDir, Vec!(3, double) normal) {
	return inDir - ((normal * inDir.dot(normal)) * 2);
}

Vec!(3, double) getInDir(double theta, double phi) {
	return -Vec!(3, double)(
		sin(theta) * cos(phi),
		sin(theta) * sin(phi),
		cos(theta)
	);
}

/// Valid if can hit ->E->W, can subsequently escape (z>0) & can hit ->E.
bool isValid(Vec!(3, double) inDir) {
	if (inDir.dot(eastNormal) >= 0)
		return false;
	Vec!(3, double) reflectDir = reflect(inDir, eastNormal);
	if (reflectDir.dot(westNormal) >= 0)
		return false;
	Vec!(3, double) outDir = reflect(reflectDir, westNormal);
	if (outDir.z <= 0)
		return false;
	if (outDir.dot(eastNormal) >= 0)
		return false;
	return true;
}

double theta_copy;
double phi_copy;
double D1Copy(double y, double z) {
	return ((-7.0671883 * z + 9.98134798 * y) * cos(phi_copy) * sin(theta_copy) + (
			0.43224778 * z - 0.61048541 * y) * cos(theta_copy)) / ((3.04505424 * cos(
			phi_copy) ^^ 2 + 1.81997208 * cos(phi_copy) * sin(phi_copy)) * sin(theta_copy) ^^ 2 + (
			-9.82233945 * cos(phi_copy) - 6.38906164 * sin(phi_copy)) * cos(
			theta_copy) * sin(
			theta_copy) - 3.04505424 * cos(theta_copy) ^^ 2);
}

double D1n(double theta, double phi, double y, double z) {
	return ((-7.0671883 * z + 9.98134798 * y) * cos(phi) * sin(theta) + (
			0.43224778 * z - 0.61048541 * y) * cos(theta)) / ((3.04505424 * cos(
			phi) ^^ 2 + 1.81997208 * cos(phi) * sin(phi)) * sin(theta) ^^ 2 + (
			-9.82233945 * cos(phi) - 6.38906164 * sin(phi)) * cos(theta) * sin(
			theta) - 3.04505424 * cos(theta) ^^ 2);
}

double D2n(double theta, double phi, double y, double z) {
	return ((-7.0671883 * z - 9.98134798 * y) * cos(phi) * sin(theta) + (
			0.43224778 * z + 0.61048541 * y) * cos(theta)) / ((3.04505424 * cos(
			phi) ^^ 2 - 1.81997208 * cos(phi) * sin(phi)) * sin(theta) ^^ 2 + (
			-9.82233945 * cos(phi) + 6.38906164 * sin(phi)) * cos(theta) * sin(
			theta) - 3.04505424 * cos(theta) ^^ 2);
}

private double internal_theta;
private double internal_phi;
private double internal_yMinFactor;
private double internal_yMaxFactor;
double internalExpression(alias Dfunc)(double y, double z) {
	if (y <= internal_yMinFactor * z || y >= internal_yMaxFactor * z)
		return 0; // Ensures bounds.
	double D = Dfunc(internal_theta, internal_phi, y, z);
	if (D <= 0)
		return 0; // Edge case, likely due to floating point error.

	// return Tnormal(internal_theta, internal_phi, z);
	return Tdni(D, internal_theta, internal_phi, z);
}

double TnormaliCropped(double y, double z) {
	if (y <= internal_yMinFactor * z || y >= internal_yMaxFactor * z)
		return 0; // Ensures bounds.
	return Tnormali(internal_theta, internal_phi, z);
}

double Tnormali(double theta, double phi, double z) {
	return 1.225283594 * exp(-1.203167728 * z ^^ 2);
}

double Tnormal() {
	return 0.7210541991;
}

// double T1ni(double phi, double theta, double y, double z) {

// }

// double T2ni(double phi, double theta, double y, double z) {

// }

double Tdni(double D, double theta, double phi, double z) {
	assert(D >= 0, D.to!string);
	return -1.225283594 * exp(-1.203167728 * z ^^ 2) *
		(
			exp(0.2822198290 * D * (sin(theta) * cos(phi) + 2.839653910 * cos(theta)) * (
				0.4716113289 * sin(theta) * cos(phi) * D - 0.1660805661 * cos(theta) * D + z))
				- 1);
}

double Yminn_z() {
	return 0.708039467;
}

double Ymaxn_z() {
	return -0.708039467;
}

double Ymid_z(double phi) {
	// return 1 / (0.6731323870 + 2.363049614 * (1 / tan(phi)));
	return 1 / (0.6731323870 + 2.363049614 * (cos(phi) / sin(phi)));
}

double backBounceNormal(double theta, double phi) {
	Vec!(3, double) inDir = getInDir(theta, phi);
	if (!isValid(inDir))
		return double.nan;

	IntegralBounds zBounds = IntegralBounds(-4, 0, 100);
	IntegralBounds yBounds = IntegralBounds(Yminn_z() * -4, Ymaxn_z() * -4, 512);
	internal_theta = theta;
	internal_phi = phi;
	internal_yMinFactor = Yminn_z();
	internal_yMaxFactor = Ymaxn_z();
	double normalFactor = integrate(&TnormaliCropped, yBounds, zBounds);
	return normalFactor;
}

double backBounce(double theta, double phi) {
	phi = -phi; // Flip sign to sidestep numeric discontinuity.

	Vec!(3, double) inDir = getInDir(theta, phi);
	// Test for invalid preconditions
	if (!isValid(inDir))
		return 0.0; // double.nan also viable but problematic with interpolation;

	// Test for guaranteed backbounce
	double d1_test = D1n(theta, phi, 0, -1);
	double d2_test = D2n(theta, phi, 0, -1);
	if (d1_test <= 0 && d2_test <= 0)
		return 1.0;

	double YminFactor = Yminn_z();
	double YmaxFactor = Ymaxn_z();
	double YmidFactor = Ymid_z(phi);
	// Note z is negative.
	bool integral1Zero = false;
	bool integral2Zero = false;
	if (YmidFactor >= YminFactor) {
		YmidFactor = YminFactor;
		integral1Zero = true;
	} else if (YmidFactor <= YmaxFactor) {
		YmidFactor = YmaxFactor;
		integral2Zero = true;
	}

	IntegralBounds zBounds = IntegralBounds(-4, 0, 100);
	IntegralBounds yBounds1 = IntegralBounds(YminFactor * -4, YmidFactor * -4, 1024); // Ensure bounds by making expression evaluate to 0 outside correct Y bounds.
	IntegralBounds yBounds2 = IntegralBounds(YmidFactor * -4, YmaxFactor * -4, 1024); // ^

	internal_theta = theta;
	internal_phi = phi;
	internal_yMinFactor = YminFactor;
	internal_yMaxFactor = YmidFactor;
	auto internalExpression1 = &internalExpression!(D1n);
	double integral1 = integral1Zero ? 0 : integrate(internalExpression1, yBounds1, zBounds);

	internal_yMinFactor = YmidFactor;
	internal_yMaxFactor = YmaxFactor;
	auto internalExpression2 = &internalExpression!(D2n);
	double integral2 = integral2Zero ? 0 : integrate(internalExpression2, yBounds2, zBounds);

	double finalIntegral = (integral1 + integral2) / Tnormal();
	return finalIntegral;
}

struct IntegralBounds {
	double min;
	double max;
	uint stepCount; // use stepCount to improve precision.
}

string _expressionWrapped(size_t count) {
	string ret = "expression(";
	foreach (i; 0 .. count - 1)
		ret ~= "coordinates[" ~ i.to!string ~ "],";
	return ret ~ "coordinates[$-1])";
}

double integrate(Arg...)(double function(Arg) expression, IntegralBounds[] bounds...) {
	assert(Arg.length == bounds.length);
	assert(Arg.length >= 1);
	debug foreach (bound; bounds)
		assert(bound.stepCount > 0);

	double integralSum = 0;
	double sampleWeight = 1;
	double[Arg.length] stepSizes;
	foreach (i; 0 .. Arg.length) {
		double diff = bounds[i].max - bounds[i].min;
		double stepSize = diff / bounds[i].stepCount;
		stepSizes[i] = stepSize;
		sampleWeight *= stepSize;
	}

	uint[Arg.length] sampleIDs; // 0 -> stepCount - 1
	double[Arg.length] coordinates;
	while (true) {
		foreach (i; 0 .. Arg.length)
			coordinates[i] = bounds[i].min + stepSizes[i] * (sampleIDs[i] + 0.5);

		double value = mixin(_expressionWrapped(Arg.length));
		integralSum += value * sampleWeight;

		foreach (i; 0 .. Arg.length) {
			uint nextID = sampleIDs[i] + 1;
			if (nextID == bounds[i].stepCount) {
				if (i + 1 == Arg.length)
					return integralSum;
				sampleIDs[i] = 0;
				continue; // increment next
			}
			sampleIDs[i] = nextID;
			break;
		}
	}
}

bool isValidAngle(double theta, double phi) {
	Vec!(3, double) inDir = getInDir(theta, phi);
	return isValid(inDir);
}

struct Bounds(T) {
	T min;
	T max; // inclusive!
	T stepSize;
}

void main(string[] args) {
	Bounds!double thetaBounds = Bounds!double(0, PI_2, PI / 256);
	Bounds!double phiBounds = Bounds!double(0, PI, PI / 256);
	// Bounds!double thetaBounds = Bounds!double(0, PI_2, PI / 128);
	// Bounds!double phiBounds = Bounds!double(0, PI, PI / 128);
	InterpContext context = new InterpContext();
	// plotFunction(context, thetaBounds, phiBounds, &isValidAngle, null);
	plotFunction(context, thetaBounds, phiBounds, &backBounce, "backBounce.bin");
	// plotFunction(context, thetaBounds, phiBounds, &backBounce, null);
	// plotFunction(context, thetaBounds, phiBounds, &backBounceNormal, null);
}

void plotFunction(T2, T, Args...)(InterpContext context, Bounds!T xBounds, Bounds!T yBounds, T2 function(T, T, Args) func, string savePath, Args args) {
	context.xBounds = *(cast(T[3]*)&xBounds);
	context.yBounds = *(cast(T[3]*)&yBounds);

	T[] xCoords = iota(xBounds.min, xBounds.max + xBounds.stepSize, xBounds.stepSize).array();
	T[] yCoords = iota(yBounds.min, yBounds.max + yBounds.stepSize, yBounds.stepSize).array();

	T[][] xs = new T[][yCoords.length];
	T[][] ys = new T[][yCoords.length];
	T2[][] zs = new T2[][yCoords.length];

	foreach (yIndex, y; yCoords) {
		xs[yIndex] = new T[xCoords.length];
		ys[yIndex] = new T[xCoords.length];
		zs[yIndex] = new T2[xCoords.length];
		foreach (xIndex, x; xCoords) {
			xs[yIndex][xIndex] = x;
			ys[yIndex][xIndex] = y;
			zs[yIndex][xIndex] = func(x, y, args);
		}
	}

	if (savePath !is null && savePath.length > 0)
		saveTable(xBounds, yBounds, zs, savePath);

	context.xs = d_to_python_numpy_ndarray(xs);
	context.ys = d_to_python_numpy_ndarray(ys);
	context.zs = d_to_python_numpy_ndarray(zs);

	context.py_stmts(script);
}

void saveTable(T, T2)(Bounds!T xBounds, Bounds!T yBounds, T2[][] data, string path) {
	File file = File(path, "w");
	ubyte[] xBoundsData = cast(ubyte[])(*(cast(T[3]*)&xBounds));
	ubyte[] yBoundsData = cast(ubyte[])(*(cast(T[3]*)&yBounds));
	string header = T.stringof ~ ' ' ~ T2.stringof ~ '\n'; // Data Types
	header ~= xBoundsData.length.to!string ~ ' ' ~ yBoundsData.length.to!string ~ '\n'; // Bounds data lengths
	header ~= data.length.to!string ~ ' ' ~ data[0].length.to!string ~ '\n'; // Dimensions
	file.write(header);

	file.rawWrite(xBoundsData);
	file.rawWrite(yBoundsData);

	ubyte[] rawData;
	foreach (T2[] segment; data)
		rawData ~= cast(ubyte[]) segment;
	file.rawWrite(rawData);
	file.close();

	// writeln(data);

	writeln(i"Data Byte Size:$(rawData.length)B");
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
