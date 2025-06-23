module backscattering;

import vdmath;
import vdmath.misc : degreesToRadians;
import integration;

import std.conv;
import std.math;
import std.math.constants : PI, PI_2;
import std.range : iota;

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

// Incorrect: does not add -dz/(df tan(alpha)) factor to local integral
// double Tdni(double D, double theta, double phi, double z) {
// 	assert(D >= 0, D.to!string);
// 	return -1.225283594 * exp(-1.203167728 * z ^^ 2) *
// 		(
// 			exp(0.2822198290 * D * (sin(theta) * cos(phi) + 2.839653910 * cos(theta)) * (
// 				0.4716113289 * sin(theta) * cos(phi) * D - 0.1660805661 * cos(theta) * D + z))
// 				- 1);
// }

double Tdni(double D, double theta, double phi, double z) {
	assert(D >= 0, D.to!string);
	return -1.225283594 * (-1 + exp(
			-0.2852077020 * (sin(theta) * cos(phi) - 3.510527281 * cos(theta)) * D * (
			0.4716113289 * D * sin(theta) * cos(phi) - 0.1660805661 * D * cos(theta) + z))) * exp(
		-1.203167728 * z ^^ 2);
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

bool isValidAngle(double theta, double phi) {
	Vec!(3, double) inDir = getInDir(theta, phi);
	return isValid(inDir);
}
