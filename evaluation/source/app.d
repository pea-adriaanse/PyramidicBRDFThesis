import std : iota;
import std.conv : to;
import std.math;
import std.algorithm : max, clamp;
import std.stdio;

import landscape;
import pyramid_shape;
import ray;
import tracing;
import vdmath;
import vdmath.misc;
import pyd.pyd;
import pyd.embedded;
import pyd.extra;

immutable python_script = import("plot.py");

immutable uint pathLength = 3;
immutable uint meshSampleCount = 1024;
immutable bool useShadowingPaul = true;
immutable uint resolution = 64;

struct Histogram(uint pathLength = pathLength) {
	static immutable uint histogramEntryCount = reflectIDCount(pathLength);

	this(ref float[histogramEntryCount] exitProbs) {
		float totalProb = 0;
		this.exitPathProb = exitProbs;
		foreach (exitProb; exitProbs)
			totalProb += exitProb;
		this.exitProb = totalProb;
	}

	this(uint[histogramEntryCount] exitPathBins, uint meshSampleCount) {
		uint exitCount = 0;
		foreach (i; 0 .. histogramEntryCount) {
			exitCount += exitPathBins[i];
			this.exitPathProb[i] = (cast(float) exitPathBins[i]) / meshSampleCount;
		}
		assert(exitCount <= meshSampleCount);
		this.exitProb = (cast(float) exitCount) / meshSampleCount;
	}

	float[histogramEntryCount] exitPathProb;
	float exitProb;
}

void main() {
	Landscape landscape = createLandscape();

	immutable uint totalSamples = resolution * (resolution + 1);

	float[resolution + 1][resolution] error; // abs? average?
	float[resolution + 1][resolution] exitProbBRDF;
	float[resolution + 1][resolution] exitProbMesh;
	float[resolution + 1] thetas;
	float[resolution] phis;

	immutable float PI_resolution = PI / resolution;
	immutable float PI_2_resolution = PI_2 / resolution;
	uint sample_i = 0;
	foreach (phi_i; 0 .. resolution) {
		float phi = phi_i * PI_resolution;
		phis[phi_i] = phi;
		foreach (theta_i; 0 .. resolution + 1) {
			float theta = clamp(theta_i * PI_2_resolution, 0, PI_2);
			thetas[theta_i] = theta;
			Vec!3 wo;
			wo.x = sin(theta) * cos(phi);
			wo.y = sin(theta) * sin(phi);
			wo.z = max(cos(theta), 0); // turns out cos(PI_2) returns -2.71051e-20L ...
			wo = wo.normalize(); // paranoia?

			Histogram!pathLength brdfHistogram = calculateBRDF!pathLength(wo, landscape.shape);
			Histogram!pathLength meshHistogram = calculateMesh!pathLength(wo, landscape);
			exitProbBRDF[phi_i][theta_i] = brdfHistogram.exitProb;
			exitProbMesh[phi_i][theta_i] = meshHistogram.exitProb;
			error[phi_i][theta_i] = absErrorSum(brdfHistogram, meshHistogram);
			sample_i += 1;
			write("\rProgress: ", sample_i, "/", totalSamples, "\t\t\t\t\r");
			
			// Do Backscattering Correction
		}
	}
	writeln();

	writeln("Error:\n",error);


	// pyd plot
	py_init();
	InterpContext context = new InterpContext();
	context.thetas_x = d_to_python_numpy_ndarray(thetas);
	context.phis_y = d_to_python_numpy_ndarray(phis);
	context.errors_z = d_to_python_numpy_ndarray(error);
	context.exitProbBRDF_z = d_to_python_numpy_ndarray(exitProbBRDF);
	context.exitProbMesh_z = d_to_python_numpy_ndarray(exitProbMesh);
	context.py_stmts(python_script);
}

float absErrorSum(uint pathLength)(Histogram!pathLength calculated, Histogram!pathLength reference) {
	float errorSum = 0;
	foreach (i; 0 .. reflectIDCount(pathLength)) {
		errorSum += abs(calculated.exitPathProb[i] - reference.exitPathProb[i]);
	}
	return errorSum;
}

auto calculateBRDF(uint pathLength)(Vec!3 wo, PyramidShape pyramidShape)
		if (pathLength >= 1) {
	immutable uint pathCount = reflectIDCount(pathLength);
	float[pathCount] occurProbs;
	Vec!3[pathCount] reflectDirs;
	float[pathCount] exitProbs;

	for (uint groupID = 0; groupID < pathCount; groupID += 4) {
		uint parentID = 0;
		uint parentFace = 0;
		float parentProb = 1;
		Vec!3 groupWo = wo;

		if (groupID >= 4) {
			parentID = (groupID / 4) - 1;
			parentFace = parentID % 4;
			parentProb = occurProbs[parentID];
			groupWo = -reflectDirs[parentID];

			if (parentProb == 0) {
				non_occurring_group: foreach (face; 0 .. 4) {
					uint pathID = groupID + face;
					occurProbs[pathID] = 0;
					exitProbs[pathID] = 0;
				}
				continue;
			}
		}

		float[4] visibleNormalDistribution;
		float cosSum = 0;
		foreach (face; 0 .. 4) {
			float cos = max(0, pyramidShape.normals[face].dot(groupWo));
			visibleNormalDistribution[face] = cos;
			cosSum += cos;
		}
		if (cosSum <= 0) {
			string groupIDStr = groupID.to!string();
			string groupPathStr = reflectIDToString(groupID);
			string groupWoStr = groupWo.toString();
			string droppedThroughputStr = parentProb.to!string();
			stderr.writeln("Dropped invalid throughput from parent (likely due to floating point error):\n\tgroupID: " ~ groupIDStr
					~ "\n\tgroupPath: " ~ groupPathStr
					~ "\n\tgroupWo:" ~ groupWoStr
					~ "\n\tthroughput: " ~ droppedThroughputStr); // problem with shadowing function! (or floating point error)
			goto non_occurring_group; // There is no clear solution. Options include offsetting throughput to parent or dropping it (as done here).
		}
		foreach (face; 0 .. 4) {
			visibleNormalDistribution[face] /= cosSum;
		}

		foreach (face; 0 .. 4) {
			uint pathID = groupID + face;
			Vec!3 outDir = reflect(-groupWo, pyramidShape.normals[face]).normalize();
			reflectDirs[pathID] = outDir;

			float prob = parentProb * visibleNormalDistribution[face];
			float shadowing = G1!useShadowingPaul(outDir, pyramidShape);
			float exitProb = prob * shadowing;
			float occurProb = prob - exitProb;

			exitProbs[pathID] = exitProb;
			occurProbs[pathID] = occurProb;
		}
	}

	Histogram!pathLength histogram = Histogram!pathLength(exitProbs); // This is where ownership support could prevent copying (when allocated on the heap*).
	return histogram;
}

float G1(bool shadowingPaul)(Vec!3 wo, ref PyramidShape pyramidShape) {
	static if (shadowingPaul)
		return shadowing_paul(wo, pyramidShape);
	else
		return shadowing_lyanne(wo, pyramidShape.normalE);
}

float shadowing_paul(Vec!3 wo, ref PyramidShape pyramidShape) {
	float dotSum = 0;
	foreach (face; 0 .. 4) {
		float dot = max(0, wo.dot(pyramidShape.normals[face]));
		dotSum += dot;
	}
	if (dotSum == 0)
		return 0;
	float cosTheta = wo.z;
	float cosAlpha = pyramidShape.normalE.z; // cos(pyramidShape.slope)
	// return 4 * cosAlpha * cosTheta / dotSum;
	return clamp(4 * cosAlpha * cosTheta / dotSum, 0, 1); // note clamp is required due to floatingpoint error
}

float shadowing_lyanne(Vec!3 wo, Vec!3 normalEast) {
	Vec!3 zeroAzimuthWo = zeroAzimuth(wo);
	if (zeroAzimuthWo.dot(normalEast) < 0)
		return 0;

	// c++ version possibly wrong!!!
	// possibly does not ensure the 1 value
	// Correction: clamping to 1 fixes this in the standard domain.! :3

	float cosTheta = wo.z;
	float cosAlpha = normalEast.z;
	float cosThetaMinAlpha = zeroAzimuthWo.dot(normalEast); // correct when clamped
	return clamp(4 * cosTheta * cosAlpha / (2 * cosTheta * cosAlpha + cosThetaMinAlpha), 0, 1);
}

Vec!3 zeroAzimuth(Vec!3 vec) {
	float x = sqrt(1.0 - (vec.z ^^ 2));
	return Vec!3(x, 0, vec.z);
}

auto calculateMesh(uint pathLength)(Vec!3 wo, const ref Landscape landscape) {
	immutable uint pathCount = reflectIDCount(pathLength);
	uint[pathCount] exitPathBins;

	// not nesessary but interesting:
	float peakOccurranceZ = tan(degreesToRadians(54.7)) / sqrt(8 * 0.6); // roughly 0.6646

	Vec!3[] sampleTargets = Landscape.createPoints(landscape.width * 2 / 3.0, meshSampleCount, -peakOccurranceZ); // z!=0 as theta=pi/2 will miss surface
	foreach (i; 0 .. meshSampleCount) {
		Vec!3 target = sampleTargets[i];
		Vec!3 origin = target + wo;
		Ray ray = Ray(origin, -wo);

		ReflectData!true reflectData = reflectRecurse!true(landscape, ray, pathLength);
		if (reflectData.exits) { // Path Exited Surface
			exitPathBins[reflectData.reflectID] += 1;
		}
	}

	Histogram!pathLength histogram = Histogram!pathLength(exitPathBins, meshSampleCount);
	return histogram;
}

Landscape createLandscape() {
	immutable density = 0.6;
	immutable slope = degreesToRadians(54.7);
	immutable uint binCount = 25;
	immutable float width = 500;
	Landscape landscape = Landscape(width, density, slope);
	landscape.createBins(binCount);
	return landscape;
}
