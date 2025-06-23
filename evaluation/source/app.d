import std : iota;
import std.algorithm : clamp, min, max;
import std.conv : to;
import std.datetime;
import std.exception : enforce;
import std.file;
import std.math;
import std.random;
import std.stdio;

import landscape;
import pyramid_shape;
import ray;
import tracing;
import vdmath;
import vdmath.misc;
import backscattering;
import plotting;
import distribution;

import pyd.pyd;
import pyd.embedded;
import pyd.extra;

immutable uint pathLength = 4;
immutable uint meshPathLength = 4;
// immutable uint meshSampleCount = 4096;
immutable uint meshSampleCount = 2048;
immutable bool useShadowingPaul = true;
immutable uint resolution = 80;
// immutable uint resolution = 64;
immutable bool unfrequentProgressBar = true; // slightly faster
immutable bool backScatteringCorrection = false;

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

struct EvaluationResults {
	float[resolution][resolution] errors; // abs? average?
	float[resolution][resolution] exitProbBRDF;
	float[resolution][resolution] exitProbMesh;
	float[resolution] thetas;
	float[resolution] phis;
	double averageError;
	double maxError;

	void save(string dir) {
		mkdirRecurse(dir);
		writeln("Saving results to directory: ", dir);
		File(dir ~ "thetas.txt", "w").writeln(thetas);
		File(dir ~ "phis.txt", "w").writeln(phis);
		File(dir ~ "errors.txt", "w").writeln(errors);
		File(dir ~ "exitProbBRDF.txt", "w").writeln(exitProbBRDF);
		File(dir ~ "exitProbMesh.txt", "w").writeln(exitProbMesh);
		File(dir ~ "averageError.txt", "w").writeln(averageError);
		File(dir ~ "maxError.txt", "w").writeln(maxError);
		File(dir ~ "settings.txt", "w").writeln(
			"pathLength=", pathLength, '\n',
			"meshPathLength=", meshPathLength, '\n',
			"meshSampleCount=", meshSampleCount, '\n',
			"useShadowingPaul=", useShadowingPaul, '\n',
			"resolution=", resolution, '\n',
			"backScatteringCorrection=", backScatteringCorrection
		);
	}
}

void main(string[] args) {
	Landscape landscape = createLandscape();
	landscape.save("landscape.ply");
	// Landscape variableLandscape = createVariableLandscape();
	// variableLandscape.save("variableLandscape.ply");

	// plotShadowingDifference(landscape);

	// float theta = degreesToRadians(85.0);
	// float phi = 0;
	// Vec!3 wo = angleToVector(theta, phi);
	// Histogram!3 brdfHist = calculateBRDF!(3, true)(wo, theta, phi, &landscape.shape);
	// Histogram!3 meshHist = calculateMesh!3(wo, landscape);
	// brdfHist.writeln;
	// meshHist.writeln;

	// plotMeshExitProb(landscape);

	EvaluationResults* results = runEvaluation(landscape);
	// EvaluationResults* results = runEvaluation(variableLandscape);
	string dir = "./results/"; // default
	if (args.length == 2)
		dir = args[1];
	results.save(dir);

	plotEvaluation2(dir);
}

float calculateMeshExitProbAngles(uint pathLength)(float theta, float phi, const Landscape* landscape) {
	Vec!3 wo = angleToVector(theta, phi);
	return calculateMeshExitProb!pathLength(wo, landscape);
}

void plotMeshExitProb(const Landscape landscape) {
	Bounds!float xBounds = Bounds!float(0, PI_2 - PI / 16, PI / 16);
	Bounds!float yBounds = Bounds!float(0, PI - PI / 16, PI / 16);

	plotFunction!(float, float, const Landscape*)(xBounds, yBounds, ZLimits(0, 1, true), &calculateMeshExitProbAngles!5, "", "results/meshExitProb.png", &landscape);
}

void plotShadowingDifference(const Landscape landscape) {
	Bounds!float xBounds = Bounds!float(0, PI_2, PI / 64);
	Bounds!float yBounds = Bounds!float(0, PI, PI / 64);
	plotFunction!(float, float, const PyramidShape*)(xBounds, yBounds, ZLimits(0, 1, false), &G1Angles!true, "", "results/shadowing/shadowing_function_paul.png", &landscape
			.shape);
	plotFunction!(float, float, const PyramidShape*)(xBounds, yBounds, ZLimits(0, 1, false), &G1Angles!false, "", "results/shadowing/shadowing_function_lyanne.png", &landscape
			.shape);
	plotFunction!(float, float, const Landscape*)(xBounds, yBounds, ZLimits(0, 1, false), &calculateMeshShadowingAngle, "", "results/shadowing/shadowing_function_mesh.png", &landscape);
	// TODO: SCHRIJF BIJ MESH OP HOE VEEL SAMPLES & ZO IN DE NAAM

	float function(float, float, const PyramidShape*) diff = (a, b, c) => (
		G1Angles!true(a, b, c) - G1Angles!false(a, b, c));

	plotFunction!(float, float, const PyramidShape*)(xBounds, yBounds, ZLimits(0, 1, true), diff, "", "results/shadowing/shadowing_function_paul_minus_lyanne.png", &landscape
			.shape);
}

EvaluationResults* runEvaluation(const Landscape landscape) {
	SysTime startTime = Clock.currTime();
	immutable uint totalSamples = resolution * (resolution + 1);

	// float[theta][phi] -> access with arr[phi][theta].
	float[resolution][resolution] errors; // abs? average?
	float[resolution][resolution] exitProbBRDF;
	float[resolution][resolution] exitProbMesh;
	float[resolution] thetas;
	float[resolution] phis;
	double averageError = 0;
	double maxError = 0;

	immutable float PI_resolution = PI / resolution;
	immutable float PI_2_resolution = PI_2 / resolution;
	uint sample_i = 0;
	foreach (phi_i; 0 .. resolution) {
		float phi = phi_i * PI_resolution;
		phis[phi_i] = phi;
		foreach (theta_i; 0 .. resolution) { // theta = PI_2 will cause problems (eg. never hit surface)
			float theta = clamp(theta_i * PI_2_resolution, 0, PI_2);
			thetas[theta_i] = theta;

			Vec!3 wo = angleToVector(theta, phi);
			if (wo.z < 0)
				wo.z = 0; // turns out cos(PI_2) returns -2.71051e-20L ...
			wo = wo.normalize(); // paranoia

			Histogram!pathLength brdfHistogram = calculateBRDF!(pathLength, backScatteringCorrection)(
				wo, theta, phi, &landscape.shape);
			Histogram!meshPathLength meshHistogram;

			try {
				meshHistogram = calculateMesh!meshPathLength(wo, landscape);
			} catch (Error e) {
				writeln(theta_i);
				writeln(theta);
				writeln(wo);
				throw e;
			}

			exitProbBRDF[phi_i][theta_i] = brdfHistogram.exitProb;
			exitProbMesh[phi_i][theta_i] = meshHistogram.exitProb;
			float error = absErrorSum(brdfHistogram, meshHistogram);
			errors[phi_i][theta_i] = error;
			averageError += error;
			maxError = max(maxError, error);

			debug if (error > 1) {
				writeln("theta: ", theta);
				writeln("phi: ", phi);
				writeln("wo: ", wo.toString());
				writeln("brdf: ", brdfHistogram.exitPathProb);
				writeln("mesh: ", meshHistogram.exitPathProb);
			}

			sample_i += 1;
			static if (!unfrequentProgressBar)
				write("\rProgress: ", sample_i, "/", totalSamples, "\t\t\t\t\r");
		}
		static if (unfrequentProgressBar)
			write("\rProgress: ", sample_i, "/", totalSamples, "\t\t\t\t\r");
	}
	averageError /= resolution ^^ 2;
	writeln('\a');

	Duration duration = Clock.currTime() - startTime;
	writeln("Finished Evaluation after Duration: ", duration, "\nAverage Error: ", averageError, "\nMax Error: ", maxError);

	EvaluationResults* results = new EvaluationResults;
	results.errors = errors;
	results.exitProbBRDF = exitProbBRDF;
	results.exitProbMesh = exitProbMesh;
	results.thetas = thetas;
	results.phis = phis;
	results.averageError = averageError;
	results.maxError = maxError;
	return results;
}

void plotEvaluation(EvaluationResults* results) {
	static immutable python_script = import("evaluation/plot.py");
	if (!py_init_called) {
		py_init();
	}
	InterpContext context = new InterpContext();
	context.thetas_x = d_to_python_numpy_ndarray(results.thetas);
	context.phis_y = d_to_python_numpy_ndarray(results.phis);
	context.errors_z = d_to_python_numpy_ndarray(results.errors);
	context.exitProbBRDF_z = d_to_python_numpy_ndarray(results.exitProbBRDF);
	context.exitProbMesh_z = d_to_python_numpy_ndarray(results.exitProbMesh);
	context.py_stmts(python_script);
}

void plotEvaluation2(string dir) {
	import std.process : execute;

	// Doesnt seem to work:
	auto res = execute(["py", "./views/evaluation/plot2.py", dir]);
	writeln(res.output);
}

float absErrorSum(uint pathLength, uint meshPathLength)(
	Histogram!pathLength calculated, Histogram!meshPathLength reference) {
	// Compile time logic for path length mismatch:
	enum bool sameLengths = (pathLength == meshPathLength);
	static if (meshPathLength > pathLength)
		alias longerHist = reference;
	else static if (meshPathLength < pathLength)
		alias longerHist = calculated;
	enum uint maxPathLength = max(pathLength, meshPathLength);
	enum uint minPathLength = min(pathLength, meshPathLength);

	// Standard error sum:
	float errorSum = 0;
	foreach (i; 0 .. reflectIDCount(minPathLength)) {
		errorSum += abs(calculated.exitPathProb[i] - reference.exitPathProb[i]);
	}
	// Mismatch addition:
	static if (!sameLengths) {
		foreach (i; reflectIDCount(minPathLength) .. reflectIDCount(maxPathLength))
			errorSum += longerHist.exitPathProb[i];
	}

	return errorSum;
}

auto calculateBRDF(uint pathLength, bool backScatteringCorrection)(Vec!3 wo, float theta, float phi, const PyramidShape* pyramidShape)
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
			stderr.writeln("Dropped invalid throughput from parent (likely due to floating point error):\n\tgroupID: ",
				groupID, "\n\tgroupPath: ", reflectIDToString(groupID), "\n\tgroupWo:", groupWo.toString(), "\n\tthroughput: ", parentProb); // problem with shadowing function! (or floating point error)
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
			float occurProb = prob - exitProb; // = prob * (1-shadowing)

			exitProbs[pathID] = exitProb;
			occurProbs[pathID] = occurProb;
		}
	}

	if (backScatteringCorrection) {
		enforce(pathLength == 3, "Not implemented for pathLengths besides 3");
		string[4] dirs2 = ["EW", "NS", "WE", "SN"];
		string[4] dirs3 = ["EWE", "NSN", "WEW", "SNS"];
		uint[4] dirs2IDs;
		uint[4] dirs3IDs;
		double[4] deltaPhi = [0, -PI_2, PI, PI_2];
		foreach (i; 0 .. 4) {
			dirs2IDs[i] = reflectStringToID(dirs2[i]);
			dirs3IDs[i] = reflectStringToID(dirs3[i]);

			double backScatteringProb = backBounce(theta, phi + deltaPhi[i]);
			double backBounceProb = backScatteringProb * exitProbs[dirs2IDs[i]];
			exitProbs[dirs2IDs[i]] -= cast(float) backBounceProb;
			exitProbs[dirs3IDs[i]] += cast(float) backBounceProb;
		}
	}

	Histogram!pathLength histogram = Histogram!pathLength(exitProbs); // This is where ownership support could prevent copying (when allocated on the heap*).
	return histogram;
}

float G1Angles(bool shadowingPaul)(float theta, float phi, const PyramidShape* pyramidShape) {
	Vec!3 wo = angleToVector(theta, phi);
	return G1!shadowingPaul(wo, pyramidShape);
}

float G1(bool shadowingPaul)(Vec!3 wo, const PyramidShape* pyramidShape) {
	static if (shadowingPaul)
		return shadowing_paul(wo, pyramidShape);
	else
		return shadowing_lyanne(wo, pyramidShape.normalE);
}

float shadowing_paul(Vec!3 wo, const PyramidShape* pyramidShape) {
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

float calculateMeshShadowingAngle(float theta, float phi, const Landscape* landscape) {
	Vec!3 wo = angleToVector(theta, phi);
	return calculateMeshShadowing(wo, landscape);
}

float calculateMeshShadowing(Vec!3 wo, const Landscape* landscape) {
	Vec!3[] sampleTargets = Landscape.createPoints(landscape.width * 2 / 3.0, meshSampleCount, 0); // using non 0 target z will cause ray to exist below surface.
	uint escapedCount = 0;
	foreach (sample; sampleTargets) {
		Ray zRay = Ray(sample, Vec!3(0, 0, -1));
		Hit zHit = traceLand(*landscape, zRay);
		assert(zHit.hit, "Surface missing at " ~ sample.toString());
		sample.z = zHit.pos.z;
		uint hitPeak = zHit.peakID;

		Ray ray = Ray(sample, wo, hitPeak);
		Hit hit = traceLand(*landscape, ray);
		if (!hit.hit)
			escapedCount += 1;
	}
	return (cast(float) escapedCount) / meshSampleCount;
}

float calculateMeshExitProb(uint pathLength)(Vec!3 wo, const Landscape* landscape) {
	Vec!3[] sampleTargets = Landscape.createPoints(landscape.width * 2 / 3.0, meshSampleCount, 0); // using non 0 target z will cause ray to exist below surface.
	uint hits = 0;
	foreach (i; 0 .. meshSampleCount) {
		Vec!3 target = sampleTargets[i];
		Vec!3 origin = target + wo;
		Ray ray = Ray(origin, -wo);

		ReflectData!false reflectData = reflectRecurse!false(*landscape, ray, pathLength);
		if (reflectData.exits && reflectData.hitSurface()) { // Path Exited Surface
			hits += 1;
		}
	}
	float exitProb = (cast(float) hits) / meshSampleCount;
	return exitProb;
}

auto calculateMesh(uint pathLength)(Vec!3 wo, const ref Landscape landscape, float maxHeight = 0) {
	immutable uint pathCount = reflectIDCount(pathLength);
	uint[pathCount] exitPathBins;

	Vec!3[] sampleTargets = Landscape.createPoints(landscape.width * 2 / 3.0, meshSampleCount, 0); // using non 0 target z will cause ray to exist below surface.
	foreach (i; 0 .. meshSampleCount) {
		Vec!3 target = sampleTargets[i];

		enforce(wo.z > 0, "Can't sample with wo.z <= 0");
		float offsetScale = 1; // constant height surface
		if (maxHeight > 0)
			offsetScale = max(1.0, maxHeight / wo.z);

		Vec!3 origin = target + wo * offsetScale; // offsetScale ensures it's above the surface.
		Ray ray = Ray(origin, -wo);

		ReflectData!true reflectData = reflectRecurse!true(landscape, ray, pathLength);
		if (reflectData.exits && reflectData.hitSurface()) { // Path Exited Surface
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
	landscape.createBins(binCount, 0, -landscape.approxHeight, true);
	return landscape;
}

Landscape createVariableLandscape() {
	immutable density = 0.6;
	immutable float slope = degreesToRadians(54.7);
	immutable uint binCount = 25;
	immutable float width = 50;
	immutable string heightsFile = "../tracing/heightCumulative.csv";
	immutable float maxVarHeight = 7.0;
	Distribution heightDistribution = new HistogramDistribution(heightsFile, maxVarHeight);
	Landscape landscape = Landscape(width, density, slope, heightDistribution);
	// no binning: not tested for variable height
	landscape.createBins(binCount, 7, -landscape.approxHeight, true); // approx height very rough in this case (likely too large).
	return landscape;
}
