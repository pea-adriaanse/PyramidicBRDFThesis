module integration;

import std.conv : to;

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
