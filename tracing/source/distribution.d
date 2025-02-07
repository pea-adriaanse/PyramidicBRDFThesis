module distribution;

import std.conv : to;
import std.random : uniform01;
import std.stdio : File;

interface Distribution {
	float sample();
}

class ConstantDistribution : Distribution {
	float height;
	this(float height) {
		this.height = height;
	}

	float sample() {
		return this.height;
	}
}

/// Models distribution using histogram
class HistogramDistribution : Distribution {
	float[] cdf = [];
	float rangeSize; // relates #buckets to range.
	float bucketSize;

	this(string file, float rangeSize) {
		this.rangeSize = rangeSize;
		foreach (char[] l; File(file).byLine()) {
			this.cdf ~= l.to!float;
		}
		this.bucketSize = rangeSize / cdf.length;
	}

	float sample() {
		float prob = uniform01();
		foreach_reverse (bucket, cumulativeProb; cdf) {
			if (cumulativeProb <= prob) {
				return (bucket + 1) * bucketSize;
			}
		}
		assert(0);
	}
}
