import std.stdio;
import reflectdist;
import std.file;
import std.conv;
import std.parallelism : parallel;
import std.range : iota;
import std.process;
import std.exception : enforce;

void main(string[] args) {
	if (args.length == 1)
		args ~= ["8", "35.7", "0", "4000", "2", "false", "test"];
	assert(args.length == 8, "Need 7 arguments: (polarAngleSplitCount polarAngleMax azimuthalAngle sampleCount reflectCount name parallel)");
	uint polarAngleSplitCount = args[1].to!uint;
	float polarAngleMax = args[2].to!float;
	assert(polarAngleMax > -90 && polarAngleMax < 90, "Need polarAngleMax between -90 and 90 (exclusively), but got: " ~ args[2]);

	string name = args[$ - 2];
	bool runParallel = args[$ - 1].to!bool;
	args = args[0] ~ args[2 .. $ - 1];

	auto res = executeShell("cd ../iridescence & cmake --build --preset cpuRelease");
	enforce(res.status == 0, "Iridescence failed to build:\n" ~ res.output);
	res = executeShell("cd ../tracing & dub build --build=release");
	enforce(res.status == 0, "Tracing failed to build:\n" ~ res.output);
	mkdirRecurse("./temp/");

	enforce(polarAngleSplitCount >= 1, "Zero polar splits will do nothing.");
	if (polarAngleSplitCount == 1) {
		args[1] = polarAngleMax.to!string;
		reflectDist(args);
		return;
	}

	void foreachBody(int i) {
		string[] parallelArgs = args.dup;
		// foreach (i; 0 .. polarAngleSplitCount) {
		float polarAngle = i * polarAngleMax / (polarAngleSplitCount - 1);
		string polarAngleStr = polarAngle.to!string;
		parallelArgs[1] = polarAngleStr; // Hijack arguments
		reflectDist(parallelArgs);
	}

	if (runParallel)
		foreach (i; parallel(iota(polarAngleSplitCount)))
			foreachBody(i);
	else
		foreach (i; 0 .. polarAngleSplitCount)
			foreachBody(i);

}
