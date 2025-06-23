module plotting;

import std.array : array;
import std.conv : to;
import std.stdio;
import std.range : iota;

import pyd.embedded;
import pyd.extra;
import pyd.pyd;

struct Bounds(T) {
	T min;
	T max; // inclusive!
	T stepSize;
}

struct ZLimits {
	float min;
	float max;
	bool automatic;
}

static immutable script = import("plotting/plot.py");
void plotFunction(T2, T, Args...)(Bounds!T xBounds, Bounds!T yBounds, ZLimits zLimits, T2 function(T, T, Args) func, string saveTablePath, string savePlotPath="temp.png", Args args) {
	if (!py_init_called) {
		py_init();
	}
	InterpContext context = new InterpContext();
	context.xBounds = *(cast(T[3]*)&xBounds);
	context.yBounds = *(cast(T[3]*)&yBounds);

	T[] xCoords = iota(xBounds.min, xBounds.max + xBounds.stepSize / 2, xBounds.stepSize).array();
	T[] yCoords = iota(yBounds.min, yBounds.max + yBounds.stepSize / 2, yBounds.stepSize).array();

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

	if (saveTablePath !is null && saveTablePath.length > 0)
		saveTable(xBounds, yBounds, zs, saveTablePath);

	context.xs = d_to_python_numpy_ndarray(xs);
	context.ys = d_to_python_numpy_ndarray(ys);
	context.zs = d_to_python_numpy_ndarray(zs);
	context.zLimits = [zLimits.min, zLimits.max];
	context.zLimitsAutomatic = zLimits.automatic;
	context.savePlotPath = savePlotPath;

	context.py_stmts(script);
}

void saveTable(T, T2)(Bounds!T xBounds, Bounds!T yBounds, T2[][] data, string path) {
	File file = File(path, "w");
	ubyte[] xBoundsData = cast(ubyte[])(*(cast(T[3]*)&xBounds));
	ubyte[] yBoundsData = cast(ubyte[])(*(cast(T[3]*)&yBounds));
	string header = T.stringof ~ ' ' ~ T2.stringof ~ '\n'; // Data Types
	header ~= xBoundsData.length.to!string ~ ' ' ~ yBoundsData.length.to!string ~ '\n'; // Bounds data lengths
	header ~= data[0].length.to!string ~ ' ' ~ data.length.to!string ~ '\n'; // Dimensions
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
