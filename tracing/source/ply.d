module ply;
import math;
import std.exception : enforce;
import std.file;
import std.path;
import std.stdio : File;

struct Ply {
	Vec!3[] vertices;
	uint[3][] faces;

	void saveToFile(string path) {
		enforce(isValidPath(path) && isValidFilename(baseName(path)), "Not a valid PLY output path: " ~ path);
		path = defaultExtension(path, ".ply");

		File file = File(path, "w");
		file.writeln("ply");
		file.writeln("format ascii 1.0");
		file.writeln("comment Created by Paul Adriaanse");

		file.writeln("element vertex ", vertices.length);
		file.writeln("property float x");
		file.writeln("property float y");
		file.writeln("property float z");

		file.writeln("element face ", faces.length);
		file.writeln("property list uint8 uint32 vertex_indices");

		file.writeln("end_header");

		foreach (v; vertices)
			file.writeln(v.x, " ", v.y, " ", v.z);

		foreach (f; faces)
			file.writeln("3 ", f[0], " ", f[1], " ", f[2]);
	}
}
