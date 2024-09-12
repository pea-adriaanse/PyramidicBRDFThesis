module ply;
import vdmath;
import std.exception : enforce;
import std.file;
import std.path;
import std.array : Appender;
import std.conv : to;
import std.format;
import std.parallelism : totalCPUs;
import core.thread;

struct Ply {
	immutable Vec!3[] vertices;
	immutable uint[3][] faces;

	static void saveToFile(immutable Vec!3[] vertices, immutable uint[3][] faces, string path) {
		saveToFile(vertices, faces, null, path);
	}

	static void saveToFile(immutable Vec!3[] vertices, immutable uint[3][] faces, immutable Vec!3[] normals, string path) {
		immutable bool useNormals = normals !is null;
		enforce(!useNormals || normals.length == vertices.length);
		enforce(isValidPath(path) && isValidFilename(baseName(path)), "Not a valid PLY output path: " ~ path);
		path = defaultExtension(path, ".ply");

		auto text = Appender!string(null);
		text ~= `ply
format ascii 1.0
comment Created by Paul Adriaanse
element vertex `;
		text ~= vertices.length.to!string;
		text ~= `
property float x
property float y
property float z
`;
		if (useNormals)
			text ~= `property float nx
property float ny
property float nz
element face `;
		text ~= faces.length.to!string;
		text ~= `
property list uint8 uint32 vertex_indices
end_header
`;

		__gshared string[] textsV;
		__gshared string[] textsF;
		textsV = new string[totalCPUs];
		textsF = new string[totalCPUs];

		class FormatThread : Thread {
			immutable Vec!3[] vertices;
			immutable uint[3][] faces;
			uint i;

			this(uint i, immutable Vec!3[] vertices, immutable uint[3][] faces) {
				this.vertices = vertices;
				this.faces = faces;
				this.i = i;
				super(&fun);
			}

			void fun() {
				Appender!string textV = Appender!string(null);
				Appender!string textF = Appender!string(null);

				if (useNormals)
					foreach (i; 0 .. vertices.length)
						formattedWrite(textV, "%.9f %.9f %.9f %.9f %.9f %.9f\n", vertices[i].x,
							vertices[i].y, vertices[i].z, normals[i].x, normals[i].y, normals[i].z);
				else
					foreach (v; this.vertices)
						formattedWrite(textV, "%.9f %.9f %.9f\n", v.x, v.y, v.z);

				foreach (f; this.faces)
					formattedWrite(textF, "3 %u %u %u\n", f[0], f[1], f[2]);

				textsV[i] = textV.data();
				textsF[i] = textF.data();
			}
		}

		Thread[] threads = new Thread[totalCPUs];
		size_t widthV = vertices.length / totalCPUs;
		size_t widthF = faces.length / totalCPUs;
		foreach (uint i; 0 .. totalCPUs) {
			size_t endV;
			size_t endF;
			if (i == totalCPUs - 1) {
				endV = vertices.length;
				endF = faces.length;
			} else {
				endV = (i + 1) * widthV;
				endF = (i + 1) * widthF;
			}
			threads[i] = new FormatThread(i, vertices[i * widthV .. endV], faces[i * widthF .. endF]).start();
		}
		foreach (i, Thread t; threads) {
			t.join();
			text ~= textsV[i];
		}
		foreach (i; 0 .. totalCPUs)
			text ~= textsF[i];

		write(path, text.data);
	}
}
