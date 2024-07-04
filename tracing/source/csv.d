module csv;
import std.array : Appender;
import std.conv : to;
import std.file;
import std.path;
import std.exception;

struct CSV {
	auto text = Appender!string(null);
	char sep;

	@disable this();

	this(char sep, bool addSep, string[] header...) {
		this.sep = sep;
		if (addSep)
			text ~= "sep=" ~ sep ~ '\n';

		if (header.length > 0) {
			foreach (i; 0 .. header.length - 1) {
				text ~= header[i];
				text ~= sep;
			}
			text ~= header[$ - 1];
			text ~= '\n';
		}
	}

	void addEntry(T...)(T args) if (args.length > 0) {
		static foreach (i; 0 .. args.length - 1) {
			text ~= args[i].to!string;
			text ~= sep;
		}
		text ~= args[$ - 1].to!string;
		text ~= '\n';
	}

	void save(string path) {
		enforce(isValidPath(path) && isValidFilename(baseName(path)), "Not a valid output path: " ~ path);
		path = defaultExtension(path, ".csv");

		write(path, text.data);
	}
}
