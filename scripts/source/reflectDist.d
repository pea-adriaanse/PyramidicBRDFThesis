module reflectdist;
import std.process;
import std.stdio;
import std.file;
import std.exception : enforce;
import std.conv : to;
import std.format : format;
import std.array : replace;
import std : sqrt;
import core.stdc.stdio : fread;
import vdmath;

void reflectDist(string[] args) {
    enforce(args.length == 5, "Need 4 (pixelX pixelY sampleCount reflectCount) arguments.");

    uint pixelX = args[1].to!uint;
    uint pixelY = args[2].to!uint;
    uint sampleCount = args[3].to!uint;
    uint reflectCount = args[4].to!uint;

    uint resolution = 400;
    float fov = 90;
    float Z = 70;

    // BRDF
    string[string] placeholders = [
        "FOV": fov.to!string, "Z": Z.to!string,
        "reflectCount": reflectCount.to!string,
        "RESOLUTION": resolution.to!string
    ];
    string cornellFileText = readText("cornellReflectDist.pbrt");
    foreach (key, value; placeholders)
        cornellFileText = cornellFileText.replace('<' ~ key ~ '>', value);
    std.file.write("cornellReflectDist2.pbrt", cornellFileText);

    string iridescenceArgs = format(
        "../scripts/cornellReflectDist2.pbrt --spp %u --pixel %u,%u", 1, pixelX, pixelY); // Don't need multiple samples for brdf distribution
    auto res2 = executeShell(
        `cd ../iridescence & cmake --build --preset cpuDebug & .\build\cpu\Debug\pbrt.exe ` ~ iridescenceArgs);
    enforce(res2.status == 0, "Iridescence Failed:\n" ~ res2.output);
    writeln(res2.output);
    writeln("Iridescence finished");

    // Mesh
    // Vec!3 distDir;
    // // Calculation different from pbrt calculation, import instead:
    // float pixelDelta = 2.0f * Z / resolution;
    // float deltaZ = -Z;
    // float deltaX = (pixelX - resolution / 2.0f + 0.5f) * pixelDelta;
    // float deltaY = (pixelY - resolution / 2.0f + 0.5f) * pixelDelta; // y orientation down
    // distDir = Vec!3(deltaX,deltaY,deltaZ).normalize();

    // Vec!3 wo;
    // File woBinary = File("../iridescence/distWo.bin", "rb");
    // fread(wo.vec.ptr, float.sizeof, 3, woBinary.getFP());
    // woBinary.close();
    // writeln("wo: ", wo);

    auto res1 = executeShell(
        "cd ../tracing & dub run --build=release --compiler=ldc2 -- reflectDist ../iridescence/distWo.bin " ~ sampleCount
            .to!string ~ ' ' ~ reflectCount.to!string); // Turns out I can also use Redirect.stdin
    enforce(res1.status == 0, "Tracing Failed:\n" ~ res1.output);
    writeln(res1.output);
    writeln("Tracing finished");

    auto res3 = executeShell("py dist.py");
    enforce(res3.status == 0, "PyPlot Failed:\n" ~ res3.output);
    writeln("PyPlot finished");
    writeln('\a'); // Ring bell
}
