module reflectdist;
import std.process;
import std.stdio;
import std.file;
import std.math.trigonometry;
import std.exception : enforce;
import std.conv : to;
import std.format : format;
import std.array : replace;
import std : sqrt;
import core.stdc.stdio : fread;
import vdmath;
import vdmath.misc;
import std.math : PI;
import std.stdio: File;

void reflectDist(string[] args) {
    Vec!3 wo;

    float polar;
    float azimuth;
    switch (args.length) {
    case 5:
        polar = degreesToRadians(args[1].to!float);
        azimuth = degreesToRadians(args[2].to!float);
        wo.x = sin(polar) * cos(azimuth);
        wo.y = sin(polar) * sin(azimuth);
        wo.z = cos(polar);
        break;
    case 6:
        wo.x = args[1].to!float;
        wo.y = args[2].to!float;
        wo.z = args[3].to!float;
        polar = acos(wo.z);
        if (polar == 0)
            azimuth = 0;
        else
            azimuth = atan2(wo.y, wo.x);
        if (azimuth < 0)
            azimuth += PI;
        break;
    default:
        assert(0, "Need 4/5 arguments: (wo.x wo.y wo.z sampleCount reflectCount) / (polarAngle azimuthalAngle sampleCount reflectCount) but got: " ~ args
                .length.to!string);
    }
    wo = wo.normalize();
    uint sampleCount = args[$ - 2].to!uint;
    uint reflectCount = args[$ - 1].to!uint;
    string identifier = format("t%.2f_f%.2f", radiansToDegrees(polar), radiansToDegrees(azimuth));

    uint resolution = 400;
    float fov = 90;
    float Z = 70;
    uint pixelX = 200;
    uint pixelY = 200;

    string woString = wo.x.to!string ~ ' ' ~ wo.y.to!string ~ ' ' ~ wo.z.to!string;

    // BRDF
    string[string] placeholders = [
        "FOV": fov.to!string, "Z": Z.to!string,
        "REFLECT_COUNT": reflectCount.to!string,
        "RESOLUTION": resolution.to!string,
        "WO_VECTOR": "[ " ~ woString ~ " ]"
    ];
    for (int usePaul = 0; usePaul < 2; usePaul++) {
        string cornellFileText = readText("cornellReflectDist.pbrt");
        placeholders["USE_PAUL"] = (usePaul ? "true" : "false");
        placeholders["DIST_OUTPUT_FILE"] = "\"../scripts/temp/dist_" ~ identifier ~ (usePaul ? "_P"
                : "_L") ~ ".csv\"";
        foreach (key, value; placeholders)
            cornellFileText = cornellFileText.replace('<' ~ key ~ '>', value);
        string fileName = "temp/cornellReflectDist_" ~ identifier ~ (usePaul ? "_P" : "_L") ~ ".pbrt";
        std.file.write(fileName, cornellFileText);

        string iridescenceArgs = format(
            "../scripts/%s --spp %u --pixel %u,%u", fileName, 1, pixelX, pixelY); // Don't need multiple samples for brdf distribution
        auto res2 = executeShell(
            `cd ../iridescence & .\build\cpu\Release\pbrt.exe ` ~ iridescenceArgs);
        enforce(res2.status == 0, "Iridescence Failed:\n" ~ res2.output);
        writeln(res2.output);
    }
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

    // File exeFile = File("../tracing/tracing.exe");
    // exeFile.lock(); // try wait for issues? (race condition . . . but it actually works)
    // exeFile.unlock();

    auto res1 = executeShell(
        "cd ../tracing & dub run --build=debug -- reflectDist " ~ woString ~ ' ' ~ sampleCount
            .to!string ~ ' ' ~ reflectCount.to!string ~ ' ' ~ identifier); // Turns out I can also use Redirect.stdin
    enforce(res1.status == 0, "Tracing Failed:\n" ~ res1.output);
    writeln(res1.output);
    writeln("Tracing finished");

    auto res3 = executeShell("py dist.py " ~ identifier ~ ' ' ~ sampleCount.to!string);
    enforce(res3.status == 0, "PyPlot Failed:\n" ~ res3.output);
    writeln("PyPlot finished");
    writeln('\a'); // Ring bell
}
