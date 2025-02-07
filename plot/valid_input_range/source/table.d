module table;

import std.conv;
import std.math.constants : PI, PI_2;
import std.algorithm.comparison : min, max;
import std.stdio;
import std.algorithm.searching;
import std.string;

struct BoundsThetaPhi {
    double min;
    double max;
    double stepSize;
}

class BRDF {

    double[] backBounceTable;
    BoundsThetaPhi thetaBounds;
    BoundsThetaPhi phiBounds;
    uint thetaDimension;
    uint phiDimension;

    double readTable(uint thetaIndex, uint phiIndex) {
        if (thetaIndex >= thetaDimension || phiIndex > phiDimension)
            assert(0, "Backbounce table of dimensions [$(thetaDimension)][$(phiDimension)] accessed at index [$(thetaIndex)][$(phiIndex)]."
                    .text);
        if (phiIndex * thetaDimension + thetaIndex >= backBounceTable.length)
            assert(0, "Index out of range."); // should be redundant.
        return backBounceTable[phiIndex * thetaDimension + thetaIndex];
    }

    double linInterpolate(double val1, double val2, double factor) {
        return val1 * factor + val2 * (1 - factor);
    }

    double linInterpolate2D(double ll, double hl, double lh, double hh, double xfactor, double yfactor) {
        double bottom = linInterpolate(ll, lh, xfactor);
        double top = linInterpolate(hl, hh, xfactor);
        return linInterpolate(bottom, top, yfactor);
    }

    // Find value by trilinear interpolation
    // Valid input: theta [-pi_2, pi_2], phi [-2pi, 4pi]
    double readTable(double theta, double phi) {
        // Transform input to [0, pi_2], [0, pi] range
        if (theta < 0) { // no negative theta
            theta = -theta;
            phi += PI;
        }
        if (phi < 0) // no negative phi
            phi += 2 * PI;
        phi %= 2 * PI; // collapse to [0,2pi>
        if (phi > PI) // Mirror phi as backBounce is symmetric.
            phi = 2 * PI - phi;

        if (theta < thetaBounds.min || theta > thetaBounds.max ||
            phi < phiBounds.min || phi > phiBounds.max) {
            writeln(i"Invalid table read coordinates ($(theta), $(phi)) outside bounds [$(thetaBounds.min), $(thetaBounds.max)],[$(phiBounds.min), $(phiBounds.max)]");
            // return std::numeric_limits<double>::signaling_NaN();
        }

        // Calculate index  coordinates.
        // Floating point error prevention by use of max/min
        double thetaDiff = max(0.0, theta - this.thetaBounds.min);
        double phiDiff = max(0.0, phi - this.phiBounds.min);
        uint thetaLow = min(cast(uint)(thetaDiff / this.thetaBounds.stepSize), this.thetaDimension - 1);
        uint phiLow = min(cast(uint)(phiDiff / this.phiBounds.stepSize), this.phiDimension - 1);

        uint thetaHigh = min(thetaLow + 1, this.thetaDimension - 1); // clamp
        uint phiHigh = min(phiLow + 1, this.phiDimension - 1); // clamp

        double ll = readTable(thetaLow, phiLow);
        double hl = readTable(thetaHigh, phiLow);
        double lh = readTable(thetaLow, phiHigh);
        double hh = readTable(thetaHigh, phiHigh);
        double xfactor = (thetaDiff - (thetaLow * this.thetaBounds.stepSize)) / this
            .thetaBounds.stepSize;
        double yfactor = (phiDiff - (phiLow * this.phiBounds.stepSize)) / this.phiBounds.stepSize;
        return linInterpolate2D(ll, hl, lh, hh, xfactor, yfactor);
    }

    // ----------------------------------

    this(string filename) {
        File tableFile = File(filename, "rb");
        string line = tableFile.readln();
        if (line[0 .. 13] != "double double")
            assert(0, "File header wrong types:\n" ~ line);

        uint thetaBoundsLength;
        uint phiBoundsLength;
        line = tableFile.readln();
        line = line.stripLeft();
        thetaBoundsLength = parse!uint(line);
        line = line.stripLeft();
        phiBoundsLength = parse!uint(line);

        uint thetaDimension;
        uint phiDimension;

        line = tableFile.readln();
        line = line.stripLeft();
        thetaDimension = line.parse!uint();
        line = line.stripLeft();
        phiDimension = line.parse!uint();
        this.thetaDimension = thetaDimension;
        this.phiDimension = phiDimension;

        BoundsThetaPhi thetaBounds;
        BoundsThetaPhi phiBounds;

        thetaBounds = tableFile.rawRead(new BoundsThetaPhi[1])[0];
        phiBounds = tableFile.rawRead(new BoundsThetaPhi[1])[0];
        if (thetaBounds.min != 0 || (thetaBounds.max - PI_2) > 1e-8)
            assert(0, "thetaBounds not of range [0,pi_2]: [" ~ thetaBounds.min.to!string ~ "," ~ thetaBounds
                    .max.to!string ~ "]");
        if (phiBounds.min != 0 || (phiBounds.max - PI) > 1e-8)
            assert(0, "phiBounds not of range [0,pi]: [" ~ phiBounds.min.to!string ~ "," ~
                    phiBounds.max.to!string ~ "]");
        this.thetaBounds = thetaBounds;
        this.phiBounds = phiBounds;

        size_t tableSize = thetaDimension * phiDimension;
        // backBounceTable.resize(tableSize);
        double[] backBounceTable = tableFile.rawRead(new double[tableSize]);
        this.backBounceTable = backBounceTable;

        // double[][] reformed = new double[][phiDimension];
        // foreach (i, ref part; reformed) {
        //     part = backBounceTable[i * thetaDimension .. i * thetaDimension + thetaDimension];
        //     part.writeln;
        // }
    }
}
