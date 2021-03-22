#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <cassert>

#include <cpl_conv.h>
#include <gdal.h>
#include <gdal_priv.h>
#include <gdal_utils.h>

#include <ogr_api.h>
#include <ogr_srs_api.h>
#include <ogr_spatialref.h>

using namespace std;

struct bbox {
    double minx, miny, maxx, maxy;
    friend ostream& operator<<(ostream & stream, const bbox& box) {
        stream << std::fixed << box.minx << "," << box.miny << "," << box.maxx << "," << box.maxy;
        return stream;
    }
};

struct scene {
    string name;
    bbox box;
    friend ostream& operator<<(ostream & stream, const scene& s) {
        stream << "Name: " << s.name << endl << s.box;
        return stream;
    }
};

// Get the bounding box for a scene
static int get_coords(const string fname, bbox &box) {
    GDALDataset *poDS = reinterpret_cast<GDALDataset *>(GDALOpen(fname.c_str(), GA_ReadOnly));
    double gt[6];
    poDS->GetGeoTransform(gt);
    OGRSpatialReference proj(poDS->GetProjectionRef());
    if (!proj.IsGeographic()) { // Only geographic projection is adjusted, otherwise assume identity
        box.minx = gt[0];
        box.miny = gt[3] + gt[5] * poDS->GetRasterYSize();
        box.maxx = gt[0] + gt[1] * poDS->GetRasterXSize();
        box.maxy = gt[3];
        delete poDS;
        return 0;
    }

    OGRSpatialReference *llproj = proj.CloneGeogCS();

    // Make sure we convert to degrees
    //        if (llproj->GetRoot()->FindChild("UNIT"))
    llproj->GetRoot()->DestroyChild(llproj->GetRoot()->FindChild("UNIT"));

    OGRCoordinateTransformation *tr = OGRCreateCoordinateTransformation(&proj, llproj);
    // Start coords
    box.minx = gt[0];
    box.miny = gt[3] + gt[5] * poDS->GetRasterYSize();
    box.maxx = gt[0] + gt[1] * poDS->GetRasterXSize();
    box.maxy = gt[3];
    tr->Transform(1, &box.minx, &box.miny);
    tr->Transform(1, &box.maxx, &box.maxy);

    // This fails for scenes that cross the 0 meridian, which is possible
    //    assert(box.minx < box.maxx && box.miny < box.maxy);
    assert(box.miny < box.maxy);

    delete tr;
    delete llproj;
    delete poDS;
    return 0;
}

// Convert a file list to a vector of scenes, which include bounding boxes
// the prefix is pre-pended to each line, in can be empty
static int get_info(const string path, const string flist, vector<scene> &results) {
    ifstream f(path + flist);
    scene s;
    if (!f.is_open())
        return 1;
    while (getline(f, s.name) && s.name.size() != 0) {
        s.name = path + s.name;
        get_coords(s.name, s.box);
        results.push_back(s);
    }
    return 0;
}

static bool contact(const bbox &b1, const bbox &b2) {

    if (b2.minx < b2.maxx)
        return !(b1.minx > b2.maxx || b1.maxx < b2.minx || b1.miny > b2.maxy || b1.maxy < b2.miny);

    // B2 is the scene, which crosses the antimeridian
    // Split it into two and check each side
    bbox b = b2;
    b.minx = -180;

    if (contact(b1, b))
        return true;
    b = b2;
    b.maxx = 180;

    return contact(b1, b);
}

static bool contact(const bbox &b, const vector<scene> &scenes) {
    for (auto s : scenes)
        if (contact(b, s.box))
            return true;
    return false;
}

// Apply geo-transform to inbox
static int convert(const double *gt, const bbox &inbox, bbox &outbox) {
    outbox.minx = gt[0] + inbox.minx * gt[1];
    outbox.maxx = gt[0] + inbox.maxx * gt[1];
    // These are swapped because gt[5] is negative
    outbox.maxy = gt[3] + inbox.miny * gt[5];
    outbox.miny = gt[3] + inbox.maxy * gt[5];
    return 0;
}

// Test against the noData
template<typename T> bool buffer_empty(const T *buffer, size_t sz, T ndv) {
    while (sz--)
        if (*buffer++ != ndv)
            return false;
    return true;
}

// The box is in longlat coordinates
// call with pOut = nullptr when done, to release static buffer
static int generateRegion(GDALDataset *pOut, const bbox &box, const vector<scene> scenes) {
    // These are semi-permanent
    static size_t bufsize = 0;
    static void *buffer = nullptr;

    if (!pOut) { // If called with no output, release buffer and quit
        CPLFree(buffer);
        buffer = nullptr;
        bufsize = 0;
        return 0;
    }

    int uerrerr;
    int xsize(pOut->GetRasterXSize()), ysize(pOut->GetRasterYSize());
    CPLString XSize, YSize;
    XSize.Printf("%d", xsize);
    YSize.Printf("%d", ysize);

    // Direct and inverse geotransforms
    double gt[6], igt[6];
    pOut->GetGeoTransform(gt);
    GDALInvGeoTransform(gt, igt);
    double ndv = pOut->GetRasterBand(1)->GetNoDataValue();

    CPLString Left(CPLOPrintf("%f", gt[0]));
    CPLString Right(CPLOPrintf("%f", gt[0] + gt[1] * xsize));
    CPLString Top(CPLOPrintf("%f", gt[3]));
    CPLString Bottom(CPLOPrintf("%f", gt[3] + gt[5] * ysize));

    vector<const char *> warpopts;
#define OPT(X) warpopts.push_back(X)
    OPT("-of"); OPT("VRT");
    OPT("-t_srs");

    if (pOut->GetProjectionRef())
        OPT(pOut->GetProjectionRef());
    else 
        OPT("+proj=longlat");

    OPT("-r"); OPT("cubic");
    OPT("-ts"); OPT(XSize); OPT(YSize);
    OPT("-te"); OPT(Left); OPT(Bottom); OPT(Right); OPT(Top);
    OPT(NULL);
#undef OPT

    GDALWarpAppOptions *wo = GDALWarpAppOptionsNew(const_cast<char **>(warpopts.data()), NULL);

    // Warp each individual input to the output coordinates
    vector<char *> vrtfnames;

    for (int i = 0; i < scenes.size(); i++) {

        GDALDatasetH poDS = GDALOpen(scenes[i].name.c_str(), GA_ReadOnly);

        vrtfnames.push_back(strdup(CPLOPrintf("/vsimem/tmp%d.vrt", i)));
        GDALDatasetH hOutDS = GDALWarp(vrtfnames[i], NULL, 1, &poDS, wo, &uerrerr);
        GDALClose(hOutDS);
    }

    GDALWarpAppOptionsFree(wo);

    // If there is a single input, use that vrt, otherwise build one that combines all inputs
    if (scenes.size() > 1) {
        vector<const char *> cropopts;
#define OPT(X) cropopts.push_back(X)
        OPT(NULL);
#undef OPT
        GDALBuildVRTOptions *vo = GDALBuildVRTOptionsNew(const_cast<char **>(cropopts.data()), NULL);
        vrtfnames.push_back(strdup("/vsimem/tmp.vrt"));

        GDALDatasetH hOutDS = GDALBuildVRT(vrtfnames.back(), vrtfnames.size() - 1, 
            NULL, vrtfnames.data(), vo, &uerrerr);

        GDALBuildVRTOptionsFree(vo);
        GDALClose(hOutDS);
    }

    // The global VRT or the single input VRT is the last in the vrtfnames
    GDALDataset *pIn = reinterpret_cast<GDALDataset *>(GDALOpen(vrtfnames.back(), GA_ReadOnly));

    // Copy the region, in 1024x1024 chunks
    bbox pixbox;
    convert(igt, box, pixbox);
    // This delta is in pixels
    int delta = 1024;
    // cerr << pixbox << " " << scenes.size() << endl;
    assert(int(pixbox.maxx - pixbox.minx) % delta == 0 && int(pixbox.maxy - pixbox.miny) % delta == 0);
    GDALDataType dt = pOut->GetRasterBand(1)->GetRasterDataType();
    size_t bsize = delta * delta * GDALGetDataTypeSizeBytes(dt);

    if (bsize > bufsize) {
        // Extend the buffer as needed. Deallocate by calling generateRegion(nullptr,...)
        bufsize = bsize;
        buffer = CPLRealloc(buffer, bsize);
    }

    CPLErr err;

    for (int y = pixbox.miny; y < pixbox.maxy; y += delta) {
        // This will count down to zero within a delta square
        // cerr << "\r     \r" << int((pixbox.maxy - y) / delta);
        // Try the whole row, check for intersections bbox
        bbox tbox = {pixbox.minx, double(y), pixbox.maxx, double(y) + delta};
        bbox lltbox;
        convert(gt, tbox, lltbox);
        if (!contact(lltbox, scenes))
            continue;

        for (int x = pixbox.minx; x < pixbox.maxx; x += delta) {
            // Check for this smaller box too, the Y is already set from above
            tbox.minx = x;
            tbox.maxx = x + delta;
            convert(gt, tbox, lltbox);
            if (!contact(lltbox, scenes))
                continue;

            err = pIn->RasterIO(GF_Read, x, y, delta, delta, buffer, delta, delta, dt,
                1, NULL, 0, 0, 0);

            if (err != CE_None) {
                cerr << "Error reading from the VRT" << endl;
                for (int i = 0; i < scenes.size(); i++)
                    cerr << i << " " << scenes[i].name << endl;
            }

            bool empty;

            switch (dt) {

#define DISPATCH(Tname, T)\
    case Tname:\
        empty = buffer_empty(reinterpret_cast<T *>(buffer), delta*delta, static_cast<T>(ndv));\
        break;

            DISPATCH(GDT_Byte, GByte);
            DISPATCH(GDT_Int16, GInt16);
            DISPATCH(GDT_UInt16, GUInt16);
            DISPATCH(GDT_Int32, GInt32);
            DISPATCH(GDT_UInt32, GUInt32);
            DISPATCH(GDT_Float32, float);

#undef DISPATCH

            default:
                assert(false); // Unimplemented
            }

            if (!empty) {
                err = pOut->RasterIO(GF_Write, x, y, delta, delta, buffer, delta, delta, dt,
                    1, NULL, 0, 0, 0);
                if (err != CE_None)
                    cerr << "Error writing to output" << endl;
            }
        }
    }
    cout << "\r"; // Go to the start of the line
    // This closes all dependent datsets
    GDALClose(pIn);

    // Clean up memory virtual files
    for (auto fn : vrtfnames) {
        VSIUnlink(fn);
        free(fn);
    }

    // Not required, but doesn't hurt either
    pOut->FlushCache();
    return 0;
}

int Usage(char **argv) {
    cerr << "Builds a single global image from a collection of inputs" << endl;
    cerr << "It reprojects the data from input to global coordinates (longlat)" << endl;
    cerr << "Best used with an MRF output, it will skip areas with no input" << endl;
    cerr << endl;
    cerr << "Usage: " << argv[0] << " [options] in_path file_of_names out" << endl;
    cerr << "Options: -delta D" << endl;
    cerr << "         D is the chunk size generated, in degrees, defaults to 1.125";
    cerr << "Options: -ndv V" << endl;
    cerr << "         V is the no data value, should match the input, defaults to 1024*1024";
    cerr << "Options: -xsize N" << endl;
    cerr << "         Y size is half of X, X defaults to 10 * 2^20." << endl;
    cerr << "         N should chosed such that N * 360 / delta / 512 is an integer" << endl;
    cerr << "Options: -ysize N" << endl;
    cerr << "         Y size, if not half of X" << endl;
    cerr << "Options: -proj PRJSTR" << endl;
    cerr << "         PRJSTR is in a format understood by gdalwarp.  Defaults to +proj=latlong" << endl;
    cerr << "Options: -bb bblx bbly bbux bbuy" << endl;
    cerr << "         Bounding box in real world coordinates.  Defaults to -180 -90 180 90" << endl;
    cerr << "Options: -te xmin ymin xmax ymas" << endl;
    cerr << "         The extent which is written. The output MRF still has global span." << endl;
    cerr << "Options: -start x y" << endl;
    cerr << "         Seed the main loops with the x and y. They have to be valid values within the extent, otherwise they get ignored" << endl;
    cerr << "Options: -stop x y" << endl;
    cerr << "         Stop when x and y of the main loops reach these values. They have to be valid values within the extent, otherwise it stops when it reaches the extent end" << endl;
    cerr << "Options: -of <GDAL_FORMAT>" << endl;
    cerr << "         Output format, defaults to MRF" << endl;
    cerr << "Options: -ot <GDAL_Type>" << endl;
    cerr << "         Output type, should match the input type, defaults to Byte" << endl;
    cerr << "Options: -co <Key>=<Val>" << endl;
    cerr << endl;
    cerr << "For DEM: -of MRF -co COMPRESS=LERC -co OPTIONS=\"LERC_PREC=0.1 DEFLATE=1 GZ=1 V1=1\"" << endl;
    cerr << "For Image: -of MRF -ndv 0 -co COMPRESS=JPEG" << endl;
    return 1;
}

int main(int argc, char **argv)
{
    GDALAllRegister();
    // Matching arrays for input names and bounding boxes in lat-lon
    vector<scene> scenes;

    // Max scenes per region without a warning
    size_t maxspr = 30;

    // Pixel count
    int pcountx = 2 * 10 * 1024 * 1024; // About 1m on mars, x size
    int pcounty = 0; // Flag, update later
    int dryrun = false; // Flag, don't produce output
    char **options = NULL;
    const char *drivername = "MRF";
    const char *proj = nullptr;
    double delta = 1.125;
    double ndv = 1024 * 1024;

    // Bounding box for the whole output, lower and upper x and y
    double bblx = -180, bbly = -90, bbux = 180, bbuy = 90;

    // The bounds of the area that will be written
    int te = 0; // Flag, values were set by user
    double lowx = bblx, lowy = bbly, highx = bbux, highy = bbuy;

    GDALDataType dt = GDT_Byte;

    // -1000 is a flag that they didn't get set by the user
    double startx = -1000, starty = -1000, stopx = -1000, stopy = -1000;

    //options = CSLAddNameValue(options, "COMPRESS", "LERC");
    //options = CSLAddNameValue(options, "OPTIONS", "LERC_PREC=0.1 DEFLATE=1 GZIP=1");

    int i = 1;
    if (argc <= i) return Usage(argv);
    while (argv[i][0] == '-') {
        if (!strcmp(argv[i], "-ndv")) {
           ndv = atof(argv[++i]);
        }
        else
        if (!strcmp(argv[i], "-bb")) {
            bblx = atof(argv[++i]);
            bbly = atof(argv[++i]);
            bbux = atof(argv[++i]);
            bbuy = atof(argv[++i]);
        }
        else
        if (!strcmp(argv[i], "-delta")) {
            delta = atof(argv[++i]);
        }
        else
        if (!strcmp(argv[i], "-dryrun")) {
            dryrun = true;
        }
        else
        if (!strcmp(argv[i], "-proj")) {
            proj = argv[++i];
        }
        else
        if (!strcmp(argv[i], "-of")) {
            drivername = argv[++i];
        }
        else
        if (!strcmp(argv[i], "-wspr")) {
            maxspr = atoi(argv[++i]);
        }
        else
        if (!strcmp(argv[i], "-te")) {
            te = 1;
            lowx = atof(argv[++i]);
            lowy = atof(argv[++i]);
            highx = atof(argv[++i]);
            highy = atof(argv[++i]);
        }
        else
        if (!strcmp(argv[i], "-start")) {
            startx = atof(argv[++i]);
            starty = atof(argv[++i]);
        }
        else
        if (!strcmp(argv[i], "-stop")) {
            stopx = atof(argv[++i]);
            stopy = atof(argv[++i]);
        }
        else
        if (!strcmp(argv[i], "-co")) {
            i++;
            char *eqpos = strchr(argv[i], '=');
            if (!eqpos) {
                cerr << "Create option syntax is key=val" << endl;
                return Usage(argv);
            }
            string key(argv[i]), val(eqpos + 1);
            key.resize(eqpos - argv[i]);
            options = CSLAddNameValue(options, key.c_str(), val.c_str());
        }
        else
        if (!strcmp(argv[i], "-xsize")) {
            pcountx = atoi(argv[++i]);
        }
        else
        if (!strcmp(argv[i], "-ysize")) {
            pcounty = atoi(argv[++i]);
        }
        else
        if (!strcmp(argv[i], "-ot")) {
            dt = GDALGetDataTypeByName(argv[++i]);
            if (dt == GDT_Unknown) {
                cerr << "Unknown data type " << argv[i] << endl;
                return Usage(argv);
            }
        }
        else {
            cerr << "Unknown option " << argv[i] << endl;
            return Usage(argv);
        }
        i++;
    }

    if (argc - i < 3)
        return Usage(argv);

    if (pcounty == 0)
        pcounty = pcountx / 2;

    if (!te) {
        lowx = bblx;
        lowy = bbly;
        highx = bbux;
        highy = bbuy;
    }

    string path(argv[i++]);
    string in_list(argv[i++]);
    string ofname(argv[i++]);

    if (*path.rbegin() != '/' and *path.rbegin() != '\\')
        path += '/';

    if (get_info(path, in_list, scenes)) {
        cerr << "Can't open " << path + in_list << endl;
        return 1;
    }

    GDALDriver *mrfdriver = reinterpret_cast<GDALDriver *>(GDALGetDriverByName(drivername));
    GDALDataset *pDS = mrfdriver->Create((path + ofname).c_str(), pcountx, pcounty, 1, dt, options);
    if (proj)
        pDS->SetProjection(proj);
    GDALRasterBand *b = pDS->GetRasterBand(1);
    // Set this same as the input
    b->SetNoDataValue(ndv);

    double gt[6] = { bblx, (bbux - bblx) / pcountx, 0 , bbuy, 0, -(bbuy - bbly)/ pcounty };
    pDS->SetGeoTransform(gt);

    // If the start point was not set, use the extent minimums
    if (startx == -1000)
        startx = lowx;
    if (starty == -1000)
        starty = lowy;

    bbox region;
    bool quit_now = false;
    for (region.miny = starty; region.miny < highy; region.miny += delta) {
        region.maxy = region.miny + delta;
        for (region.minx = lowx; region.minx < highx; region.minx += delta) {
            if (lowx != startx) { // Only done in the first xy loop
                region.minx = startx;
                startx = lowx;
            }

            // Check the stop condition
            quit_now = region.miny == stopy && region.minx == stopx;
            if (quit_now) break; // the x loop

            region.maxx = region.minx + delta;
            // Select scenes that intersect this region
            vector<scene> subscenes;
            for (auto s : scenes)
                if (contact(region, s.box))
                    subscenes.push_back(s);


            // Generate each chunk that has input scenes
            if (!subscenes.empty()) {
                cerr << "Starting: " << region.minx << " " << region.miny << endl;
                // This could be a variable value to warn on
                if (subscenes.size() > maxspr) {
                    cerr << "Warning: " << subscenes.size() << " scenes used for this region" << endl;
                    if (dryrun) {
                        cerr << "Scenes: " << endl;
                        for (auto s : subscenes)
                            cerr << s.name << endl;
                    }
                }

                if (!dryrun)
                    generateRegion(pDS, region, subscenes);
            }
        }
        if (quit_now) break; // the y loop
    }

    delete pDS;
    // Call to deallocate region buffer, only first parameter counts
    generateRegion(nullptr, region, scenes);
}
