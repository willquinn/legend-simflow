//usr/bin/env root -x -q -l ${0}\(\""${0}"\",\""${*}"\"\); exit $?

#include <getopt.h>

struct mini_shroud {
    double height;
    double radius;
    double center_x;
    double center_y;
    double center_z;
};

// volume of lar-cylinder-K42 simulated region (estimated by MaGe): 3.72495e+06 cm^3
// TODO: update values
double in_fraction = 0.00195451;
double out_fraction = 0.99804549;

std::vector<mini_shroud> mss = {
    {956  , 51,  109   ,  188.79, -25},
    {960  , 48,  198.22,   95.44, -27},
    {652  , 44,  220   ,   0    , 127},
    {898  , 48,  198.22,  -95.44,   4},
    {543  , 54,  109   , -188.79, 181.5},
    {795  , 51, -109   , -188.79,  55.5},
    {959.5, 48, -198.22,  -95.44, -26.75},
    {966.5, 44, -220   ,    0   , -30.25},
    {951.5, 48, -198.22,   95.44, -22.75},
    {913  , 51, -109   ,  188.79,  -3.5},
};


bool point_is_inside_MS(double x, double y, double z, mini_shroud ms) {
    // translate coordinate system
    double tr_x = x - ms.center_x;
    double tr_y = y - ms.center_y;
    double tr_z = z - ms.center_z;

    // height inside MS AND radius inside MS
    return (abs(tr_z) <= ms.height/2) and (sqrt(tr_x * tr_x + tr_y * tr_y) <= ms.radius);
}


bool point_is_inside_any_MS(double x, double y, double z) {
    for (auto& ms: mss) {
        if (point_is_inside_MS(x, y, z, ms)) {
            return true;
        }
    }
    return false;
}


void usage() {
    std::cerr << "\n"
              << "USAGE: split_K42_vertices [options] input_file output_file_in output_file_out\n"
              << "\n"
              << "options:\n"
              << "  --help|-h : print this help message and exit\n"
              << "\n";
    gSystem->Exit(1);
}


void split_K42_vertices(std::string prog = "split_K42_vertices", std::string args = "") {

    // this is for getopt to work
    args = prog + " " + args;

    int argc = 0;
    char** argv = new char*[50];

    // get all arguments
    std::istringstream iss(args);
    std::string word;
    while (iss >> word) {
        char* tmp = new char[50];
        strcpy(tmp, word.c_str());
        argv[argc] = tmp;
        argc++;
    }

    const char* const short_opts = ":h";
    const option long_opts[] = {
        {"help", no_argument, nullptr, 'h'},
        {nullptr, no_argument, nullptr, 0}
    };

    // read in with getopt
    int opt = 0;
    while ((opt = getopt_long(argc, argv, short_opts, long_opts, nullptr)) != -1) {
        switch (opt) {
            case 'h': // -h or --help
            case '?': // Unrecognized option
            default:
                usage();
        }
    }

    // get extra arguments
    std::vector<std::string> extra_args;
    for(; optind < argc; optind++){
        extra_args.emplace_back(argv[optind]);
    }

    if (extra_args.size() != 3) {
        usage();
    }

    auto infile = extra_args[0];
    auto outfile_in = extra_args[1];
    auto outfile_out = extra_args[2];

    // open file
    TFile file(infile.c_str(), "READ");
    TTree * fTree = dynamic_cast<TTree*>(file.Get("fTree"));

    // create output files and copy original tree twice
    TFile file_of_out(outfile_out.c_str(), "RECREATE");
    TTree * fTree_out = dynamic_cast<TTree*>(fTree->CloneTree(0));
    TFile file_of_in(outfile_in.c_str(), "RECREATE");
    TTree * fTree_in  = dynamic_cast<TTree*>(fTree->CloneTree(0));

    // get vertex position
    MGTMCEventSteps* evt = nullptr;
    fTree->SetBranchAddress("eventPrimaries", &evt);

    // FIXME: this way of changing the number of simulated events does not work, somehow
    // unsigned long fNEvents = 0;
    // fTree->SetBranchAddress("fNEvents", &fNEvents);

    int nevents = fTree->GetEntries();
    int n_inside = 0;

    // loop events and write them in the right tree
    for(int e = 0; e < nevents; e++) {
        fTree->GetEntry(e);
        auto v = evt->GetStep(0)->GetPositionVector();

        if (point_is_inside_any_MS(v.X(), v.Y(), v.Z())) {
            // fNEvents *= in_fraction;
            fTree_in->Fill();
            n_inside++;
        }
        else {
            // fNEvents *= out_fraction;
            fTree_out->Fill();
        }
    }

    std::cout << "INFO: total entries: " << nevents << " of which "
              << n_inside << " inside the mini-shrouds" << std::endl;

    file_of_out.Write();
    file_of_in.Write();

    return;
}
