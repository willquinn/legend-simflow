//usr/bin/env root -l -x ${0}\(\""${0}"\",\""${*}"\"\); exit $?

#include "getopt.h"

void plot_mage_vertices(std::string prog = "plot_mage_vertices", std::string args = "") {

    const auto usage = [&prog](){
        std::cerr << "\nUSAGE: " + prog + " [-o OUTPUT_FILE] FILE [FILE...]\n";
        gSystem->Exit(1);
    };

    // this is for getopt to work
    args = prog + " " + args;

    int argc = 0;
    char** argv = new char*[200];

    // get all arguments
    std::istringstream iss(args);
    std::string word;
    while (iss >> word) {
        char* tmp = new char[200];
        strcpy(tmp, word.c_str());
        argv[argc] = tmp;
        argc++;
    }

    const char* const short_opts = "o:bh";
    const option long_opts[] = {
        {"output",     required_argument, nullptr, 'o'},
        {"batch",      no_argument,       nullptr, 'b'},
        {"help",       no_argument,       nullptr, 'h'},
        {nullptr,      no_argument,       nullptr, 0  }
    };

    // defaults
    std::string output;
    bool batch = false;

    // read in with getopt
    int opt = 0;
    while ((opt = getopt_long(argc, argv, short_opts, long_opts, nullptr)) != -1) {
        switch (opt) {
            case 'o':
                output = optarg;
                break;
            case 'h': // -h or --help
            case 'b':
                batch = true;
                break;
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

    if (extra_args.empty()) usage();

    // do it!

    // do not display any of the standard histogram decorations
    gStyle->SetOptStat(false);
    gStyle->SetOptFit(0);
    gROOT->SetBatch(batch);

    auto tree = new TChain("fTree");
    for (const auto& f : extra_args) tree->Add(f.c_str());

    auto c = new TCanvas("c", "MaGe event vertices", 1000, 1000);
    c->Divide(2, 2);
    c->cd(1);
    tree->Draw("eventPrimaries.fSteps.fY:eventPrimaries.fSteps.fX");
    c->cd(2);
    tree->Draw("eventPrimaries.fSteps.fZ:eventPrimaries.fSteps.fY");
    c->cd(3);
    tree->Draw("eventPrimaries.fSteps.fZ:eventPrimaries.fSteps.fX");
    c->cd(4);
    tree->Draw("eventPrimaries.fSteps.fZ:eventPrimaries.fSteps.fY:eventPrimaries.fSteps.fX");

    if (!output.empty()) c->SaveAs(output.c_str());

    if (batch) gSystem->Exit(0);
}
