#include <TTree.h>
#include <TFile.h>
#include <vector>
#include <stdio.h>

int branch_cutter(const char* in_file, const char* out_file)
{
    //Access the input file and the TTree to be trimmed
    TFile *input = TFile::Open(in_file);
    TTree *input_tree;
    input->GetObject("data", input_tree);

    //Open a new file that the new TTree will be written to
    TFile *output = TFile::Open(out_file, "RECREATE");

    //Clear the undesired branches
    input_tree->SetBranchStatus("prevEventFlags", 0);
    input_tree->SetBranchStatus("prevEventTime", 0);
    input_tree->SetBranchStatus("timeLast", 0);
    input_tree->SetBranchStatus("mh_n", 0);
    input_tree->SetBranchStatus("channelIdx", 0);
    input_tree->SetBranchStatus("tot", 0);
    input_tree->SetBranchStatus("t1coarse", 0);
    input_tree->SetBranchStatus("t1fine", 0);
    input_tree->SetBranchStatus("t2coarse", 0);
    input_tree->SetBranchStatus("t2fine", 0);
    input_tree->SetBranchStatus("qfine", 0);
    input_tree->SetBranchStatus("qT1", 0);
    input_tree->SetBranchStatus("qT2", 0);
    input_tree->SetBranchStatus("tacID", 0);
    input_tree->SetBranchStatus("qcoarse", 0);
 
    //Write the trimmed TTree to the output file
    TTree *output_tree = input_tree->CloneTree(-1,"fast");
    output->Write();

    output->Close();
    input->Close();

    delete output;
    delete input;

    return 0;
}

int main(int argc, char* argv[])
{
    if (argc != 3) {
	fprintf(stderr, "Two arguments are required, an input file name, and an output file name.");
        exit(1);
    }
    branch_cutter(argv[1], argv[2]);
    return 0;
}
