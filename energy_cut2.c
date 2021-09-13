#include <TFile.h>
#include <TTree.h>
#include <stdio.h>
#include <vector>
#include <TH1D.h>
#include <stdlib.h>
#include <algorithm>


std::vector<unsigned int>	*brChannelID;
std::vector<float>	        *brEnergy;
std::vector<float>          *brQT1;


int energy_cut2(const char* in_file, const char* out_file, double max_energy)
{
    int i;

    //Open the data from the input root file/TTree
    TFile *foo = new TFile(in_file);
    TTree *T;
    T = (TTree*) foo->Get("data");
    
    //Access data and assign addresses
    unsigned int n = T->GetEntries();

    T->SetBranchAddress("energy", &brEnergy);
    
    //Create the output file for the energy cut by max hit energy
    
    TFile *nfmax = new TFile(out_file, "RECREATE");
    TTree *newtree_max;
    newtree_max = T->CloneTree(0);

    //Loop through each event
    for (i = 0; i < n; i++) {
        if (i % 100==0) {
            printf("\r%d/%d",i,n);
            fflush(stdout);
        }
        T->GetEntry(i);
        float max_hit_energy = *std::max_element(brEnergy->begin(), brEnergy->end());
        
        if (max_hit_energy >= max_energy) {
            newtree_max->Fill();
        }
    }
    printf("\r%d/%d\n",i,n);
    fflush(stdout);
    
    //Write out the file with energy cuts applied
    nfmax->Write();
    nfmax->Close();
    delete nfmax;

    foo->Close();
    delete foo;

    return 0;
}

int main(int argc, char* argv[]) {
    if (argc != 4) {
        printf("Three arguments are expected. Input file string, output file strings, and a double for the energy of the maximum hit event.");
        printf("%d arguments were supplied.\n", argc);
        exit(1);
    }
    energy_cut2(argv[1], argv[2], atof(argv[3]));
    return 0;
}

