#include <TFile.h>
#include <TTree.h>
#include <stdio.h>
#include <vector>
#include <TH1D.h>
#include <stdlib.h>
#include <algorithm>


std::vector<unsigned int>	*brChannelID;
std::vector<float>	        *brEnergy;


int energy_cut2(const char* in_file, const char* out_file_max, double max_energy) {
    
    //Open the data from the input root file/TTree
    TFile *foo = new TFile(in_file);
    TTree *T;
    T = (TTree*) foo->Get("data");
    
    //Access data and assign addresses
    unsigned int n = T->GetEntries();

    T->SetBranchAddress("channelID", &brChannelID);
    T->SetBranchAddress("energy", &brEnergy);
    
    
    unsigned int channelID;
    float energy;
    int count_sum = 0;
    int count_avg = 0;
    
    //Create the output file for the energy cut by max hit energy
    
    TFile *nfmax = new TFile(out_file_max, "RECREATE");
    TTree *newtree_max;
    newtree_max = T->CloneTree(0);

    int hits = 0;

    //Loop through each event
    for (int i=0; i<n; i++) {
        if (i % 100==0) {
            printf("\r%d/%d",i,n);
            fflush(stdout);
        }
        T->GetEntry(i);
        float event_energy=0;
        float max_hit_energy = *std::max_element(brEnergy->begin(), brEnergy->end());
        
        //Loop through the hits in an event
        for (int j=0; j<brChannelID->size(); j++) {
            event_energy+=brEnergy->at(j);
            hits++;
        }
        
        if (event_energy/(brChannelID->size()) >= avg_energy) {
            newtree_avg->Fill();
        }
    }
    h->Draw();
    
    //Write out the file with energy cuts applied
    nfmax->Write();
    */
    return 0;
}

int main(int argc, char* argv[]) {
    if (argc != 4) {
        printf("Three arguments are expected. Input file string, output file strings, and a double for the energy of the maximum hit.m event.);
        printf("%d arguments were supplied.\n", argc);
        exit(1);
    }
    energy_cut2(argv[1], argv[2], atof(argv[3]));
    return 0;
}

