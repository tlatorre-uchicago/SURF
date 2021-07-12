#include <TFile.h>
#include <TTree.h>
#include <stdio.h>
#include <vector>
#include <TH1D.h>
#include <stdlib.h>
#include <algorithm>


std::vector<unsigned int>	*brChannelID;
std::vector<float>	        *brEnergy;


int energy_cut2(const char* in_file, const char* out_file_sum, const char* out_file_max, const char* out_file_avg, double sum_energy, double max_energy, double avg_energy) {
    
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
    
    //Create the output files for the energy cut by total event energy
    //and by average hit energy in an event
    TFile *nfsum = new TFile(out_file_sum, "RECREATE");
    TTree *newtree_sum;
    newtree_sum = T->CloneTree(0);
    
    TFile *nfmax = new TFile(out_file_max, "RECREATE");
    TTree *newtree_max;
    newtree_max = T->CloneTree(0);
    
    TFile *nfa = new TFile(out_file_avg, "RECREATE");
    TTree *newtree_avg;
    newtree_avg = T->CloneTree(0);
    
    int hits = 0;
    double hits_eliminated_sum = 0;
    double hits_eliminated_avg = 0;
    
    TH1D *h = new TH1D("h", "energy cut using mins", 500,0,1500);
    
    //Loop through each event
    for (int i=0; i<n; i++) {
        if (i % 100==0) {
            printf("\r%d/%d",i,n);
            fflush(stdout);
        }
        T->GetEntry(i);
        float event_energy = 0;
        float max_hit_energy = *std::max_element(brEnergy->begin(), brEnergy->end());
        
        //Loop through the hits in an event
        for (int j=0; j<brChannelID->size(); j++) {
            event_energy+=brEnergy->at(j);
            hits++;
        }
        
        //If the total energy exceeds the input minimum, add the hits
        //to the new TTree file
        if (event_energy >= sum_energy) {
            newtree_sum->Fill();
            h->Fill(event_energy);
        }
        if (event_energy/(brChannelID->size()) >= avg_energy) {
            newtree_avg->Fill();
        }
        
        if (max_hit_energy >= max_energy) {
            newtree_max->Fill();
        }
    }
    h->Draw();
    
    //Write out the files with energy cuts applied
    nfsum->Write();
    nfmax->Write();
    nfa->Write();
    
    /*
    printf("\nThe fraction of hits eliminated by min sum cut: %f\n", hits_eliminated_min/hits);
    printf("The fraction of hits eliminated by avg cut: %f\n", hits_eliminated_avg/hits);
    printf("There are %d events.\n There were %d bad min sum events and %d bad avg events", n, count_min, count_avg);
    */
    return 0;
}

int main(int argc, char* argv[]) {
    if (argc != 7) {
        printf("Six arguments are expected. Input file string, "
    "three output file strings, and three doubles for minimum event "
    " energy, in events maximum hit energy per event cutoffs respectively, and average hit energy   in an event.\n");
        printf("%d arguments were supplied.\n", argc);
        exit(1);
    }
    energy_cut2(argv[1], argv[2], argv[3], argv[4], atof(argv[5]), atof(argv[6]), atof(argv[7]));
    return 0;
}

