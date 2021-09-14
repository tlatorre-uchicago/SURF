#include <vector>
#include <stdio.h>
#include <TH1D.h>
#include <TTree.h>
#include <TList.h>
#include <TFile.h>
#include <iostream>
#include <sstream>
#include <algorithm>

/*
ROOT macro to read in a file produced by ./convert_raw_event and analyze it.
Date: 2021-6-24.
*/

int brChannelIdx[128];
/* Time is in ps. */
std::vector<long long>        	*brTime;
std::vector<unsigned int>	*brChannelID;
/* In arbitrary ADC units. */
std::vector<float>	        *brEnergy;


void all_sipm_pair_time_differences(const char* in_file, const char* out_file, double energy_threshold) {
    
    /*Retrieve the filename
    char filename[100];    
    cout << "Please Enter a raw_to_event filename: ";
    cin >> filename;
    */
    
    TFile *foo = new TFile(in_file);
    TTree *T;
    T = (TTree *) foo->Get("data");

    T->SetBranchAddress("time", &brTime);
    T->SetBranchAddress("channelID", &brChannelID);
    T->SetBranchAddress("energy", &brEnergy);
    
    int time_first;
    int time_second;
    
    // Number of events in the TTree. 
    unsigned int n = T->GetEntries();
    
    //Make 16 histograms, one for each channel pair
    int pairs[16][2];
    TH1D *tdh[16];   
    char histoname[256];
    char histotitle[256];
    for (int k=0; k<16; k++) {
        sprintf(histoname, "h%d", k);
        sprintf(histotitle, "Time Difference for channels %d and %d", k + 64, 31 - k + 64);
        tdh[k] = new TH1D(histoname, histotitle, 300, -2000, 2000);
    }
    
    //Create a list of opposite channelIDs
    for (int p=0; p<16; p++) {
        pairs[p][0] = p + 64;
        pairs[p][1] = 31 - p + 64;
    }
    
    
    for (int i=0; i<n; i++) {
        T->GetEntry(i);
        if (i % 100 ==0) printf("\r%d/%d", i, n);
        
        for (int m=0; m<16; m++) {
            //Check if the events find both channels of the mth pair in a given event
            if (std::find(brChannelID->begin(), brChannelID->end(), pairs[m][0]) != brChannelID->end()) {
                if (std::find(brChannelID->begin(), brChannelID->end(), pairs[m][1]) !=brChannelID->end()) {
                    
                    //If both channels in the pair show up, use an iterator
                    //to retrieve their positions and times.
                    
                    int elem_first = pairs[m][0];
                    auto itr_first = std::find(brChannelID->begin(), brChannelID->end(), elem_first);
                    int index_first = std::distance(brChannelID->begin(), itr_first);
                    time_first = brTime->at(index_first);
                    int energy_first = brEnergy->at(index_first);
                                    
                    int elem_second = pairs[m][1];
                    auto itr_second = std::find(brChannelID->begin(), brChannelID->end(), elem_second);
                    int index_second = std::distance(brChannelID->begin(), itr_second);
                    time_second = brTime->at(index_second);
                    int energy_second = brEnergy->at(index_second);
                    //Add the difference in times of the hits to that channel pairs histogram
                    if (energy_first > energy_threshold && energy_second > energy_threshold) {
                    tdh[m]->Fill(time_first - time_second);
                    }
                }
            }
        }
    }
    
    TFile *f = new TFile(out_file, "RECREATE");
    
    for (int h=0; h<16; h++) {
        tdh[h]->Fit("gaus");
        tdh[h]->Write();
        delete tdh[h];
    } 

    f->Close();
    delete f;
}

int main(int argc, char* argv[])
{
    if (argc != 4) {
        printf("Three arguments required: an input filename and output filename, and energy threshold double for included hits.");
        exit(1);
    }

    all_sipm_pair_time_differences(argv[1], argv[2], atof(argv[3]));

    return 0;
}
