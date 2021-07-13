#include <TFile.h>
#include <TTree.h>
#include <stdio.h>
#include <TH2D.h>
#include <vector>

float		brStep1;
float		brStep2;
unsigned short  brPrevEventFlags;
long long       brPrevEventTime;
double          brTimeLast;
double          brTimeLastTmp;
long long 	brStepBegin;
long long 	brStepEnd;

/* Should be the number of hits in the event. But the number here disagrees
 * with the size of the vectors? Maybe there is a limit. */
unsigned short	brN;
int             brChannelIdx[128];
/* Time is in ps. */
std::vector<long long>        	*brTime;
std::vector<unsigned int>	*brChannelID;
std::vector<float>		*brToT;
std::vector<unsigned short>	*brT1Coarse;
std::vector<unsigned short>	*brT1Fine;
std::vector<unsigned short>	*brT2Coarse;
std::vector<unsigned short>	*brT2Fine;
std::vector<unsigned short>	*brQCoarse;
std::vector<unsigned short>	*brQFine;
/* In arbitrary ADC units. */
std::vector<float>	        *brEnergy;
std::vector<float>	        *brQT1;
std::vector<float>	        *brQT2;
std::vector<unsigned short>	*brTacID;

//Making the empty histogram
TH2D h("h","Channels with coinciding events",
        32,64,96,   // X axis
        32,64,96);  // Y axis
        
int twod_channel_coincidence(const char* filename) {

    TFile foo(filename);
    TTree *T = (TTree*) foo.Get("data");
    
    //Get the addresses for the time and channel number
    long long time;   
    T->GetEntry(0);

    unsigned int channelID;

    T->SetBranchAddress("step1", &brStep1);
    T->SetBranchAddress("step2", &brStep2);
    T->SetBranchAddress("prevEventFlags", &brPrevEventFlags);
    T->SetBranchAddress("prevEventTime", &brPrevEventTime);
    T->SetBranchAddress("timeLast", &brTimeLast);
    T->SetBranchAddress("mh_n", &brN);
    T->SetBranchAddress("channelIdx", brChannelIdx);
    T->SetBranchAddress("tot", &brToT);
    T->SetBranchAddress("t1coarse", &brT1Coarse);
    T->SetBranchAddress("t1fine", &brT1Fine);
    T->SetBranchAddress("t2coarse", &brT2Coarse);
    T->SetBranchAddress("t2fine", &brT2Fine);
    T->SetBranchAddress("qcoarse", &brQCoarse);
    T->SetBranchAddress("qfine", &brQFine);
    T->SetBranchAddress("time", &brTime);
    T->SetBranchAddress("channelID", &brChannelID);
    T->SetBranchAddress("energy", &brEnergy);
    T->SetBranchAddress("qT1", &brQT1);
    T->SetBranchAddress("qT2", &brQT2);
    T->SetBranchAddress("tacID", &brTacID);
    
    long long prev_time = time;
    
    unsigned int n = T->GetEntries();
    float energy;
    /*
    For consecutive evnts that are within the desired time window,
    increment the length of the event. If the next hit is outside the
    time window, end the event and add each permuutation of channelID
    pairs to the histogram data.
    */
    
    int hits[128];
    
    for (int i=1; i<n; i++) {
        T->GetEntry(i);
        if (i % 1000 == 0) {printf("\r%d/%d", i, n); }
        for (int j=0; j<brChannelID->size(); j++) {
            hits[j] = brChannelID->at(j);
        }
        for (int k=0; k<brChannelID->size(); k++) {
            for (int m=0; m<k; m++) {
                h.Fill(hits[m],hits[k]);
                h.Fill(hits[k],hits[m]);
            }
        }
    }
    //Take the maximum bin as two opposite sipms and print out its location.
    Int_t MaxBin = h.GetMaximumBin();
    Int_t x,y,z;
    h.GetBinXYZ(MaxBin, x, y, z);
    printf("The bins with maximum coincidence are (%d, %d)", x+63, y+63);
    h.Draw("colz");
    
    return 0;    
}

int main(int argc, const char* argv[]) {
    if (argc != 2) {
        printf("One argument required: input filename");
    }
    twod_channel_coincidence(argv[1]);
    return 0;
}
