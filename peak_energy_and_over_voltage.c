#include <TH1D.h>
#include <TSpectrum.h>
#include <TFile.h>
#include <stdio.h>
#include <TTree.h>
#include <vector>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TF1.h>


std::vector<float>		*brEnergy;
std::vector<unsigned int>	*brChannelID;
float brStep1;
float brStep2;


int peak_energy_and_over_voltage(int channel_num, double vthresh, const char* in_file, const char* out_file) {    

    
    TFile *foo = new TFile(in_file);
    TTree *T;
    T = (TTree *) foo->Get("data");
    T->SetBranchAddress("channelID",&brChannelID);
    T->SetBranchAddress("energy",&brEnergy);
    T->SetBranchAddress("step1",&brStep1);	
    T->SetBranchAddress("step2",&brStep2);


    float real_vth = 10000.0*(vthresh+1)+100.0*(10+1)+10.0+1.0;
    float unique_ovs[128];
    int num_ovs = 1;
    int indices[128];
    float over_voltages[128];
    float vths[128];
    unsigned int n = T->GetEntries();
    T->GetEntry(0);
    float current_ov = brStep1;
    float current_vth = brStep2;
    int indices_count = 1;
    over_voltages[0] = current_ov;
    vths[0] = current_vth;
    unique_ovs[0] = current_ov;
    indices[0] = 0;
    //Loop through all events
    for (int i=0; i<n; i++) {
        T->GetEntry(i);
	if (i % 1000 == 0) {
	    printf("\r%d/%d", i, n);
	}
        //Check if the over voltage changes
        if (brStep1 != current_ov) {
            //If it does, save what event it happens at and the what the new
            //over voltage and threshold voltages are
    	    indices[indices_count] = i;
            over_voltages[indices_count] = brStep1;
            vths[indices_count] = brStep2;
            indices_count++;
            current_ov = brStep1;
            current_vth = brStep2;
            int unique = 0;
            for (int j=0; j<num_ovs; j++) {
                if (unique_ovs[j] == current_ov) {
                    unique++;
                }
            }
            if (unique == 0) {
                unique_ovs[num_ovs] = current_ov;
                num_ovs++;
            }
        }
    }


    float start_indices[num_ovs];
    float end_indices[num_ovs];
    int taken_indices = 0;
    for (int i=0; i<indices_count; i++) {
	printf("The %dth vth is %f\n", i, vths[i]);
	printf("The desired vth is %f\n", real_vth);
	if (vths[i]==real_vth) {
	    start_indices[taken_indices] = indices[i];
	    end_indices[taken_indices] = indices[i+1];
	    taken_indices++;
	}
    }
    printf("\n%d good indices found", taken_indices);
    printf("\n There are %d over voltage changes", indices_count);
    printf("\n There are %d unique over voltages", num_ovs);


    TH1D *tdh[num_ovs];
    char histoname[256];
    char histotitle[256];
    for (int i=0; i<num_ovs; i++) {
	sprintf(histoname, "h%d", i);
	sprintf(histotitle, "Energy for Channel %d, ov%f, vth%f", channel_num, unique_ovs[i], vthresh);
	tdh[i] = new TH1D(histoname, histotitle, 250, 0, 250);
    }


    for (int i=0; i<num_ovs; i++) {
	int start_index = start_indices[i];
	int end_index =end_indices[i];
	for (int j=start_index; j<end_index; j++) {
	    T->GetEntry(j);
	    if (std::find(brChannelID->begin(), brChannelID->end(), channel_num) != brChannelID->end()) {
		auto itr = std::find(brChannelID->begin(), brChannelID->end(), channel_num);
		int index = std::distance(brChannelID->begin(), itr);
		float energy = brEnergy->at(index);
		tdh[i]->Fill(energy);
	    }
	}
    }


    double x[num_ovs];
    double y[num_ovs];


    TApplication *theApp = new TApplication("App", 0 , 0);
    TCanvas *c1 = new TCanvas("c1", "Graph", 0, 0, 700, 500);
    /*
    for (int i=0; i<num_ovs; i++) {
	TSpectrum *pfinder = new TSpectrum();
	double sigma = 10;
	double threshold = 0.99;
	double num_peaks = pfinder->Search(tdh[i], sigma, "", threshold);
	Double_t *xpeaks = pfinder->GetPositionX();
	printf("The peak of the %dth over voltage is at %f", i+1, xpeaks[0]);
	y[i] = xpeaks[0];
	x[i] = unique_ovs[i];
    }
    */


    double ex[num_ovs];
    double ey[num_ovs];
    for (int i=0; i<num_ovs; i++) {
	TF1 *f1 = new TF1("f1", "gaus");
	tdh[i]->Fit("f1");
	x[i] = unique_ovs[i];
	y[i] = f1->GetParameter(1);
	ey[i] = f1->GetParError(1);
	ex[i] = 0;
	printf("The peak is centered at %f\n", y[i]);
	printf("The peak error is %f\n", ey[i]);
    }

    
    foo->Close();
    TFile *bar = new TFile(out_file, "RECREATE");


    TGraphErrors *gr = new TGraphErrors(num_ovs, x, y, ex, ey);
    char graph_title[128];
    sprintf(graph_title, "1.27 MeV peaks at various over voltages and vth = %2.0f", vthresh);
    gr->SetTitle(graph_title);
    gr->Draw("ALP");
    gr->GetXaxis()->SetTitle("Over Voltage (V)");
    gr->GetYaxis()->SetTitle("Energy (ADC)");
    gr->Write();

    c1->Modified();
    c1->Update();
    theApp->Run();


    return 0;
}

int main(int argc, char* argv[]) {
    if (argc != 5) {
	printf("There should be four arguments, a channel number, a threshold voltage, an input file and an output file");
	return 1;
    }
    peak_energy_and_over_voltage(atoi(argv[1]), atof(argv[2]), argv[3], argv[4]);
    return 0;
}
