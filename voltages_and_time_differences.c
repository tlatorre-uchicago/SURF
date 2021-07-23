#include <vector>
#include <stdio.h>
#include <TH1D.h>
#include <TTree.h>
#include <TList.h>
#include <TFile.h>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <TH1D.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TApplication.h>
#include <TMultiGraph.h>
#include <TGraphErrors.h>
#include <TF1.h>


std::vector<float>		*brEnergy;
std::vector<unsigned int>	*brChannelID;
std::vector<long long>		*brTime;
float brStep1;
float brStep2;


int voltages_and_time_differences(const char* in_file, const char* out_file, int bar, float energy_threshold) {
    //Get the channel numbers from the input bar
    int first_channel;
    int second_channel;
    if (bar <= 7) {
	first_channel = 64 + 14 - 2*bar;
    }	    
    else if (bar > 7 && bar <= 15) {
	second_channel = 64 + 30 - 2*(bar - 8); 
    }
    else {
	printf("The bar should be an integer in [0,1,2,...,14,15]\n");
	return 1;
    }
    first_channel = 159 - second_channel;
    second_channel = 159 - first_channel;    
    printf("First channel is %d and second channel is %d\n", first_channel, second_channel);


    TFile *foo = new TFile(in_file);
    TTree *T;
    T = (TTree*) foo->Get("data");
    
    
    T->SetBranchAddress("step1", &brStep1);
    T->SetBranchAddress("step2", &brStep2);
    T->SetBranchAddress("energy", &brEnergy);
    T->SetBranchAddress("time", &brTime);
    T->SetBranchAddress("channelID", &brChannelID);


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
    int num_vths = indices_count/num_ovs;
    float unique_vths[num_vths];
    for (int i=0; i<num_vths; i++) {
	unique_vths[i] = vths[i*num_ovs];
    }


    for (int i=0; i<indices_count; i++) {
	int diff = floor(vths[i]);
	vths[i] = (vths[i] - (diff % 100000))/(10000);
    }
    for (int i=0; i<num_vths; i++) {
	int diff = floor(unique_vths[i]);
	unique_vths[i] = (unique_vths[i] - (diff % 100000))/(10000);
    }


    //Make a histogram for each pair of ov and threshold 
    TH1D *tdh[indices_count];
    char histoname[256];
    char histotitle[256];
    for (int i=0; i<indices_count; i++) {
	sprintf(histoname, "ov_%2.0f_vth_%2.0f", over_voltages[i], vths[i]);
	sprintf(histotitle, "Time differences for over voltage %2.0f and threshold voltage %2.0f", over_voltages[i], vths[i]);
	tdh[i] = new TH1D(histoname, histotitle, 300, -2000, 2000);	
    }


    int first_time;
    int second_time;
    float first_energy;
    float second_energy;
    int start_index;
    int end_index;
    for (int i=0; i<indices_count; i++) {
	if (i<indices_count - 1) {
	    start_index = indices[i];
	    end_index = indices[i+1];
	}
	if (i == indices_count - 1) {
	    start_index = indices[indices_count - 1];
	    end_index = n;
	}
	printf("\nChecking stage %d/%d\n", i + 1, indices_count);
	for (int j=start_index; j<end_index; j++) {
	    T->GetEntry(j);
	    if (j%100 == 0) {
                printf("\r%d/%d", j - start_index, end_index - start_index);
	    }
	    //Now check if an event has a time_difference
	    //Does the event have a hit in the first channel
	    if (std::find(brChannelID->begin(), brChannelID->end(), first_channel) != brChannelID->end()) {
		//Does the event also have a hit in the second channel
		if (std::find(brChannelID->begin(), brChannelID->end(), second_channel) != brChannelID->end()) {
		    //Find the time and energy from the hit at the first channel
		    auto first_itr = std::find(brChannelID->begin(), brChannelID->end(), first_channel);
		    int first_index = std::distance(brChannelID->begin(), first_itr);
		    first_energy = brEnergy->at(first_index);
		    first_time = brTime->at(first_index);



		    //Find the time and energy from the hit at the second channel
		    auto second_itr = std::find(brChannelID->begin(), brChannelID->end(), second_channel);
		    int second_index = std::distance(brChannelID->begin(), second_itr);
		    second_energy = brEnergy->at(second_index);
		    second_time = brTime->at(second_index);


		    //Check if both of the energies exceed the energy threshold
		    //If they do, add them to the correct histogram
		    if (first_energy > energy_threshold && second_energy > energy_threshold) {
    			tdh[i]->Fill(first_time - second_time);
		    }
		}
	    }
	}
	//Update the range of events to those with the next over voltage
    }

    float errors[128];
    float time_resolutions[128];
    TFile *f = new TFile(out_file, "RECREATE");
    for (int i=0; i<indices_count; i++) {
	TF1 *f1 = new TF1("f1", "gaus");
	tdh[i]->Fit("f1");
	time_resolutions[i] = f1->GetParameter(2);
	errors[i] = f1->GetParError(2);
	tdh[i]->Write();
	printf("The time resolution is %f.\n", time_resolutions[i]);
    }
 
    
    f->Close();


    TApplication *theApp = new TApplication("App",0,0);
    TCanvas* c1 = new TCanvas("c1", "MultiGraph", 0, 0, 700, 500);
    c1->Connect("Closed()", "TApplication", theApp, "Terminate()");
    TMultiGraph *mg = new TMultiGraph();
    for (int i=0; i<num_vths; i++) {
	double x[num_ovs];
	double y[num_ovs];
	double ex[num_ovs];
	double ey[num_ovs];
	int start_at = i*num_ovs;
	int end_at = i*num_ovs + num_ovs;
	for (int j=0; j<num_ovs; j++) {
	    x[j] = over_voltages[i*num_ovs+j];
            y[j] = time_resolutions[i*num_ovs+j];
	    ex[j] = 0.000001;
	    ey[j] = errors[i*num_ovs+j];
	}
	TGraphErrors *gr = new TGraphErrors(num_ovs, x, y, ex, ey);
	char graphname[256];
	sprintf(graphname, "Threshold voltage: %2.0f", unique_vths[i]);
	gr->SetName(graphname);
	gr->SetTitle(graphname);
	int color = floor(120*i/num_vths) + 1;
	gr->SetLineColor(color);
	printf("Setting line color to %d\n", color);
	gr->SetMarkerStyle(20+i);
	mg->Add(gr);
    }
    mg->Draw("ALP");
    mg->GetXaxis()->SetTitle("Over Voltage (V)");
    mg->GetYaxis()->SetTitle("Time Resolution (Ps)");


    c1->Modified();
    c1->Update();
    c1->BuildLegend();
    theApp->Run();
    theApp->Terminate();


    return 0;
}


int main(int argc, char* argv[]) {
    if (argc != 5 ) {
	printf("Four arguments required: input file, output file, bar number and energy threshold.");
	return 1;
    }
    voltages_and_time_differences(argv[1], argv[2], atoi(argv[3]), atof(argv[4]));
    return 0;
}
