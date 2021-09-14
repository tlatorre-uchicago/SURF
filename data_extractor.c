#include <string>
#include <stdio.h>
#include <TTree.h>
#include <TFile.h>
#include <vector>
#include <TH1D.h>
#include <TF1.h>
#include <TF2.h>
#include <TSpectrum.h>
#include <TMath.h>
#include <set>
#include <fstream>


using namespace std;


std::vector<float>		*brEnergy;
std::vector<unsigned int>	*brChannelID;
std::vector<long long >		*brTime;
float brStep1;
float brStep2;

//Function to find the median of a histogram
//This will be used as a starting point for the mean parameter
//of the Gaussian fits at the end of main
double Median(const TH1D * h1) {
    int n = h1->GetXaxis()->GetNbins();
    std::vector<double> x(n);
    h1->GetXaxis()->GetCenter( &x[0] );
    const double * y = h1->GetArray();
    return TMath::Median(n, &x[0], &y[1]);
}


int main(int argc, char* argv[]) {
    int i = 0;
    double ov;
    int bar;
    char *filenames[256];
    double tunes[256];
    char outfile[256];
    double xlow;
    double xhigh;

    //Parse the input to find the number of files, selected over voltage and bar number
    int nfiles = 0;
    for (i=1; i<argc; i++) {
	if (strcmp(argv[i], "--ov") == 0) {
	    i++;
	    ov = atof(argv[i]);
	    printf("Checking for over_voltage = %f \n", ov);
	}
	else if (strcmp(argv[i], "--bar") == 0) {
	    i++;
	    bar = atoi(argv[i]);
	    printf("checking bar %d\n", bar);
	}
	else if (strcmp(argv[i], "--outfile") == 0) {
	    i++;
	    strcpy(outfile, argv[i]);
	}
	else {
	    filenames[nfiles] = argv[i];
	    std::string str = filenames[nfiles];
	    std::string str2 = ("tune_");
	    std::size_t found = str.find(str2);
	    int end_length = str.length();
	    str.replace(found + 7, end_length, "");
	    str.replace(0, found + 5, "");
	    double tune = stof(str);
	    tunes[nfiles] = tune;
	    nfiles++;
	}
    }

    //Get the different channels from the input bar
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


    //Initialize 2d data and error arrays
    double x[512];	//laser energy
    double y[512];	//vth
    double z[512];	//Resolution
    double ex[512];	
    double ey[512];
    double ez[512];


    //Store these up here because they will need to be accessed for the multigraph
    int num_points = 0;
    int num_vths = 0;
    std::set<float> global_vths = {};


    for (int j=0; j<nfiles; j++) {
	printf("Checking file %d/%d\n", j+1, nfiles);
	TFile *foo = new TFile(filenames[j]);
	TTree *T;
	T = (TTree *) foo->Get("data");
	T->SetBranchAddress("channelID",&brChannelID);
        T->SetBranchAddress("energy",&brEnergy);
        T->SetBranchAddress("step1",&brStep1);
        T->SetBranchAddress("step2",&brStep2);
	T->SetBranchAddress("time",&brTime);


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
        for (int k=0; k<n; k++) {
	    if (k % 1000 == 0) {
		printf("\r%d/%d", k, n);
	    }
	    T->GetEntry(k);
	    //Check if the over voltage changes
	    if (brStep1 != current_ov) {
		//If it does, save what event it happens at and the what the new
		//over voltage and threshold voltages are
		indices[indices_count] = k;
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


	num_vths = indices_count/num_ovs;
	float unique_vths[num_vths];
	for (int k=0; k<num_vths; k++) {
	    unique_vths[k] = vths[k*num_ovs];
	}

	//Convert vths to readable form
	for (int k=0; k<num_vths; k++) {
	    int diff = floor(unique_vths[k]);
	    unique_vths[k] = (unique_vths[k] - (diff % 10000))/10000 - 1;
	    global_vths.insert(unique_vths[k]);
	}

	//Find the event indices that contain the events with chosen ov
	int start_indices[num_vths];
	int end_indices[num_vths];
	int good_indices = 0;
	for (int k=0; k<indices_count; k++) {
	    if (over_voltages[k] == ov) {
		start_indices[good_indices] = indices[k];
		if (k == indices_count-1) {
		    end_indices[good_indices] = n;
		}
		else {
		    end_indices[good_indices] = indices[k+1];

		}
		good_indices++;
	    }
	}
	
	//Construct timing difference histograms
	TH1D *tdh[good_indices];
	char histoname[256];
	char histotitle[256];
	for (int k=0; k<good_indices; k++) {
	    sprintf(histoname, "tdh_vth_%f_tune_%f", unique_vths[k], tunes[j]);
	    sprintf(histotitle, "Time differences for vth %f and tune %f", unique_vths[k], tunes[j]);
	    tdh[k] = new TH1D(histoname, histotitle, 500, -1500, 1500);
	}

	//Construct energy histograms
	TH1D *epeaks[good_indices];
	for (int k=0; k<good_indices; k++) {
	    sprintf(histoname, "energy_vth_%f_tune_%f", unique_vths[k], tunes[j]);
	    sprintf(histotitle, "Energy distribution for vth %f and tune %f", unique_vths[k], tunes[j]);
	    epeaks[k] = new TH1D(histoname, histotitle, 300, 0, 1500);
	}


	printf("\nFilling histograms\n");

	int first_time;
        int second_time;
	float first_energy;
	float second_energy;
	int start_index;
	int end_index;
	printf("There are %d good indices\n", good_indices);
	for (int k=0; k<good_indices; k++) {
	    start_index = start_indices[k];
	    end_index = end_indices[k];
	    printf("\rHistogram %d/%d", k, good_indices);
	    for (int m=start_index; m<end_index; m++) {
		T->GetEntry(m);
		if (std::find(brChannelID->begin(), brChannelID->end(), first_channel) != brChannelID->end()) {
		  
		    auto itr = std::find(brChannelID->begin(), brChannelID->end(), first_channel);
		    int energy_index = std::distance(brChannelID->begin(), itr);
		    float histo_energy = brEnergy->at(energy_index);
		    if (histo_energy >= 10) {
		        epeaks[k]->Fill(histo_energy, 1);
		    }
		    
			
		    if (std::find(brChannelID->begin(), brChannelID->end(), second_channel) != brChannelID->end()) {
			auto first_itr = std::find(brChannelID->begin(), brChannelID->end(), first_channel);
			int first_index = std::distance(brChannelID->begin(), first_itr);
			first_energy = brEnergy->at(first_index);
			first_time = brTime->at(first_index);


			auto second_itr = std::find(brChannelID->begin(), brChannelID->end(), second_channel);
			int second_index = std::distance(brChannelID->begin(), second_itr);
			second_energy = brEnergy->at(second_index);
			second_time = brTime->at(second_index);


			tdh[k]->Fill(first_time - second_time);
			//epeaks[k]->Fill(second_energy);
		    }
		}
	    }
	}


	for (int k=0; k<good_indices; k++) {
	    //Get the means of the energy plots
	    int low = floor(Median(epeaks[k]) - 20);
	    int high = floor(Median(epeaks[k]) + 21);
	    TF1 *f1 = new TF1("f1", "gaus");
	    f1->SetRange(low, high);
	    f1->SetParameter(1, Median(epeaks[k]));
	    epeaks[k]->Fit("f1", "R");
	    x[k + num_points] = f1->GetParameter(1) / 74.74 * 1.27;
	    ex[k + num_points]  = f1->GetParError(1) / 74.74 * 1.27;


	    //Get the vth for each of these points
	    y[k + num_points] = unique_vths[k];
	    ey[k + num_points] = 0; 


	    //Get the standard deviations of the timing difference plots
	    TF1 *f2 = new TF1("f2", "gaus");
	    tdh[k]->Fit("f2");
	    z[k + num_points] = f2->GetParameter(2);
	    ez[k + num_points] = f2->GetParError(2);
	}
	num_points += good_indices;
	foo->Close();
    }
    
    
    ofstream output;
    output.open(outfile);
    output << "Energy (MeV)" << "," << "VTH1 (ADC units)" << "," << "Time Width Resolution (Ps)" << "," << "Energy Uncertainty" << "," << "VTH1 Uncertainty" << "," << "Time Width Uncertainty (Ps)" << "," << "Over Voltage (V)" << endl;
    for (int i=0; i<num_points; i++) {
        output << x[i] << "," << y[i] <<"," << z[i] << "," << ex[i] << "," << ey[i] << "," << ez[i] << "," << ov << endl;
    }
    
    return 0;
}
