#include <string>
#include <stdio.h>
#include <TTree.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TGraph2D.h>
#include <TGraph2DErrors.h>
#include <TMultiGraph.h>
#include <vector>
#include <TH1D.h>
#include <TF1.h>
#include <TF2.h>
#include <TSpectrum.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TMath.h>
#include <gsl/gsl_sf_psi.h>
#include <set>
#include <TLegend.h>

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


struct s {
    float value;
    int index;
};


int cmp(const void *a, const void *b) {
    struct s *a1 = (struct s *)a;
    struct s *a2 = (struct s *)b;
    if ((*a1).value > (*a2).value)
	return 1;
    else if ((*a1).value < (*a2).value) 
	return -1;
    else
	return 0;
}


Double_t poly_func(Double_t *val, Double_t *par) {
    Double_t x = val[0];
    Double_t y = val[1];
    double n = x*par[1];
    double k = y*par[2];
    if ((n < 0) || (k < 0) || (n < k)) {
	return 1e99;
    }
    Double_t f = pow(pow(par[3],2) + 43e6*par[0]*(gsl_sf_psi_n(1, n-k+1) - gsl_sf_psi_n(1, n + 1)), 0.5);
    return f;
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

    //Create application so that plot can draw from the terminal
    TApplication *theApp = new TApplication("App",0 ,0);
    TCanvas *c1 = new TCanvas("c1", "Graph", 0, 0, 700, 500);
    c1->Connect("Closed()", "TApplication", theApp, "()");

    //Initialize 2d data and error arrays
    double x[512];	//laser energy
    double y[512];	//vth
    double z[512];	//Resolution
    double ex[512];	
    double ey[512];
    double ez[512];


    double p0, p1, p2, p3, p0err, p1err, p2err, p3err;

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
		//printf("%dth start index is %d\n",good_indices ,start_indices[good_indices]);
		//printf("%dth end index is %d\n",good_indices ,end_indices[good_indices]);
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
	    printf("The %dth energy peak (vth=%f) is %f\n\n", k+1, unique_vths[k], x[k+num_points]);
	    printf("The median energy is %f\n", Median(epeaks[k]));


	    //Get the vth for each of these points
	    y[k + num_points] = unique_vths[k];
	    ey[k + num_points] = 0; 


	    //Get the standard deviations of the timing difference plots
	    TF1 *f2 = new TF1("f2", "gaus");
	    tdh[k]->Fit("f2");
	    z[k + num_points] = f2->GetParameter(2);
	    ez[k + num_points] = f2->GetParError(2);
	    printf("The %dth timing resolution (vth=%f) is %f\n\n", k+1, unique_vths[k], z[k+num_points]);
	}


	num_points += good_indices;
	foo->Close();
    }

    //Find the data extrema so appropriate bounds can be set for plotting
    xlow = x[0];
    xhigh = x[0];
    double ylow = y[0];
    double yhigh = y[0];
    double zlow = z[0];
    double zhigh = z[0];

    for (int i=1; i<num_points; i++) {
	if (x[i] < xlow) {
	    xlow = x[i];
	}	
	if (x[i] > xhigh) {
	    xhigh = x[i];
	}
	if (y[i] < ylow) {
	    ylow = y[i];
	}	
	if (y[i] > yhigh) {
	    yhigh = y[i];
	}
	if (z[i] < zlow) {
	    zlow = z[i];
	}	
	if (z[i] > zhigh) {
	    zhigh = z[i];
	}
    }

    //Make an output file for the TGraph2D
    TFile *foo = new TFile(outfile, "RECREATE");
    //Populate the graph with previously constructed arrays
    TGraph2DErrors *gr = new TGraph2DErrors(num_points, x, y, z, ex, ey, ez);
    TAxis *xaxis = gr->GetXaxis();
    TAxis *yaxis = gr->GetYaxis();
    TAxis *zaxis = gr->GetZaxis();

    xaxis->SetTitleOffset(2.5);
    yaxis->SetTitleOffset(2);
    zaxis->SetTitleOffset(2);   
    gr->SetMinimum(zlow);
    gr->SetMaximum(zhigh);

    //Edit plot labels
    gr->SetTitle("Laser Intensity and Threshold Voltage vs. Timing Resolution;Laser Energy (MeV);Threshold Voltage (V);Timing Resolution (Ps)");
    
    //Apply user defined fit
    //TF2 *f2d = new TF2("f2d", "sqrt([0]**2 + ([1]*(y-[3])/(x**[2]))**2)");
    TF2 *fp = new TF2("fp", poly_func, xlow, xhigh, ylow, yhigh, 4, 2);
    fp->SetMinimum(zlow);
    fp->SetMinimum(zhigh);
    
    fp->SetParameters(100,1000,0.1,60);
    fp->SetParNames("p0","p1","p2","p3");

    gr->Fit("fp", "RV");
    gr->Draw("surf2");
    /*f2d->SetMinimum(zlow);
    f2d->SetMaximum(zhigh); 
    f2d->SetRange(xlow, ylow, xhigh, yhigh);
    gr->Fit("f2d", "U", "R"); //U option uses user defined fit
    p0 = f2d->GetParameter(0);
    p1 = f2d->GetParameter(1);
    p2 = f2d->GetParameter(2);
    p3 = f2d->GetParameter(3);
    p0err = f2d->GetParError(0);
    p1err = f2d->GetParError(1);
    p2err = f2d->GetParError(2);
    p3err = f2d->GetParError(3);
    */

    double p0poly = fp->GetParameter(0);
    double p1poly = fp->GetParameter(1);
    double p2poly = fp->GetParameter(2);
    double p3poly = fp->GetParameter(3);
    double p0polyerr = fp->GetParError(0);
    double p1polyerr = fp->GetParError(1);
    double p2polyerr = fp->GetParError(2);
    double p3polyerr = fp->GetParError(3);
    fp->Draw("same");

    //f2d->Draw("SAME");
    gr->Write();
    c1->Modified();
    c1->Update();
    //theApp->Run();
    //theApp->();


    float unique_vths[num_vths];
    int vcount = 0;
    for (auto elem: global_vths) { //Iterate through the global_vths set (which is ordered the same way that data is acquired) and input values into an array which the TMultiGraph can use
	unique_vths[vcount] = elem;
	vcount++;
    }




    struct s smgx[num_points];
    for (int i=0; i<num_points; i++) {
        smgx[i].value = x[i];
        smgx[i].index = i;
    }
    qsort(smgx, num_points, sizeof(smgx[0]), cmp);


    double sorted_mgx[num_points];
    double sorted_mgy[num_points];
    double sorted_mgex[num_points];
    double sorted_mgey[num_points];
    double sorted_vth[num_points];
    for (int i=0; i<num_points; i++) {
	int t = smgx[i].index;
        sorted_mgx[i] = x[t]; 
        sorted_mgy[i] = z[t];
        sorted_mgex[i] = ex[t];
        sorted_mgey[i] = ez[t];
        sorted_vth[i] = y[t];
    }
    
    //Make the multigraph with energy vs resolution at each vth
    TMultiGraph *mg = new TMultiGraph();
    for (int i=0; i<num_vths; i++) {
	int points = num_points/nfiles;
	double mgx[points];
	double mgy[points];
	double mgex[points];
	double mgey[points];
	int points_added = 0;
	
	char graph_title[128];
	char data_title[128] = "Data";
	sprintf(graph_title, "Data vs. Fit for DAC Threshold %2.0f", unique_vths[i]);
	
	float measured_pointsx[nfiles];
	float measured_errx[nfiles];
	float measured_pointsy[nfiles];
	float measured_erry[nfiles];
	float predicted_pointsx[1000];
	float predicted_errx[1000];
	float predicted_pointsy[1000];
	float predicted_erry[1000];
	for (int k=0; k<num_points; k++) {
	    if (sorted_vth[k] == unique_vths[i]) {//Check if we have a correct vth
		mgx[points_added] = sorted_mgx[k];
		mgex[points_added] = sorted_mgex[k];
		mgy[points_added] = sorted_mgy[k];
		mgey[points_added] = sorted_mgey[k];

		measured_pointsx[points_added] = mgx[points_added];
		measured_errx[points_added] = mgex[points_added];
	        measured_pointsy[points_added] = mgy[points_added];
		measured_erry[points_added] = mgey[points_added];

		points_added++;

        	for (int j=0; j<1000; j++) {
		    double xpos = xlow + (xhigh-xlow)*j/1000;
		    double x = xpos;
    	            double y = unique_vths[i];
    		    double n = x*p1poly;
    		    double k = y*p2poly;
		    double fit_point = pow(pow(p3poly,2) + pow(p0poly*gsl_sf_psi_n(1, n-k+1) - gsl_sf_psi_n(1, n+1),2), 0.5);
	    	    predicted_pointsx[j] = xpos;
		    predicted_errx[j] = 0.00001;
		    predicted_pointsy[j] = fit_point;
		    predicted_erry[j] = 0.00001;
		}
	    }
	}
	TGraphErrors *measured = new TGraphErrors(nfiles, measured_pointsx, measured_pointsy, measured_errx, measured_erry);

	TGraphErrors *predicted = new TGraphErrors(1000, predicted_pointsx, predicted_pointsy, predicted_errx, predicted_erry);
	
	
	measured->SetName("measured");
	measured->SetTitle(data_title);
	measured->SetDrawOption("AP");
	measured->SetMarkerStyle(20);
	predicted->SetName("predicted");
	predicted->SetTitle("Fit");
	predicted->SetLineColor(0);
	predicted->SetMarkerStyle(0);
	TMultiGraph *g = new TMultiGraph();
	g->Add(measured);
	g->Add(predicted);
	g->SetTitle(graph_title);
	g->GetXaxis()->SetTitle("Laser Energy (MeV)");
	g->GetYaxis()->SetTitle("Time Difference Resolution");
	g->Write();
	TLegend *leg = new TLegend(0.1,0.7, 0.3, 0.9);
	leg->SetNColumns(2);
	leg->SetHeader("Test");
	leg->SetFillColor(0);
	leg->AddEntry("measured", "Data", "l");
	leg->AddEntry("predicted", "Fit", "l");
	leg->SetName("My_leg");
	leg->Write();


	TGraphErrors *mgr = new TGraphErrors(points_added, mgx, mgy, mgex, mgey);   //In these x is energy and y is time resolution
	char mgraphname[256];
	sprintf(mgraphname, "Threshold Voltage %f", unique_vths[i]);
	mgr->SetName(mgraphname);
	mgr->SetTitle(mgraphname);
	int color = floor(100*i/num_vths) + 1;
	mgr->SetLineColor(color);
	gr->SetMarkerStyle(20+i);
	mg->Add(mgr);
	mg->Draw();	
    }

    printf("The parameters are p0:%f, p1:%f, p2:%f, p3:%f", p0poly ,p1poly, p2poly, p3poly);
    mg->GetXaxis()->SetTitle("Laser Energy (MeV)");
    mg->GetYaxis()->SetTitle("Time Difference Resolution (Ps)");
    mg->Write();
    c1->BuildLegend();
    foo->Close();
    

    return 0;
}
