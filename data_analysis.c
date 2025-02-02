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
#include <utility>
#include <stdexcept>
#include <sstream>
#include <fstream>
#include <iostream>


using namespace std;


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
    double n = x*par[3];
    double k = (y-par[2])*par[1];
    if ((n < 0) || (k < 0) || (n < k)) {
	//printf("n is %f and k is %f\n", n, k);
	return 1e9;
    }
    Double_t f = pow(pow(par[4],2) + par[0]*(gsl_sf_psi_n(1, n-k+1) - gsl_sf_psi_n(1, n + 1)), 0.5);
    //printf("f is %f, n is %f and k is %f\n",f, n, k);
    //printf("gsl terms are: %f\n", gsl_sf_psi_n(1, n-k+1) - gsl_sf_psi_n(1, n + 1));
    return f;
}

static void usage(void) {
    fprintf(stderr,
"data_analysis"
"\n"
"  --infile <filename>   Input filename.\n"
"  --outfile <filename>  Output filename.\n"
"  --help                Output this help and exit.\n"
"\n");
    exit(1);
}


int main(int argc, char* argv[])
{
    char *infile = NULL;
    char *outfile = NULL;
    int nfiles = 0;
    for (int i=1; i<argc; i++) {
	if (strcmp(argv[i], "--infile") == 0) {
	    i++;
            infile = argv[i];
	    printf("Extracting data from %s\n", infile);
	} else if (strcmp(argv[i], "--outfile") == 0) {
	   i++;
	   outfile = argv[i];
	   printf("Outputting to %s\n", outfile);
	} else if (strcmp(argv[i], "--help") == 0) {
            usage();
	}
    }

    if (infile == NULL || outfile == NULL)
        usage();

    //Initialize 2d data and error arrays
    double xtemp[512];	//laser energy
    double ytemp[512];	//vth
    double ztemp[512];	//Resolution
    double extemp[512];	
    double eytemp[512];
    double eztemp[512];

    FILE *myFile = fopen(infile,"r");

    if (!myFile) {
        fprintf(stderr, "error opening file '%s'\n", infile);
        exit(1);
    }

    char line[1024];
    int counter = 0; 
    int lines = 0;
    while (fgets(line,1024,myFile)) {
        if (lines == 0) {
            lines += 1;
            continue;
        }

        char *tok = strtok(line, ",");
        counter = 0;
        while (tok) { 
            double dat = atof(tok);

            switch (counter % 7) {
            case 0:
                xtemp[lines - 1] = dat;
                break;
            case 1:
                ytemp[lines - 1] = dat;
                break;
            case 2:
                ztemp[lines - 1] = dat;
                break;
            case 3:
                extemp[lines - 1] = dat;
                break;
            case 4:
                eytemp[lines - 1] = dat;
                break;
            case 5:
                eztemp[lines - 1] = dat;
                break;
            }

            counter++;
            tok = strtok(NULL, ",");
        }

        lines += 1;
    }

    fclose(myFile);

    int num_points = lines - 1;

    double x[num_points];
    double y[num_points];
    double z[num_points];
    double ex[num_points];
    double ey[num_points];
    double ez[num_points]; 


    std::set<float> global_vths = {};


    for (int i=0; i<num_points; i++) {
       	x[i] = xtemp[i];
	y[i] = ytemp[i];
	z[i] = ztemp[i];
	ex[i] = extemp[i];
	ey[i] = eytemp[i];
	ez[i] = eztemp[i];
	global_vths.insert(y[i]);
    }

    if (num_points % global_vths.size()) {
        fprintf(stderr, "number of points %i is not a multiple of the number of thresholds %zu\n", num_points, global_vths.size());
        fprintf(stderr, "unable to determine the number of files\n");
        exit(1);
    }

    nfiles = num_points/global_vths.size();

    float unique_vths[128];
    int vcount = 0;
    for (auto elem: global_vths) {
	unique_vths[vcount] = elem;
	vcount++;
    }

     
    double p0, p1, p2, p3, p0err, p1err, p2err, p3err;

    //Find the data extrema so appropriate bounds can be set for plotting
    double xlow = x[0];
    double xhigh = x[0];
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
    double fitx[128];
    double fity[128];
    double fitz[128];
    double fitex[128];
    double fitey[128];
    double fitez[128];
    int fit_points = 0;
    for (int i=0; i<num_points; i++) {
	if (y[i] != 2) {
	    fitx[fit_points] = x[i];
	    fity[fit_points] = y[i];
	    fitz[fit_points] = z[i];
	    fitex[fit_points] = ex[i];
	    fitey[fit_points] = ey[i];
	    fitez[fit_points] = ez[i];
	    fit_points++;
	}
    }
    printf("fit points = %i\n", fit_points);
    TGraph2DErrors *gr = new TGraph2DErrors(fit_points, fitx, fity, fitz, fitex, fitey, fitez);
    TAxis *xaxis = gr->GetXaxis();
    TAxis *yaxis = gr->GetYaxis();
    TAxis *zaxis = gr->GetZaxis();

    xaxis->SetTitleOffset(2.5);
    yaxis->SetTitleOffset(2);
    zaxis->SetTitleOffset(2);   
    gr->SetMinimum(zlow);
    gr->SetMaximum(zhigh);

    //Edit plot labels
    gr->SetTitle("Laser Intensity and Threshold Voltage vs. Timing Resolution;Effective Energy (MeV);Threshold Voltage (V);Timing Resolution (Ps)");
    
    //Apply user defined fit
    TF2 *fp = new TF2("fp", poly_func, xlow, xhigh, ylow, yhigh, 5);
    fp->SetMinimum(zlow);
    fp->SetMinimum(zhigh);
    
    //fp->SetParameters(1, 10, 5, 1500, 75);
    fp->SetParLimits(1, 0.1, 100);
    fp->SetParLimits(2, 0, 10);
    fp->SetParLimits(3, 1000, 2000);
    fp->SetParLimits(4,40,100);
    fp->SetParNames("p0","p1", "p2", "p3", "p4");
    


    gr->Fit("fp", "RV");
    gr->Draw("surf2");


    double p0poly = fp->GetParameter(0);
    double p1poly = fp->GetParameter(1);
    double p2poly = fp->GetParameter(2);
    double p3poly = fp->GetParameter(3);
    double p4poly = fp->GetParameter(4);
    double p0polyerr = fp->GetParError(0);
    double p1polyerr = fp->GetParError(1);
    double p2polyerr = fp->GetParError(2);
    double p3polyerr = fp->GetParError(3);
    double p4polyerr = fp->GetParError(4);
    fp->Draw("same");

    //f2d->Draw("SAME");
    gr->Write();
    //c1->Modified();
    //c1->Update();
    //theApp->Run();
    //theApp->Terminate();

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
    
    TLegend *l = new TLegend(0.7, 0.5, 0.87, 0.87);
    l->SetFillColor(0);

    //Make the multigraph with energy vs resolution at each vth
    TMultiGraph *mg = new TMultiGraph();
    for (int i=0; i<vcount; i++) {
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
		    double fit_point = pow(pow(p3poly,2) + p0poly*(gsl_sf_psi_n(1, n-k+1) - gsl_sf_psi_n(1, n+1)), 0.5);
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
	measured->SetLineColor(30);
	measured->SetMarkerStyle(20);
	predicted->SetName("predicted");
	predicted->SetTitle("Fit");
	predicted->SetLineColor(0);
	predicted->SetMarkerStyle(0);
	TMultiGraph *g = new TMultiGraph();
	g->Add(measured);
	g->Add(predicted);
	g->SetTitle(graph_title);
        /* Can't set the axis titles directly due to
         * https://root-forum.cern.ch/t/tmultigraph-getxaxis-settitle-version-4-01/3928/5,
         * so instead we set them via the SetTitle() command. */
	g->SetTitle(";Effective Energy (MeV);Time Difference Resolution (ps)");
	//g->GetXaxis()->SetTitle("Laser Energy (MeV)");
	//g->GetYaxis()->SetTitle("Time Difference Resolution");
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
	sprintf(mgraphname, "DAC Threshold %.0f", unique_vths[i]);
	mgr->SetName(mgraphname);
	mgr->SetTitle(mgraphname);
	int color = floor(100*i/vcount) + 1;
	mgr->SetLineColor(color);
	gr->SetMarkerStyle(20+i);
        if ((int) unique_vths[i] == 4 || (int) unique_vths[i] == 10 || (int) unique_vths[i] == 20 || (int) unique_vths[i] == 30 || (int) unique_vths[i] == 40 || (int) unique_vths[i] == 50) {
            mg->Add(mgr);
            l->AddEntry(mgr, mgraphname, "lp");
        }
    }

    printf("The parameters are p0:%f, p1:%f, p2:%f, p3:%f, p4:%f", p0poly ,p1poly, p2poly, p3poly, p4poly);
    mg->SetTitle(";Effective Energy (MeV);Time Difference Resolution (ps)");
    //mg->GetXaxis()->SetTitle("Laser Energy (MeV)");
    //mg->GetYaxis()->SetTitle("Time Difference Resolution (Ps)");
    TCanvas *c = new TCanvas("c1");
    mg->Draw("apl");
    l->Draw();
    mg->Write();
    l->Write();
    c->SaveAs("hist.pdf","pdf");

    printf("\n\nThere are %d vths\n", vcount);
    //c1->BuildLegend();

    delete fp;
    delete gr;
    foo->Close();
    delete foo;
    

    return 0;
}
