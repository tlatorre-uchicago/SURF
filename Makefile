#all: energy_cut2 all_sipm_pair_time_differences twod_channel_coincidence branch_cutter voltages_and_time_differences peak_energy_and_over_voltage energy_threshold_resolution_fit polygamma_fit
all: energy_cut2 all_sipm_pair_time_differences twod_channel_coincidence branch_cutter voltages_and_time_differences peak_energy_and_over_voltage polygamma_fit data_analysis data_extractor

OPTIMIZATION?=-O0
CFLAGS=$(OPTIMIZATION) -Wall -g $(shell root-config --cflags)

LDLIBS=$(shell root-config --libs) $(shell gsl-config --libs)

energy_cut2: energy_cut2.c
	g++ -o $@ $^ $(CFLAGS) $(LDLIBS)

all_sipm_pair_time_differences: all_sipm_pair_time_differences.c
	g++ -o $@ $^ $(CFLAGS) $(LDLIBS)

twod_channel_coincidence: twod_channel_coincidence.c
	g++ -o $@ $^ $(CFLAGS) $(LDLIBS)

branch_cutter: branch_cutter.c
	g++ -o $@ $^ $(CFLAGS) $(LDLIBS)

voltages_and_time_differences: voltages_and_time_differences.c
	g++ -o $@ $^ $(CFLAGS) $(LDLIBS)

peak_energy_and_over_voltage: peak_energy_and_over_voltage.c
	g++ -o $@ $^ $(CFLAGS) $(LDLIBS)

#energy_threshold_resolution_fit: energy_threshold_resolution_fit.c
#	g++ -o energy_threshold_resolution_fit energy_threshold_resolution_fit.c  `root-config --cflags --libs` -lSpectrum

polygamma_fit: polygamma_fit.c
	g++ -o $@ $^ $(CFLAGS) $(LDLIBS)

data_analysis: data_analysis.c
	g++ -o $@ $^ $(CFLAGS) $(LDLIBS)

data_extractor: data_extractor.c
	g++ -o $@ $^ $(CFLAGS) $(LDLIBS)

clean:
	rm -f *.o energy_cut2 all_sipm_pair_time_differences twod_channel_coincidence branch_cutter voltages_and_time_differences peak_energy_and_over_voltage polygamma_fit

.PHONY: clean
