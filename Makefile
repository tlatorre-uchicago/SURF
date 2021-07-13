all: energy_cut2

energy_cut2: energy_cut2.c
	g++ -o energy_cut2 energy_cut2.c `root-config --cflags --libs`

all_sipm_pair_time_differences: all_sipm_pair_time_differences.c
	g++ -o all_sipm_pair_time_differences all_sipm_pair_time_differences.c `root-config --cflags --libs`
