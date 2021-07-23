all: energy_cut2 all_sipm_pair_time_differences twod_channel_coincidence branch_cutter voltages_and_time_differences

energy_cut2: energy_cut2.c
	g++ -o energy_cut2 energy_cut2.c `root-config --cflags --libs`

all_sipm_pair_time_differences: all_sipm_pair_time_differences.c
	g++ -o all_sipm_pair_time_differences all_sipm_pair_time_differences.c `root-config --cflags --libs`

twod_channel_coincidence: twod_channel_coincidence.c
	g++ -o twod_channel_coincidence twod_channel_coincidence.c `root-config --cflags --libs`

branch_cutter: branch_cutter.c
	g++ -o branch_cutter branch_cutter.c `root-config --cflags --libs`

voltages_and_time_differences: voltages_and_time_differences.c
	g++ -o voltages_and_time_differences voltages_and_time_differences.c `root-config --cflags --libs`
