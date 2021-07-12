all: energy_cut2

energy_cut2: energy_cut2.c
	g++ -o energy_cut2 energy_cut2.c `root-config --cflags --libs`
