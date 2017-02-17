monoHZZ4L study for RGS paper based on Les Houches work 10-19 June 2015 HBP & Nicola
	pp -> 4 mu, or 4e, or 2e2mu

(1) SETUP
    
Edit setup.(c)sh as needed then

    source setup.(c)sh
    
Then
    make

(2) CREATE REDUCED NTUPLES FROM DELPHES FILES

To run the analysis do
    
    cd work
	source build.sh
	./mkNtuple -x cross-section -l integrated-luminosity Delphes-file
	
    
	    
    
