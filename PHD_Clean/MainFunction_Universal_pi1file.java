/**
 * 
 */
package PHD_Clean;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.data.H3F;
//import org.jlab.io.base.DataEvent;
import org.jlab.jnp.hipo4.io.HipoReader;
import org.jlab.jnp.hipo4.io.HipoWriter;
import org.jlab.jnp.physics.EventFilter;
import org.jlab.jnp.physics.LorentzVector;
import org.jlab.jnp.physics.Particle;
import org.jlab.jnp.physics.PhysicsEvent;
import org.jlab.jnp.reader.DataManager;


import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo.data.HipoEvent;
import org.jlab.jnp.hipo.data.HipoGroup;
import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.SchemaFactory;
import org.jlab.jnp.hipo4.data.Schema.SchemaBuilder;
import org.jlab.jnp.hipo4.data.Schema;
/**
 * @author gangelini
 * Function to compute Pion BSA and Multiplicity raw counts (Counts in Phi, counts integrated in phi)
 */





public class MainFunction_Universal_pi1file {
	static LorentzVector smear(LorentzVector Vector){ 
		 Random r = new Random();

	    double inM = Vector.mass();
	    double sP  = Vector.p();
	    double sTh = Vector.theta();
	    double sPh = Vector.phi();

	
	    double sThD = Math.toDegrees(sTh);
	    double momS1 = 0.0184291 -0.0110083*sThD + 0.00227667*sThD*sThD -0.000140152*sThD*sThD*sThD + 3.07424e-06*sThD*sThD*sThD*sThD;
	    double momS2 = 0.02*sThD;
	    double momR  = 0.01 * Math.sqrt(  Math.pow(momS1*sP, 2) + Math.pow(momS2, 2));
	    momR *= 2.0; // <- only for data

	    double theS1 = 0.004*sThD + 0.1;
	    double theS2 = 0;
	    double theR  = Math.sqrt(Math.pow(theS1*Math.sqrt(sP*sP+0.13957*0.13957)/(sP*sP),2) + Math.pow(theS2,2));
	    theR *= 2.5; // <- only for data

	    double phiS1 = 0.85-0.015*sThD;
	    double phiS2 = 0.17-0.003*sThD;
	    double phiR  = Math.sqrt(Math.pow(phiS1*Math.sqrt(sP*sP+0.13957*0.13957)/(sP*sP),2) + Math.pow(phiS2,2) );
	    phiR *= 3.5; // <- only for data

	    sPh +=Math.toRadians(phiR)  * r.nextGaussian();
	    sTh += Math.toRadians(theR)  * r.nextGaussian();
	    sP  += momR  * r.nextGaussian() * Vector.p() ; 
	    LorentzVector FinalVector = new LorentzVector();
	    double px = sP* Math.sin(sTh)* Math.cos(sPh);
		double py = sP* Math.sin(sTh)* Math.sin(sPh);
		double pz = sP*Math.cos(sTh);
		FinalVector.setPxPyPzE(px, py, pz, Math.sqrt( sP*sP + inM*inM ) );
		return FinalVector;
	}
	
	// MultiBins is an object that allow to define custom binning up to 6Dimensions
	private static MultiBins InvariantMassBins,PionCounts,Counts_Phi,Helicity0,Helicity1,MissingM,MCParticles,Electron_counts ;

	// Analysis
	//0/1 (MC no, MC yes),0/1 (debug no, debug yes), 0/1  (Multi dimensional no=0 , yes =1 ), reading folder , writing folder 
	public static void main(String[] args) throws FileNotFoundException {
		
//ReadParameters Parameters = new ReadParameters("/Users/gangelini/testParameters.txt");	

	ReadParameters Parameters = new ReadParameters("Parameters.txt");
//ReadParameters Parameters = new ReadParameters("/lustre19/expphy/volatile/clas12/gangel/PHD_Analysis/Parameters.txt");	

int ContaStep0=0;
int ContaStep1=0;
int ContaStep2=0;
int ContaStep3=0;
int ContaStep4=0;
int ContaStep5=0;
int ContaStep6=0;
int ContaStep7=0;


		// Select the particle you want to analyze
		int Contaunelettrone =0 ;
		int Contatuttielettroni=0;
		//	int ContaTuttiIpioni=0;
		int ContaiPioniBuoni=0;
		//	int contapionebanca=0;

		int PionID=Parameters.GetPID();
		boolean MC=Parameters.GetMCBoolean();
		boolean debug=Parameters.Get_DebugBoolean();
		int goodeventsStop=10000000;
	
        String workdir =Parameters.Get_InputDir(); 
		String workdirout =Parameters.Get_OutputDir(); 
        int valueBinStatus=Parameters.GetIndexAnalysis();


    	 String ListRun45[]= {"skim4_005032.hipo","skim4_005036.hipo","skim4_005038.hipo","skim4_005039.hipo","skim4_005040.hipo","skim4_005041.hipo","skim4_005043.hipo","skim4_005045.hipo","skim4_005046.hipo","skim4_005047.hipo","skim4_005051.hipo","skim4_005052.hipo","skim4_005053.hipo","skim4_005116.hipo","skim4_005117.hipo","skim4_005119.hipo","skim4_005120.hipo","skim4_005124.hipo","skim4_005125.hipo","skim4_005126.hipo","skim4_005127.hipo","skim4_005128.hipo","skim4_005129.hipo","skim4_005130.hipo",
    	 "skim4_005137.hipo","skim4_005138.hipo","skim4_005139.hipo","skim4_005153.hipo","skim4_005158.hipo","skim4_005159.hipo","skim4_005160.hipo","skim4_005162.hipo","skim4_005163.hipo","skim4_005164.hipo","skim4_005165.hipo","skim4_005166.hipo","skim4_005167.hipo","skim4_005168.hipo","skim4_005169.hipo","skim4_005180.hipo","skim4_005181.hipo","skim4_005182.hipo","skim4_005183.hipo","skim4_005189.hipo","skim4_005190.hipo","skim4_005191.hipo","skim4_005193.hipo","skim4_005194.hipo",
    	 "skim4_005195.hipo","skim4_005196.hipo","skim4_005197.hipo","skim4_005198.hipo","skim4_005199.hipo","skim4_005200.hipo","skim4_005201.hipo","skim4_005202.hipo","skim4_005203.hipo","skim4_005204.hipo","skim4_005205.hipo","skim4_005206.hipo","skim4_005208.hipo","skim4_005202.hipo","skim4_005203.hipo","skim4_005204.hipo","skim4_005205.hipo","skim4_005206.hipo","skim4_005208.hipo","skim4_005211.hipo","skim4_005212.hipo","skim4_005215.hipo","skim4_005216.hipo","skim4_005219.hipo",
    	 "skim4_005220.hipo","skim4_005221.hipo","skim4_005222.hipo","skim4_005223.hipo","skim4_005225.hipo","skim4_005229.hipo","skim4_005230.hipo","skim4_005231.hipo","skim4_005232.hipo","skim4_005233.hipo","skim4_005234.hipo","skim4_005235.hipo","skim4_005237.hipo","skim4_005238.hipo","skim4_005239.hipo","skim4_005247.hipo","skim4_005248.hipo","skim4_005249.hipo","skim4_005250.hipo","skim4_005252.hipo","skim4_005253.hipo","skim4_005257.hipo","skim4_005258.hipo","skim4_005259.hipo",
    	 "skim4_005261.hipo","skim4_005262.hipo","skim4_005302.hipo","skim4_005303.hipo","skim4_005304.hipo","skim4_005305.hipo","skim4_005306.hipo","skim4_005307.hipo","skim4_005310.hipo","skim4_005311.hipo","skim4_005315.hipo","skim4_005316.hipo","skim4_005317.hipo","skim4_005318.hipo","skim4_005319.hipo","skim4_005320.hipo","skim4_005323.hipo","skim4_005324.hipo","skim4_005325.hipo","skim4_005333.hipo","skim4_005334.hipo","skim4_005335.hipo","skim4_005336.hipo","skim4_005339.hipo",
    	 "skim4_005340.hipo","skim4_005341.hipo","skim4_005346.hipo","skim4_005347.hipo","skim4_005349.hipo","skim4_005351.hipo","skim4_005354.hipo","skim4_005355.hipo","skim4_005367.hipo","skim4_005414.hipo","skim4_005415.hipo","skim4_005416.hipo"};
 

    	 String ListRun50[] = {"skim4_005341.hipo","skim4_005342.hipo","skim4_005343.hipo","skim4_005344.hipo","skim4_005345.hipo","skim4_005356.hipo","skim4_005357.hipo","skim4_005358.hipo","skim4_005359.hipo","skim4_005360.hipo","skim4_005361.hipo","skim4_005362.hipo","skim4_005366.hipo"};
    
    	 String  ListRun55[]= { "skim4_005368.hipo","skim4_005369.hipo","skim4_005370.hipo","skim4_005371.hipo","skim4_005372.hipo","skim4_005373.hipo","skim4_005374.hipo","skim4_005375.hipo","skim4_005376.hipo","skim4_005377.hipo","skim4_005378.hipo","skim4_005379.hipo","skim4_005380.hipo","skim4_005381.hipo","skim4_005382.hipo","skim4_005383.hipo","skim4_005386.hipo","skim4_005390.hipo","skim4_005391.hipo","skim4_005392.hipo","skim4_005393.hipo","skim4_005394.hipo","skim4_005398.hipo",
    	 "skim4_005399.hipo","skim4_005400.hipo","skim4_005401.hipo","skim4_005402.hipo","skim4_005403.hipo","skim4_005404.hipo","skim4_005406.hipo","skim4_005407.hipo"};
   	
 

   	       SchemaBuilder schemaBuilder = new SchemaBuilder("SIDIS::ePiX",120,1);    
   		   schemaBuilder.addEntry("runnumber", "I", "run number  id");
   		   schemaBuilder.addEntry("event", "L", "event  id");
   		   schemaBuilder.addEntry("helicity", "B", "event helicity");
           schemaBuilder.addEntry("pid", "I", "Pion id");
           schemaBuilder.addEntry("e_p", "D", "Trigger electron momentum");
           schemaBuilder.addEntry("e_theta", "D", "Trigger electron theta");
           schemaBuilder.addEntry("e_phi", "D", "Trigger electron phi ");
           schemaBuilder.addEntry("Q2", "D", "Q2 of the event");
           schemaBuilder.addEntry("W", "D", "W of the event");
           schemaBuilder.addEntry("x", "D", "x of the event");
           schemaBuilder.addEntry("y", "D", "y of the event");
           schemaBuilder.addEntry("epsilon", "D", "epsilon of the event");
           schemaBuilder.addEntry("pi_p", "D", "Pion momentum");
           schemaBuilder.addEntry("pi_theta", "D", "Pion theta");
           schemaBuilder.addEntry("pi_phi", "D", "Pion phi");
           schemaBuilder.addEntry("z", "D", "x z of the pion");
           schemaBuilder.addEntry("pt", "D", " PT  of the pion (with respect virtual photon)");
           schemaBuilder.addEntry("phi", "D", "Phi (trento) of the pion");
           schemaBuilder.addEntry("xf", "D", "x Feynman of the pion");
           schemaBuilder.addEntry("mX", "D", " Missing mass ehX");
           schemaBuilder.addEntry("eta", "D", "Rapidity in Breit frame");
           
           schemaBuilder.addEntry("pi_chi2pid", "D", "Pion_Chi2PID");
           schemaBuilder.addEntry("PCAL_lv", "D", "PCAL lv");
           schemaBuilder.addEntry("PCAL_lw", "D", "PCAL lw");
           schemaBuilder.addEntry("DC1_x", "D", "DC Layer 1 x hit");
           schemaBuilder.addEntry("DC1_y", "D", "DC Layer 1 x hit");
           schemaBuilder.addEntry("DC2_x", "D", "DC Layer 2 x hit");
           schemaBuilder.addEntry("DC2_y", "D", "DC Layer 2 y hit");
           schemaBuilder.addEntry("DC3_x", "D", "DC Layer 3 x hit");
           schemaBuilder.addEntry("DC3_y", "D", "DC Layer 3 y hit");
           schemaBuilder.addEntry("DC2_sector", "I", "Sector in DC layer2");
           
           schemaBuilder.addEntry("generated_id", "I", "Generated PID");
           schemaBuilder.addEntry("generated_pi_p", "D", "Generated momentum");
           schemaBuilder.addEntry("generated_pi_theta", "D", "Generated theta");
           schemaBuilder.addEntry("getenerated_pi_phi", "D", "Generated phi");
           schemaBuilder.addEntry("generated_z", "D", "Generated Z");
           schemaBuilder.addEntry("generated_pt", "D", "Generated PT");
           schemaBuilder.addEntry("generated_phi", "D", "Generated Phi Trento");
           schemaBuilder.addEntry("generated_xf", "D", "Generated XF");
           schemaBuilder.addEntry("generated_mX", "D", "Generated Missing Mass ehX");
           schemaBuilder.addEntry("generated_eta", "D", "Generated Eta Bright Frame ");
           schemaBuilder.addEntry("status", "I", "Status of the cuts");
           
           
       
           
           
   		  Schema schemaW  = schemaBuilder.build();
   		  HipoWriter  writeroutput = new HipoWriter();
          writeroutput.getSchemaFactory().addSchema(schemaW); // Should I use schema or factory? 
          writeroutput.open("Skimmed_Hipo.hipo");
         // writeroutput.open(workdirout+"/Skimmed_Hipo.hipo");

           boolean pi0 =false;
           if(Parameters.GetPID()==111)pi0=true;
		
		   boolean PolygonalBinClas=true; // Define if for x and Q2 we arew using the polygon class 
		   boolean SingleHadron= false;

		   boolean Q21D=false; boolean xB1D=false; boolean z1D=false; boolean PT1D=false;
		   boolean Qx2D=false; boolean QxzPt=false;


		if (valueBinStatus==0) {Q21D=true;PolygonalBinClas=false;}
		if (valueBinStatus==1) {xB1D=true;PolygonalBinClas=false;}
		if (valueBinStatus==2) {z1D=true;PolygonalBinClas=false;}
		if (valueBinStatus==3) {PT1D=true;PolygonalBinClas=false;}
		if (valueBinStatus==4) {Qx2D=true;PolygonalBinClas=true;}
		if (valueBinStatus==5) {QxzPt=true;PolygonalBinClas=true;}

// ATTENTION Z MIN HERE
		int Q2_Bins=1; double Q_min= Parameters.GetQ2LowLimit(); double Q_max= Parameters.GetQ2HighLimit(); // 1 - 12 		  
		int Xb_Bins=1; double Xb_min= Parameters.GetxBLowLimit() ; double Xb_max= Parameters.GetxBHighLimit();  // 0 -0.9
		int Z_Bins= 7;  double Z_min= Parameters.GetZLowLimit();   double Z_max= Parameters.GetZHighLimit();  // 0.2 - 0.9
		int PT_Bins=1; double PT_min= 0;  double PT_max= 2;// 0.1 - 4
		int Mass_Bins= 100; double min_pion= 0.05;  double max_pion=0.25;	 
		int Phi_Bins= 12; double Phi_min=0.0; double Phi_max=360.0;

		//Common Bins Phi and Missing mass
		int bincount_phi = 12;
		double bin_phi[] = {0.0, 30.0, 60.0, 90.0, 120.0, 150.0, 180.0, 210.0, 240.0, 270.0, 300.0, 330.0, 360.0};
		// To be incremented to 100 by using Axis class and self division:
		int binmass=300;
		double []bin_mass= new double[300];
		for(int i=0; i<300; i++) {
			double binsize =0.25/299 ;
			double minvalue=0.05;
			double conto=minvalue+binsize*i;
			bin_mass[i]=conto;
		}

		if (valueBinStatus==0) {
			int bincount_q2 = 1;
			double bin_q2[] = {1.0, 11.0};
			//int bincount_q2 = 13;
			//double bin_q2[] = {1.3, 1.65, 1.85, 2.05, 2.3, 2.6, 3.0, 3.5, 4.1, 4.9, 5.8, 6.8, 8.0, 11.0};
			int bincount_x = 1;
			double bin_x[] = {Xb_min,Xb_max};
			int bincount_z = 1;
			double bin_z[] = {Z_min,Z_max};
			int bincount_pt = 1;
			double bin_pt[] = {PT_min,PT_max};

			InitializeMultiBins( PolygonalBinClas,  bincount_q2, bincount_x,  bincount_z,  bincount_pt,  binmass,  bincount_phi,
					bin_q2,bin_x,  bin_z,  bin_pt, bin_mass,  bin_phi);
		}

		if (valueBinStatus==1) {
			int bincount_q2 = 1;
			double bin_q2[] = {Q_min,Q_max};
			int bincount_x = 11;
			double bin_x[] = {0.09, 0.145, 0.18, 0.215, 0.25, 0.29, 0.33, 0.37, 0.42, 0.48, 0.55, 0.70};
			int bincount_z = 1;
			double bin_z[] = {Z_min,Z_max};
			int bincount_pt = 1;
			double bin_pt[] = {PT_min,PT_max};

			InitializeMultiBins( PolygonalBinClas,  bincount_q2, bincount_x,  bincount_z,  bincount_pt,  binmass,  bincount_phi,
					bin_q2,bin_x,  bin_z,  bin_pt, bin_mass,  bin_phi);
		}

		if (valueBinStatus==2) {
			int bincount_q2 = 1;
			double bin_q2[] = {Q_min,Q_max};
			int bincount_x = 1;
			double bin_x[] = {Xb_min,Xb_max};
			int bincount_z = 13;
			double bin_z[] = {0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75,0.80,0.85};
			int bincount_pt = 1;
			double bin_pt[] = {PT_min,PT_max};
			InitializeMultiBins( PolygonalBinClas,  bincount_q2, bincount_x,  bincount_z,  bincount_pt,  binmass,  bincount_phi,
					bin_q2,bin_x,  bin_z,  bin_pt, bin_mass,  bin_phi);
		}
		if (valueBinStatus==3) {
			int bincount_q2 = 1;
			double bin_q2[] = {Q_min,Q_max};
			int bincount_x = 1;
			double bin_x[] = {Xb_min,Xb_max};
			int bincount_z = 1;
			double bin_z[] = {Z_min,Z_max};
			int bincount_pt = 14;
			double bin_pt[] = {0.0, 0.12, 0.22, 0.295, 0.36, 0.425, 0.5, 0.595, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.7};

			InitializeMultiBins( PolygonalBinClas,  bincount_q2, bincount_x,  bincount_z,  bincount_pt,  binmass,  bincount_phi,
					bin_q2,bin_x,  bin_z,  bin_pt, bin_mass,  bin_phi);
		}			


		if(valueBinStatus==4) {
			int bincount_q2 = 6;
			double bin_q2[] = {Q_min,Q_max};
			int bincount_x = 7;
			double bin_x[] = {Xb_min,Xb_max};
			int bincount_z = 20; 
			double bin_z[] = {0.20,0.23,0.26,0.29,0.32,0.35,0.38,0.41,0.44,0.47,0.50,0.53,0.56,0.59,0.62,0.65,0.68,0.71,0.74,0.77,0.80};
			int bincount_pt = 1; 
			double bin_pt[] = {PT_min,PT_max};

			InitializeMultiBins( PolygonalBinClas,  bincount_q2, bincount_x,  bincount_z,  bincount_pt,  binmass,  bincount_phi,
					bin_q2,bin_x,  bin_z,  bin_pt, bin_mass,  bin_phi);

		}

		
		if(valueBinStatus==5) {
			System.out.println(" I am here ");
			int bincount_q2 = 6;
			double bin_q2[] = {Q_min,Q_max};
			int bincount_x = 7;
			double bin_x[] = {Xb_min,Xb_max};
			//int bincount_z = 7; 
			//double bin_z[] = {0.15, 0.2, 0.24, 0.29, 0.36, 0.445, 0.55, 0.7};
		//	int bincount_pt =7; 
			//double bin_pt[] = {0.05, 0.2, 0.3, 0.4, 0.5, 0.6, 0.75, 1.0};
			int bincount_z = 12;
			double bin_z[] = {0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.75,0.80};
		//	int bincount_z = 5;
		//	double bin_z[] = {0.20,0.30, 0.40, 0.50, 0.60, 0.70};
			int bincount_pt =28; 
			double bin_pt[] = {0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95, 1.0,1.05,1.1,1.15,1.2,1.25,1.3,1.35,1.4};
		//	int bincount_pt =11; 
		//	double bin_pt[] = {0.0,0.12,0.24,0.36,0.48,0.60,0.72,0.84,0.96,1.08, 1.20,1.32};
		
			InitializeMultiBins( PolygonalBinClas,  bincount_q2, bincount_x,  bincount_z,  bincount_pt,  binmass,  bincount_phi,
					bin_q2,bin_x,  bin_z,  bin_pt, bin_mass,  bin_phi);

		}

		//This is PT2 
		if(PolygonalBinClas==true) {
			// We will create extra bins , THe overflow and udnerfloow will be in Q2 =4 and xb 7.
			//The only utalize are 2 bins of Q2 and 4 of xB 
			Q2_Bins=4; Q_min= 1;     Q_max= 12; // 1 - 12 		  
			Xb_Bins=7; Xb_min= 0.01; Xb_max= 0.9 ;  // 0 -0.9
			Z_Bins= 1; Z_min= 0.2;   Z_max= 0.8;  // 0.2 - 0.9
			PT_Bins=1;  PT_min= 0;   PT_max= 2;// 0.1 - 4
		}



		// Used to look into invariant mass distribution
		H1F XMass = new H1F("XMass",200,0,2);

 
	//	Histos Histo = new Histos(workdirout); // To save on file electron and photons distributions
		long count =0 ; // Total Event counters 
		long good_e_count=0; // Selected electrons
		long all_electrons=0; // Total Electrons reconstructed in events 
		long twophotons=0; //Events with 2 photons candidates to pi0 
		long ContaElDIS=0;
		long GoodRangeEl=0;

		FiducialCuts SIDIS_Cuts = new FiducialCuts(Parameters);  // This creates a fiducial cut with standardized cuts for my analysis
		SIDIS_Cuts.GetMap().show(); // print the map of cuts ssociate with each PID. It removes Data in region of detector that cannot be trusted. 

		int helicity=-9; //  Helicity inizialization  
		double beamEnergy = 10.603; // Set the energy of the beam
		LorentzVector beam   = new LorentzVector(0.0,0.0,beamEnergy,beamEnergy); // Getting the Lorentz Vector associated to the beam
		LorentzVector target = new LorentzVector(0.0,0.0,0.0,0.93827); // Defining Lorentz Vector for the proton target
		
		//Filtering the events properly
		EventFilter filter   = new EventFilter();
		if (Parameters.GetPID()==111 ) filter.setFilter("11:22:22:X+:X-:Xn"); // Search in the events for electron , 2 photons , and any other particle 
		else if (Parameters.GetPID()==-211 ) filter.setFilter("11:-211:X+:X-:Xn");
		else if (Parameters.GetPID()==211) filter.setFilter("11:211:X+:X-:Xn");
		System.out.println(" I am filtering the events using the scheme: ");System.out.println(filter);

		// Inizializing the Particles and Analysis tools 
		Pi0_Particle New_Pi0s = new Pi0_Particle(); // Pi0Particle 
	
		// Routine for the SIDIS analysis:
		AnalysisSIDIS Analysis = new AnalysisSIDIS(workdirout,PolygonalBinClas,Z_Bins, Z_min, Z_max, PT_Bins, PT_min, PT_max, SingleHadron); // Create an Analysis of Pi0
		
		// Analysis on the MC generated particles: 
		MCAnalysis LUND_pi = new MCAnalysis(PolygonalBinClas,InvariantMassBins, Z_Bins, Z_min, Z_max, PT_Bins, PT_min, PT_max, Phi_Bins); // Create an Analysis of MonteCarlo files
		LUND_pi.setBeam(beam);LUND_pi.setTarget(target); //Used to compute the PT, Phi, etc. from lund momenta

		long evento=0;
		int rnmb=0;
		long eventid = 0;
		
		File file = new File("input.hipo");
	//	File file = new File("/Volumes/My Passport for Mac/run5045/skim4_005045.hipo");
				if (file.isFile()) { 
					
					System.out.println(" Reading  " + file.getName());
					HipoReader reader = new HipoReader(); reader.open(file.getName());
					Event event=new Event(); 
					SchemaFactory schema = reader.getSchemaFactory();
					Bank RRun = new Bank(schema.getSchema("RUN::config"));
					Bank REvent = new Bank(schema.getSchema("REC::Event"));
					Bank RCal = new Bank(schema.getSchema("REC::Calorimeter"));
					Bank RPart = new Bank(schema.getSchema("REC::Particle"));
					Bank RTraj = new Bank(schema.getSchema("REC::Traj"));
					Bank RTrack = new Bank(schema.getSchema("REC::Track"));
					Bank RCherenkov = new Bank(schema.getSchema("REC::Cherenkov"));	
					
					

					
					double LorentzProduct=0, xB=0, Pq=0, Pl=0, y_var=0,epsilon=0,epsilon2= 0,gamfactor=0,epsilonnum=0,epsilonden=0;
					double e_mom=0, e_Theta =0,startime=0; 
					while(reader.hasNext()==true){
						eventid++;
						//System.out.println("---------------------");
						count++;
						boolean checkEle = false;
						reader.nextEvent(event);
						event.read(REvent); event.read(RCal); event.read(RPart);event.read(RTraj); event.read(RRun);event.read(RCherenkov);
						event.read(RTrack);
						if(Parameters.GetMCBoolean()) { evento = eventid;}
						else {evento =(long) RRun.getInt("event", 0);}
					
						rnmb=RRun.getInt("run", 0);
						if(count%100000 == 0) System.out.println(count + " events");
						if(good_e_count%50000 == 0) System.out.println(good_e_count + " good electrons");	
						helicity= REvent.getByte("helicity",0);
						startime=(double)REvent.getFloat("startTime", 0);
						double p_pi =0;
						int eventlist=0;
						if(RPart.getRows()>0 && RCal.getRows()>0 && RTraj.getRows()>0) // Be sure I have all the banks I need
						{
							PhysicsEvent phsEvent = DataManager.getPhysicsEvent(beamEnergy,RPart) ;
							int TriggerParticleID = RPart.getInt("pid",0); //Get the ID of the trigger
							//Check on the HTTC signal: 
							boolean HTTCSignal=false;
							double nphe =0;
							for(int i=0; i<RCherenkov.getRows(); i++) {
								if( RCherenkov.getInt("pindex", i)==0 && RCherenkov.getInt("detector", i)==15 && RCherenkov.getFloat("nphe", i) >2 ) {
									nphe=RCherenkov.getFloat("nphe", i);
									HTTCSignal=true;
								}
							}
							//RPart.show();
							//RTrack.show();
							//RTraj.show();
							//System.out.println("--------");
							double chi2pidelectron=(double)RPart.getFloat("chi2pid", 0);
							int statuselectron=(int)RPart.getShort("status", 0);
							double min_mom=0, Q2=0, multifact=0;
							int indicepione=-9;
							
							
							//==> First selection : Filter, 1 electron only; Trigger is electron , HTTC NPE>2 
						//	if( chi2pidelectron<3 && statuselectron<0 && filter.isValid(phsEvent)==true && TriggerParticleID==11 && phsEvent.countByPid(11)==1 && HTTCSignal==true) {
								//if( chi2pidelectron<3 && statuselectron<0 && filter.isValid(phsEvent)==true && TriggerParticleID==11  && HTTCSignal==true) {
							
							//if( TriggerParticleID==11 && (RPart.getShort("status", 0))>-4000 && (RPart.getShort("status", 0))<-2000 ) {
								//ContaStep1++;
								//Read the electron:
								Particle electron = phsEvent.getParticleByPid(11, 0);
								// Extend the electron by defining a Reconstructed particle, menaing attaching to it the Part,Trajectory and Calorimeter hits.
								ParticleREC electronRec = new ParticleREC(electron,RPart,RTraj,RCal);
								all_electrons++;// I found an electron 		
								LorentzVector vecW2 = new LorentzVector(0,0,0,0);
								LorentzVector vecQ2 = new LorentzVector(0,0,0,0);
								LorentzVector vettoreE = LorentzVector.from(electron.vector()); // electron.vector() returns the LorentzVector of the particle electron
								LorentzVector vecE= new LorentzVector(0,0,0,0);
									//Without smearing
								vecE=vettoreE;
								
								//if(MC=true) { vecE = smear(vettoreE); }
								//else { vecE = vettoreE;} // electron.vector() returns the LorentzVector of the particle electron
								vecW2.add(target);vecW2.add(beam);vecW2.sub(vecE);vecQ2.add(beam); vecQ2.sub(vecE);
								Pq = -target.px()*vecQ2.px()-target.py()*vecQ2.py()-target.pz()*vecQ2.pz()+target.e()*vecQ2.e();
								Pl = -target.px()*beam.px()-target.py()*beam.py()-target.pz()*beam.pz()+target.e()*beam.e();
								xB = -vecQ2.mass2()/(2*Pq); 
								y_var = Pq/Pl;
								Q2= -vecQ2.mass2();
								epsilon2 = ( 1 - y_var - Q2 / (4*10.603*10.603) ) / ( 1 - y_var + 0.5*y_var*y_var + Q2 / (4*10.603*10.603) );
								gamfactor=(2*0.93827*xB)/Math.sqrt(-vecQ2.mass2());
								multifact= (0.25)*Math.pow(y_var, 2)* Math.pow(gamfactor, 2);
								epsilonnum=1-y_var- multifact ;
								epsilonden=1-y_var+(0.5)*Math.pow(y_var, 2)+multifact;
								epsilon=epsilonnum/epsilonden;
								
								e_mom=electron.p();
								e_Theta = Math.toDegrees(electron.theta());
								double e_TrackChi2=RTrack.getFloat("chi2", 0);
								double e_TackNDF=(double)RTrack.getShort("NDF", 0);
								double e_TrackRedChi2 = e_TrackChi2/e_TackNDF;
								
								//Histo.SetUnSkim_el(electronRec,-vecQ2.mass2(),xB,vecW2.mass(),y_var);

								if (vecW2.mass()>= Parameters.GetW2Cut() && -vecQ2.mass2()>= Parameters.GetQ2Cut())ContaElDIS++;

								//===> 2 Cut: Kinemaitcs and electron angle and momentum
								if( TriggerParticleID==11) {
								int sector_e=0;
								if (electronRec.hasDetector(71))sector_e= electronRec.getSector(71);
								else if (electronRec.hasDetector(71)==false && electronRec.hasDetector(72)) sector_e= electronRec.getSector(72);
								boolean GoodVertex=false;
								double Minvertex=0,Maxvertex=0, Sigma=Parameters.Get_el_VertexSigma();
								if(Parameters.GetMCBoolean()==true) {
								Minvertex =(-3.2)-5.0+Sigma*(0.48+(0.008-0.005*electron.p())*Math.toDegrees(electron.theta())) ;
								Maxvertex= (-3.2)+5.0+Sigma*(0.48+(0.008-0.005*electron.p())*Math.toDegrees(electron.theta())) ;
								}
								else {
									if (sector_e==1 )
									{
										Minvertex = (-3.13)-5.2-(0.95+(-0.009-0.003*electron.p())*Math.toDegrees(electron.theta()));
										Maxvertex = (-3.13)+5.2+(0.95+(-0.009-0.003*electron.p())*Math.toDegrees(electron.theta()));
									}
									else if (sector_e==2)
									{
										Minvertex = (-2.9)-5.2-(0.95+(0.019-0.009*electron.p())*Math.toDegrees(electron.theta()));
										Maxvertex = (-2.9)+5.2+(0.95+(0.019-0.009*electron.p())*Math.toDegrees(electron.theta()));
									}
									else if (sector_e==3)
									{
										Minvertex = (-3.3)-5.2-(0.95+(0.016-0.009*electron.p())*Math.toDegrees(electron.theta()));
										Maxvertex = (-3.3)+5.2+(0.95+(0.016-0.009*electron.p())*Math.toDegrees(electron.theta()));
									}
									else if (sector_e==4)
									{
										Minvertex = (-2.91)-5.2-(0.95+(0.027-0.011*electron.p())*Math.toDegrees(electron.theta()));
										Maxvertex = (-2.91)+5.2+(0.95+(0.027-0.011*electron.p())*Math.toDegrees(electron.theta()));
									}
									else if (sector_e==5)
									{
										Minvertex = (-3.15)-5.2-(0.95+(0.006-0.006*electron.p())*Math.toDegrees(electron.theta()));
										Maxvertex = (-3.15)+5.2+(0.95+(0.006-0.006*electron.p())*Math.toDegrees(electron.theta()));
									}
									else if (sector_e==6)
									{
										Minvertex = (-3.18)-5.2-(0.95+(-0.007-0.004*electron.p())*Math.toDegrees(electron.theta()));
										Maxvertex = (-3.18)+5.2+(0.95+(-0.007-0.004*electron.p())*Math.toDegrees(electron.theta()));
									}
								}
								if(electron.vz()>Minvertex && electron.vz()<Maxvertex) GoodVertex=true;
								
								if(vecW2.mass()>= Parameters.GetW2Cut() && -vecQ2.mass2()>= Parameters.GetQ2Cut() && y_var<= Parameters.GetYCut() && GoodVertex && Math.toDegrees(electron.theta())>Parameters.Get_el_minTheta()&& Math.toDegrees(electron.theta())<Parameters.Get_el_maxTheta())
								{    	
								   //Wrong Phi to be fixed 
									boolean Phicut=false;
									if(Parameters.GetPhiCutBoolean()) {
									double ephi=(electron.phi()*57.3+30+180)%60.0;
									if(ephi>(0+Parameters.GetPhiCutRange()) && ephi<(60-Parameters.GetPhiCutRange())) Phicut=true;
									}
									else Phicut=true;
								//	if(Phicut)System.out.println(" Phi cut is true");
								
								//	if(Phicut==false) System.out.println(" Phi cut "+ Phicut);
									Contatuttielettroni++;
									if(phsEvent.countByPid(11)>1) Contaunelettrone++;	
									
									//Checking that the electron passes the fiducial cuts:
								
									boolean GoodElectron=false;
									 if(Parameters.Get_el_FC()==false ) GoodElectron=true;
									 else if(Parameters.Get_el_FC()==true && SIDIS_Cuts.Status(electronRec)==true )  GoodElectron=true;
									
									
									// FC and Phi Shaving
									if(GoodElectron && Phicut) {
									//To reactivate 
										//if(electronRec.p()>3 && electronRec.p()<6 && Math.toDegrees(electronRec.theta())>15 && Math.toDegrees(electronRec.theta())<25)	GoodRangeEl++;
										
									 // Histo.SetSkim_el(electronRec,-vecQ2.mass2(),xB,vecW2.mass(),y_var, LUND_pi.getQ2(), LUND_pi.getXB(), LUND_pi.getW(), LUND_pi.getY(),e_TrackChi2,e_TackNDF,chi2pidelectron);

										if (MC==true){
											//to reactivate 
										Bank MCParticle = new Bank(schema.getSchema("MC::Particle"));
										Bank MCLund = new Bank(schema.getSchema("MC::Lund"));
										event.read(MCParticle); event.read(MCLund);
										if(pi0==true) LUND_pi.add(Parameters,MCParticle,MCLund, 111,evento); 
										else  LUND_pi.add(Parameters,MCParticle,MCLund, PionID,evento); 
										LUND_pi.setParent(0);
									}
										Electron_counts.FillH2(-vecQ2.mass2(),xB);// Filling the electron counts 
										good_e_count++; //Good electron SISDIS within all cuts	
										//Check the number of pions in the event
										
									if(pi0==true) {
										// Define Photons
										Photons photons = new Photons(electronRec,RPart,RCal);
										//Histo.SetUnSkim_ph(photons);	
										if (photons.getStatus(2)==true) {
											//System.out.println(" ho due fotoni");
									   // Histo.SetSkim_ph(photons);
										New_Pi0s.add(photons);
										// Method 1 consider the combination of all photons
										// Method 2 consider only the two most energetiec photons 
										New_Pi0s.setMethod(Parameters.Get_Gamma_Method());
										New_Pi0s.setEnergyCut(Parameters.Get_Gamma_EnergyCut()); //minimum energy for the photons to consider
										// 2nd argument of the next fuction : Set a minimum angle cut for defying pi0s (if > ok, if < not considered) 
										New_Pi0s.getPi0s(beam,target,vecE,vecQ2,vecW2,xB,-vecQ2.mass2(),3);
										ContaiPioniBuoni= ContaiPioniBuoni+ New_Pi0s.getNrPi0s();
										New_Pi0s.resetPhotons();
										if(New_Pi0s.getNrPi0s()>0) {
										Analysis.Analyze(PolygonalBinClas,beam,target,electronRec,helicity,XMass,InvariantMassBins,MissingM,3,xB,-vecQ2.mass2(), New_Pi0s, PionCounts,Counts_Phi, Helicity0, Helicity1,Electron_counts,vecW2.mass(),y_var,epsilon, evento);
										}	
										// leggere 
										
									}
									}	
									else if(pi0!=true) {
										
										ArrayList<Double> TuplHadron_momentum = new ArrayList<Double>();
										ArrayList<Double> TuplHadron_theta = new ArrayList<Double>();
										ArrayList<Double> TuplHadron_phiclas = new ArrayList<Double>();
										ArrayList<Double> TuplHadron_z = new ArrayList<Double>();
										ArrayList<Double> TuplHadron_pt = new ArrayList<Double>();
										ArrayList<Double> TuplHadron_phi = new ArrayList<Double>();
										ArrayList<Double> TuplHadron_xf = new ArrayList<Double>();
										ArrayList<Double> TuplHadron_MX = new ArrayList<Double>();
										ArrayList<Double> TuplHadron_eta = new ArrayList<Double>();
										ArrayList<Double> expectedZ= new ArrayList<Double>();
										ArrayList<Double> expectedPT= new ArrayList<Double>();
										ArrayList<Double> expectedTheta= new ArrayList<Double>();
										ArrayList<Double> expectedPhiClas12= new ArrayList<Double>();
										ArrayList<Double> expectedPhiTrento= new ArrayList<Double>();
										ArrayList<Double> expectedMomentum= new ArrayList<Double>();
										ArrayList<Double> expectedxF= new ArrayList<Double>();
										ArrayList<Double> expectedeta= new ArrayList<Double>();
										ArrayList<Double> expectedmx= new ArrayList<Double>();

										ArrayList<Integer> expectedPid= new ArrayList<Integer>();
										
										
										ArrayList<Double> TuplHadron_DC1x= new ArrayList<Double>();
										ArrayList<Double> TuplHadron_DC1y= new ArrayList<Double>();
										ArrayList<Double> TuplHadron_DC2x= new ArrayList<Double>();
										ArrayList<Double> TuplHadron_DC2y= new ArrayList<Double>();
										
										ArrayList<Double> TuplHadron_DC3x= new ArrayList<Double>();
										ArrayList<Double> TuplHadron_DC3y= new ArrayList<Double>();
										ArrayList<Integer> TuplHadron_DC2sector= new ArrayList<Integer>();
										
										
										ArrayList<Double> TuplHadron_chi2pid= new ArrayList<Double>();

								
										
										
										int npioni = phsEvent.countByPid(PionID);
										for (int kk=0; kk<npioni ; kk++) {
											Particle pion_positive2 = phsEvent.getParticleByPid(PionID, kk);
											//Defined the Particle REC with the info from the other banks attached
											LorentzVector Pione_b_smearing = new LorentzVector();// PHM is the lorentz vector of the pion
											Pione_b_smearing.setPxPyPzE(pion_positive2.px(), pion_positive2.py(), pion_positive2.pz(), pion_positive2.e());

											LorentzVector Pion4v= new LorentzVector(0,0,0,0);
											//Without smearing
											//Pion4v=Pione_b_smearing;
											
										if(MC=true) { Pion4v = smear(Pione_b_smearing); }
										else { Pion4v = Pione_b_smearing;} // electron.vector() returns the LorentzVector of the particle electron
											
										pion_positive2.setP(Pion4v.p());
										pion_positive2.setTheta(Pion4v.theta());
											
											
											ParticleREC Pione = new ParticleREC(pion_positive2,RPart,RTraj,RCal);
											Pione.changePid(PionID); //Be sure that the PID is set to 211
											// Looping over the REC::PARTICLE bank to identify the pion I am looking at:
											//
											// For the purpose of Chi2PID 
											boolean Chi2PIDBoolean=false;
											double H_TrackChi2=(double)RTrack.getFloat("chi2", Pione.index);
											double H_TackNDF=(double)RTrack.getShort("NDF", Pione.index);
											
											double H_TrackRedChi2 = H_TrackChi2/H_TackNDF;
											
											double Chi2Pid2= (double) RPart.getFloat("chi2pid",Pione.index);
										
											if(Parameters.Get_Hadron_Ch2PIDBoolean()) {
											double chi2pid_max=0;
											double chi2_pid_limit_lower= -3.0*0.88;
											if(Pione.p() < 2.44){ chi2pid_max = 3;}
											else{ chi2pid_max = (0.00869 + 14.98587 * Math.exp(-Pione.p() / 1.18236) + 1.81751 * Math.exp(-Pione.p() / 4.86394));}
											chi2pid_max=0.88*chi2pid_max;
											if(Chi2Pid2 > chi2_pid_limit_lower && Chi2Pid2 < chi2pid_max)Chi2PIDBoolean=true;
											else Chi2PIDBoolean=false;
											}
											else  { Chi2PIDBoolean=true;}
										//	if(Chi2PIDBoolean==false) System.out.println("Chi2PIDBoolean false");
											// Fiducial Cuts : 
											//for now only
											Chi2PIDBoolean=true;
											
											boolean GoodPi=false;
											
											if(Parameters.Get_Hadron_FC() )GoodPi=SIDIS_Cuts.Status(Pione); //Boolean for Fiducial cuts 
											else GoodPi=true;
								    		//if(GoodPi==false) System.out.println("Fiducial Cuts Pions are false");

											int Sigma_vertexH= Parameters.Get_Hadron_VertexSigma();
											double H_Minvertex=0,H_Maxvertex=0;
											if (Parameters.GetMCBoolean()==true) {
												H_Minvertex =-2.7-0.01*Math.toDegrees(Pione.theta())-2.5-Sigma_vertexH*(0.7-0.02*Math.toDegrees(Pione.theta()));
												H_Maxvertex= -2.7-0.01*Math.toDegrees(Pione.theta())+2.5+Sigma_vertexH*(0.7-0.02*Math.toDegrees(Pione.theta()));
											}
											else {
												H_Minvertex =-4.22+0.034*Math.toDegrees(Pione.theta())-2.5-Sigma_vertexH*(1.9-0.05*Math.toDegrees(Pione.theta()));
												H_Maxvertex= -4.22+0.034*Math.toDegrees(Pione.theta())+2.5+Sigma_vertexH*(1.9-0.05*Math.toDegrees(Pione.theta()));
											}
											
											//&&New_Pi.getMissMass()>Cut_MissingMass_H
												if(Pione.vz()<(H_Maxvertex) && Pione.vz()>(H_Minvertex) && Chi2PIDBoolean  && GoodPi && H_TackNDF>Parameters.Get_Hadron_TrackNDF_Min() && H_TackNDF<Parameters.Get_Hadron_TrackNDF_Max() ) {	
												
													//System.out.println("Good pion");
												//++;
												// For this good Pion generated trhe class ChargedParitcle, meaning get Z,PT,PHI, Pseudorapidity etc. 
												ChargedPi_Particle New_Pi= new ChargedPi_Particle(); // Charged Pi class
												New_Pi.getChargedPi(phsEvent,pion_positive2, beam,target,vecE,vecQ2,vecW2,xB,-vecQ2.mass2(),0.0);

												
												
												//Be sure the Z cut applies only at the integrated Z cases 
												boolean ZCUT=false;		    		 
												if(valueBinStatus==2 || valueBinStatus==5)ZCUT=true;
												else if(valueBinStatus!=2 && valueBinStatus !=5 && New_Pi.getZ()>Z_min)ZCUT=true;
												
												// Status Forward 
												boolean StatusForward=false;
												if(Parameters.Get_Hadron_Forward() && (RPart.getShort("status", Pione.index)%4000/2000)>0)StatusForward=true;
												else if (Parameters.Get_Hadron_Forward()==false) StatusForward=true;
												//&&New_Pi.getMissMass()>Cut_MissingMass_H

														if(StatusForward ) {		
	
														if( Parameters.GetMCBoolean()==true) {
													//	if(New_Pi.getZ()>0.6 && xB>0.4 && Q2>4 && New_Pi.getPt()>0.4)System.out.println(" ----- > "+ New_Pi.getMissMass()+ " theta "+New_Pi.getTheta()+" phi "+ New_Pi.getPhi()+" momentum" + New_Pi.getMomentum() );
													ContaiPioniBuoni++; // Count the good pions
												//	System.out.println(" ------ " );

													//System.out.println(New_Pi.getMomentum());
													TuplHadron_momentum.add(New_Pi.getMomentum());
													TuplHadron_theta.add(New_Pi.getTheta());
													TuplHadron_phiclas.add(New_Pi.getPhiClas12());
													TuplHadron_z.add(New_Pi.getZ());
													TuplHadron_pt.add(New_Pi.getPt());
													TuplHadron_phi.add(New_Pi.getPhi());
													TuplHadron_xf.add(New_Pi.getXF());
													TuplHadron_MX.add(New_Pi.getMissMass());
													TuplHadron_eta.add(New_Pi.getEtaBar());
															
													TuplHadron_DC1x.add(Pione.getDetectorHits(61).x());
													TuplHadron_DC1y.add(Pione.getDetectorHits(61).y());								
													TuplHadron_DC2x.add(Pione.getDetectorHits(62).x());
													TuplHadron_DC2y.add(Pione.getDetectorHits(62).y());								
													TuplHadron_DC3x.add(Pione.getDetectorHits(63).x());
													TuplHadron_DC3y.add(Pione.getDetectorHits(63).y());
													TuplHadron_DC2sector.add(Pione.getSector(62));
								
													TuplHadron_chi2pid.add(Chi2Pid2);
													
													// Looking at the reconstructed particles 
													double small=99; double smallEn=99; 
													int smallindex=-3; int enindex=-3; 
													for(int k =0 ; k<LUND_pi.Get_H_theta().size(); k++ ) {
													if(LUND_pi.Get_H_vz().get(k)>-6 && LUND_pi.Get_H_vz().get(k)>-6 && LUND_pi.Get_H_PID().get(k)!=-PionID ) {
															double normThetaDiff=Math.abs((LUND_pi.Get_H_theta().get(k) - New_Pi.getTheta())/LUND_pi.Get_H_theta().get(k));
															double normMomDiff=Math.abs((LUND_pi.Get_H_mom().get(k) - New_Pi.getMomentum())/LUND_pi.Get_H_mom().get(k));
															double normPhiDiff=Math.abs((LUND_pi.Get_H_phi().get(k) - New_Pi.getPhiClas12())/LUND_pi.Get_H_phi().get(k));
															double difference = normThetaDiff+normMomDiff+normPhiDiff;
															double endiff = Math.abs(LUND_pi.Get_H_En().get(k)-New_Pi.getEnergy() );
														//	System.out.println(" Difference in Momentum " + normMomDiff  + " PID " +LUND_pi.Get_H_PID().get(k) );
															if (endiff<smallEn) { smallEn=endiff;  enindex=k;}

															if (difference<small) {small = difference; smallindex=k;}	
													  }
													}
										
													if(smallindex==-3) {
														smallindex=enindex;
													}
														
														expectedxF.add(LUND_pi.Get_H_XF().get(smallindex)) ;
														expectedZ.add(LUND_pi.Get_H_z().get(smallindex));
														expectedPT.add(LUND_pi.Get_H_PT().get(smallindex));
														expectedTheta.add(LUND_pi.Get_H_theta().get(smallindex));
														expectedPhiClas12.add(LUND_pi.Get_H_phi().get(smallindex));
														expectedPhiTrento.add(LUND_pi.Get_H_phiTrento().get(smallindex));
														expectedPid.add(LUND_pi.Get_H_PID().get(smallindex));
														expectedMomentum.add(LUND_pi.Get_H_mom().get(smallindex));
														expectedmx.add(LUND_pi.Get_H_MX().get(smallindex));
														expectedeta.add(LUND_pi.Get_H_etaB().get(smallindex));
														
														//System.out.println(LUND_pi.Get_H_mom().get(smallindex));
													//Histo.SetSkimmedHadron(Pione.getSector(63),New_Pi,Chi2Pid2,H_TrackChi2,H_TackNDF,Pione,electronRec);
												
											
												}
													  	
												
												else {
													ContaiPioniBuoni++; // Count the good pions
													TuplHadron_momentum.add(New_Pi.getMomentum());
													TuplHadron_theta.add(New_Pi.getTheta());
													TuplHadron_phiclas.add(New_Pi.getPhiClas12());
													TuplHadron_z.add(New_Pi.getZ());
													TuplHadron_pt.add(New_Pi.getPt());
													TuplHadron_phi.add(New_Pi.getPhi());
													TuplHadron_xf.add(New_Pi.getXF());
													TuplHadron_MX.add(New_Pi.getMissMass());
													TuplHadron_eta.add(New_Pi.getEtaBar());
															
													TuplHadron_DC1x.add(Pione.getDetectorHits(61).x());
													TuplHadron_DC1y.add(Pione.getDetectorHits(61).y());								
													TuplHadron_DC2x.add(Pione.getDetectorHits(62).x());
													TuplHadron_DC2y.add(Pione.getDetectorHits(62).y());								
													TuplHadron_DC3x.add(Pione.getDetectorHits(63).x());
													TuplHadron_DC3y.add(Pione.getDetectorHits(63).y());
													TuplHadron_DC2sector.add(Pione.getSector(62));
								
													TuplHadron_chi2pid.add(Chi2Pid2);
													expectedxF.add(0.0) ;
													expectedZ.add(0.0);
													expectedPT.add(0.0);
													expectedTheta.add(0.0);
													expectedPhiClas12.add(0.0);
													expectedPhiTrento.add(0.0);
													expectedPid.add(0);
													expectedMomentum.add(0.0);
													expectedmx.add(0.0);
													expectedeta.add(0.0);
												}
													}
												} // Loop over invariant mass and other cuts
											}//  ALL K 
										
									
									


										
										
										
										
										if(TuplHadron_momentum.size()>0) {
											Bank tuple = new Bank(schemaW,TuplHadron_momentum.size());

										for(int kk=0; kk<TuplHadron_momentum.size() ; kk++) {
										tuple.putInt("runnumber", kk, rnmb);
										tuple.putLong("event", kk, evento);
										tuple.putByte("helicity", kk, (byte) helicity);
										tuple.putInt("pid", kk, PionID);
										tuple.putDouble("e_p", kk, electronRec.p());
										tuple.putDouble("e_theta", kk,  electronRec.theta());
										tuple.putDouble("e_phi", kk,  electronRec.phi());
										tuple.putDouble("Q2", kk,  Q2);
										tuple.putDouble("W", kk,  vecW2.mass());
										tuple.putDouble("x", kk,  xB);
										tuple.putDouble("y", kk,  y_var);
										tuple.putDouble("epsilon",kk,  epsilon);
										
										tuple.putDouble("PCAL_lv",kk, electronRec.getDetectorHits(71).y());
										tuple.putDouble("PCAL_lw",kk,  electronRec.getDetectorHits(71).z());
								           
								           
										tuple.putDouble("pi_p", kk, TuplHadron_momentum.get(kk));
										tuple.putDouble("pi_theta", kk, TuplHadron_theta.get(kk));
										tuple.putDouble("pi_phi", kk, TuplHadron_phiclas.get(kk));
										tuple.putDouble("z", kk, TuplHadron_z.get(kk));
										tuple.putDouble("pt", kk,TuplHadron_pt.get(kk));
										tuple.putDouble("phi", kk, TuplHadron_phi.get(kk));
										tuple.putDouble("xf", kk, TuplHadron_xf.get(kk));
										tuple.putDouble("mX", kk, TuplHadron_MX.get(kk));
										tuple.putDouble("eta", kk, TuplHadron_eta.get(kk));	
										tuple.putDouble("pi_chi2pid", kk, TuplHadron_chi2pid.get(kk));	
										tuple.putDouble("DC1_x", kk, TuplHadron_DC1x.get(kk));	
										tuple.putDouble("DC1_y", kk, TuplHadron_DC1y.get(kk));	
										tuple.putDouble("DC2_x", kk, TuplHadron_DC2x.get(kk));	
										tuple.putDouble("DC2_y", kk, TuplHadron_DC2y.get(kk));	
										tuple.putDouble("DC3_x", kk, TuplHadron_DC3x.get(kk));	
										tuple.putDouble("DC3_y", kk, TuplHadron_DC3y.get(kk));	
										
										tuple.putInt("DC2_sector", kk, TuplHadron_DC2sector.get(kk));	
								      
										
										
										tuple.putInt("generated_id", kk,expectedPid.get(kk));
										tuple.putDouble("generated_pi_p", kk,expectedMomentum.get(kk));
										tuple.putDouble("generated_pi_theta", kk,expectedTheta.get(kk));
										tuple.putDouble("getenerated_pi_phi", kk,expectedPhiClas12.get(kk));
										tuple.putDouble("generated_z", kk,expectedZ.get(kk) );
										tuple.putDouble("generated_pt", kk,expectedPT.get(kk));
										tuple.putDouble("generated_phi", kk,expectedPhiTrento.get(kk));
										tuple.putDouble("generated_xf", kk,expectedxF.get(kk) );
										tuple.putDouble("generated_mX", kk,expectedmx.get(kk));
										tuple.putDouble("generated_eta", kk,expectedeta.get(kk));

										tuple.putInt("status", kk, 100);
							
										
										//System.out.println(" Id "+kk+ " mom " + TuplHadron_momentum.get(kk)  +" theta" + TuplHadron_theta.get(kk) +  "phi" + TuplHadron_phiclas.get(kk)+ "z " + TuplHadron_z.get(kk)+ " pt " + TuplHadron_pt.get(kk) + "phitT "+ TuplHadron_phi.get(kk)  );
										}
										
										//	System.out.println("Salvo la tupla, evento numero "+ evento);
										Event  eventoW = new Event();
									    eventoW.reset();
									     eventoW.write(tuple);
								        writeroutput.addEvent(eventoW);	
								      
										}
										else {											
											Bank tuple = new Bank(schemaW,1);
											tuple.putInt("runnumber", 0, rnmb);
											tuple.putLong("event", 0, evento);
											tuple.putByte("helicity", 0, (byte) helicity);
											tuple.putInt("pid", 0, 11);
											tuple.putDouble("e_p", 0, electronRec.p());
											tuple.putDouble("e_theta", 0,  electronRec.theta());
											tuple.putDouble("e_phi", 0,  electronRec.phi());
											tuple.putDouble("Q2", 0,  Q2);
											tuple.putDouble("W", 0,  vecW2.mass());
											tuple.putDouble("x", 0,  xB);
											tuple.putDouble("y", 0,  y_var);
											tuple.putDouble("epsilon",0,  epsilon);
											
											tuple.putDouble("PCAL_lv",0, electronRec.getDetectorHits(71).y());
											tuple.putDouble("PCAL_lw",0,  electronRec.getDetectorHits(71).z());
									           
									           
											tuple.putDouble("pi_p", 0, 0);
											tuple.putDouble("pi_theta", 0, 0);
											tuple.putDouble("pi_phi", 0, 0);
											tuple.putDouble("z", 0, 0);
											tuple.putDouble("pt", 0,0);
											tuple.putDouble("phi", 0, 0);
											tuple.putDouble("xf", 0, 0);
											tuple.putDouble("mX", 0, 0);
											tuple.putDouble("eta", 0, 0);	
											
											tuple.putDouble("pi_chi2pid", 0, 0);	
											tuple.putDouble("DC1_x", 0, 0);	
											tuple.putDouble("DC1_y", 0, 0);	
											tuple.putDouble("DC2_x", 0, 0);	
											tuple.putDouble("DC2_y", 0, 0);	
											tuple.putDouble("DC3_x", 0, 0);	
											tuple.putDouble("DC3_y", 0, 0);	
											
											tuple.putInt("DC2_sector", 0, 0);	
									      			
											tuple.putInt("generated_id", 0,0);
											tuple.putDouble("generated_pi_p", 0,0);
											tuple.putDouble("generated_pi_theta", 0,0);
											tuple.putDouble("getenerated_pi_phi", 0,0);
											tuple.putDouble("generated_z", 0,0);
											tuple.putDouble("generated_pt", 0,0);
											tuple.putDouble("generated_phi", 0,0);
											tuple.putDouble("generated_xf", 0,0 );
											tuple.putDouble("generated_mX", 0,0);
											tuple.putDouble("generated_eta", 0,0);
											
											Event  eventoW = new Event();
										    eventoW.reset();
										     eventoW.write(tuple);
									        writeroutput.addEvent(eventoW);	
										}
										//Histo.ElectronResolutions(LUND_pi.e_p,LUND_pi.e_theta,LUND_pi.e_phi,LUND_pi.e_Q2,LUND_pi.e_xB,LUND_pi.e_W,electronRec,-vecQ2.mass2(),xB,vecW2.mass());

									} //if pi != pi0
									
										
										}// SIDIS CUT 
									} // Good Electric status Fiducial 
								} // DIS electron
							//} // Filter 

						if (debug==true) {
							//if(good_e_count==200000 ) {
							if(good_e_count>=goodeventsStop ) {
								System.out.println("Two photons events with photon cuts " + twophotons);   
								System.out.println("All count"+ count+ " while " + good_e_count + " are FINAL good electrons");	
								System.out.println(all_electrons + " FINAL all electron filtered ");
								SimpleDateFormat dateFormat = new SimpleDateFormat("MMMdd_hh_mm");
								Date now = new Date();
								String time = dateFormat.format(now);
								String time2; 
							       if(MC==true) {
							    	   time2 = now+"MC";
									LUND_pi.SaveHistograms_MultiBins(workdirout,time2); // produce the MC histograms
							    		LUND_pi.getNr("z");
										LUND_pi.Histo(time2,workdirout);

							       }

								else  {time2 = time;}
								//	Histo.Print_ChargedHadron(time2); // this print a comparison Generated- Reconstructed Hadron in MC 

								System.out.println("Printing histos ");		
								if (pi0==true) {
									System.out.println(" --- I am analyzing Neutral Pions  " );
									Analysis.CountPi0s(1,3,XMass,InvariantMassBins,  MissingM, PionCounts,Counts_Phi, Helicity0,Helicity1, Electron_counts,time2, MC,workdirout);
								}
								else {
									System.out.println(" --- I am analyzing Charged Pions  " );
									Analysis.CountPions(XMass, Electron_counts,time2, MC,workdirout);
								}

								System.out.println(" ============================================= ");
								System.out.println(" Events read " + count);
								System.out.println(" Elettroni DIS for normalizaiton " +ContaElDIS);//printing the invariant mass for pi0 
								System.out.println(" Elettroni DIS in good range for normalizaiton " +GoodRangeEl);//printing the invariant mass for pi0 
								System.out.println(" Elettroni DIS " +good_e_count);//printing the invariant mass for pi0 
								System.out.println(" All Electrons " +all_electrons);//printing the invariant mass for pi0
								System.out.println(" Good Pions are" + ContaiPioniBuoni);
						        System.out.println("ContaStep"+ContaStep0); 
								System.out.println("ContaStep1 "+ContaStep1);
						        System.out.println("ContaStep2 "+ContaStep2);
						        System.out.println("ContaStep3 "+ContaStep3);
						        System.out.println("ContaStep4 "+ContaStep4);
						          
								System.out.println(" ============================================= ");
								LUND_pi.getNr("z");
								LUND_pi.Histo(time2,workdirout);
								SIDIS_Cuts.Print(time2,workdirout);
								//Histo.PrintUnSkim_el(time2); Histo.PrintSkim_el(time2);Histo.PrintUnSkim_ph(time2); Histo.PrintSkim_ph(time2);
								System.out.println(" ------ BININNG USED IN Z AND PT ARE ----");
								System.out.println(" Z " + Z_Bins);
								System.out.println(" Pt " + PT_Bins);
								  writeroutput.close();


								return;
							} // if number of events = X
						}   // if debug is true 
						
						}// If I have the rows
					//	}
				//	}
					 
					}// Reader has next file			
				} // A good File 
			
	
			//outP.close();
			//outP2.close();
			//outPel.close();
			//outPip.close();
		//}
		System.out.println("Two photons events with photon cuts " + twophotons);   
		System.out.println("All count"+ count+ " while " + good_e_count + " are FINAL good electrons");	
		System.out.println(all_electrons + " FINAL all electron filtered ");
		SimpleDateFormat dateFormat = new SimpleDateFormat("MMMdd_hh_mm");
		Date now = new Date();
		String time = dateFormat.format(now);
		String time2; 
	       if(MC==true) {
	    	   time2 = now+"MC";
			LUND_pi.SaveHistograms_MultiBins(workdirout,time2); // produce the MC histograms
	    		LUND_pi.getNr("z");
				LUND_pi.Histo(time2,workdirout);

	       }

		else  {time2 = time;}
			//Histo.Print_ChargedHadron(time2); // this print a comparison Generated- Reconstructed Hadron in MC 

		System.out.println("Printing histos ");		
		if (pi0==true) {
			System.out.println(" --- I am analyzing Neutral Pions  " );
			Analysis.CountPi0s(1,3,XMass,InvariantMassBins,  MissingM, PionCounts,Counts_Phi, Helicity0,Helicity1, Electron_counts,time2, MC,workdirout);
		}
		else {
			System.out.println(" --- I am analyzing Charged Pions  " );
			Analysis.CountPions(XMass, Electron_counts,time2, MC,workdirout);
		}

		System.out.println(" ============================================= ");
		System.out.println(" Events read " + count);
		System.out.println(" Elettroni DIS for normalizaiton " +ContaElDIS);//printing the invariant mass for pi0 
		System.out.println(" Elettroni DIS in good range for normalizaiton " +GoodRangeEl);//printing the invariant mass for pi0 
		System.out.println(" Elettroni DIS " +good_e_count);//printing the invariant mass for pi0 
		System.out.println(" All Electrons " +all_electrons);//printing the invariant mass for pi0
		System.out.println(" Good Pions are" + ContaiPioniBuoni);
        System.out.println("ContaStep"+ContaStep0); 
		System.out.println("ContaStep1 "+ContaStep1);
        System.out.println("ContaStep2 "+ContaStep2);
        System.out.println("ContaStep3 "+ContaStep3);
        System.out.println("ContaStep4 "+ContaStep4);
          
		System.out.println(" ============================================= ");
		LUND_pi.getNr("z");
		LUND_pi.Histo(time2,workdirout);
		SIDIS_Cuts.Print(time2,workdirout);
		//Histo.PrintUnSkim_el(time2); Histo.PrintSkim_el(time2);Histo.PrintUnSkim_ph(time2); Histo.PrintSkim_ph(time2);
		System.out.println(" ------ BININNG USED IN Z AND PT ARE ----");
		System.out.println(" Z " + Z_Bins);
		System.out.println(" Pt " + PT_Bins);
		  writeroutput.close();


	}
	
/**
	 * @param fromBank the bank containing the index variable
	 * @param idxVarName the name of the index variable
	 * @return map with keys being the index in toBank and values the indices in fromBank
	 */
	private static  Map<Integer,List<Integer>> loadMapByIndex(	Bank rec_Calorimeter,String idxVarName) {
		Map< Integer,List<Integer> > map = new HashMap <Integer, List<Integer> >();
		if (rec_Calorimeter!=null) {
			for (int iFrom=0; iFrom<rec_Calorimeter.getRows(); iFrom++) {
				//System.out.println(" IFROM + " + iFrom);
				final int iTo = rec_Calorimeter.getInt(idxVarName,iFrom);
				if (!map.containsKey(iTo)) map.put(iTo,new ArrayList<Integer>()); 
				//	System.out.println("iTO e' " + iTo);
				map.get(iTo).add(iFrom);

			}
		}
		else {System.out.println(" Banca e' nulla? Non dovrebbe! ");
		};
		return map;
	}

	private static void InitializeMultiBins(boolean PolygonalBinClas, int bincount_q2,int bincount_x, int bincount_z, int bincount_pt, int binmass, int bincount_phi,
			double[] bin_q2,double[] bin_x, double[] bin_z, double[] bin_pt, double[] bin_mass, double[] bin_phi) {
	
		int Q2_Bins=bincount_q2;
		int Xb_Bins=bincount_x;
		int Z_Bins=bincount_z;
		int PT_Bins=bincount_pt;
		int Phi_Bins=bincount_phi;
		int Mass_Bins=binmass;
		System.out.println("BIn Q2 " + Q2_Bins+ " binq2 min "+ bin_q2[0]);
		InvariantMassBins = new MultiBins(Q2_Bins,Xb_Bins,Z_Bins,PT_Bins,Mass_Bins);
		InvariantMassBins.SetName_1stVariable("Q2");InvariantMassBins.SetName_2ndVariable("Xb");
		InvariantMassBins.SetName_3rdVariable("z");InvariantMassBins.SetName_4thVariable("Pt");
		InvariantMassBins.SetName_5thVariable("Invariant Mass");
		System.out.println( " Multidimensional Analysis code ");
		System.out.println(" Bin1: " + InvariantMassBins.GetName_1stVariable() + " Bin2: " + InvariantMassBins.GetName_2ndVariable());
		System.out.println(" Bin3: " + InvariantMassBins.GetName_3rdVariable() + " Bin4: " + InvariantMassBins.GetName_4thVariable());
		InvariantMassBins.SetUnevenBins(bin_q2,bin_x,bin_z,bin_pt,bin_mass);
		InvariantMassBins.GenerateHistograms("InvariantMass");
		if(PolygonalBinClas==true) {
			InvariantMassBins.InizializeClasPolygons();
		}
		PionCounts = new MultiBins(Q2_Bins,Xb_Bins,Z_Bins,PT_Bins,0);
		PionCounts.SetName_1stVariable("Q2");PionCounts.SetName_2ndVariable("Xb");
		PionCounts.SetName_3rdVariable("z");PionCounts.SetName_4thVariable("Pt");	 System.out.println( " Multidimensional Analysis code ");
		System.out.println(" Bin1: " + PionCounts.GetName_1stVariable() + " Bin2: " + PionCounts.GetName_2ndVariable());
		System.out.println(" Bin3: " + PionCounts.GetName_3rdVariable() + " Bin4: " + PionCounts.GetName_4thVariable());
		PionCounts.SetUnevenBins(bin_q2,bin_x,bin_z,bin_pt,bin_mass);
		PionCounts.GenerateHistograms("PionCounts");
		if(PolygonalBinClas==true) { PionCounts.InizializeClasPolygons();}
		// Counting particle in Phi Bins 
		Counts_Phi = new MultiBins(Q2_Bins,Xb_Bins,Z_Bins,PT_Bins,Phi_Bins);
		Counts_Phi.SetName_1stVariable("Q2");Counts_Phi.SetName_2ndVariable("Xb");
		Counts_Phi.SetName_3rdVariable("z");Counts_Phi.SetName_4thVariable("Pt"); Counts_Phi.SetName_5thVariable("Phi");
		Counts_Phi.SetUnevenBins(bin_q2,bin_x,bin_z,bin_pt,bin_phi);
		Counts_Phi.GenerateHistograms("Counts_Phi");
		if(PolygonalBinClas==true) { Counts_Phi.InizializeClasPolygons();}

		// Counting particle with Helicity0 as function of Phi bins. 
		Helicity0 = new MultiBins(Q2_Bins,Xb_Bins,Z_Bins,PT_Bins,Phi_Bins);
		Helicity1 = new MultiBins(Q2_Bins,Xb_Bins,Z_Bins,PT_Bins,Phi_Bins);
		Helicity0.SetName_1stVariable("Q2");Helicity0.SetName_2ndVariable("Xb");
		Helicity0.SetName_3rdVariable("z");Helicity0.SetName_4thVariable("Pt");
		Helicity1.SetName_1stVariable("Q2");Helicity1.SetName_2ndVariable("Xb");
		Helicity1.SetName_3rdVariable("z");Helicity1.SetName_4thVariable("Pt");
		Helicity0.SetName_5thVariable("Phi");Helicity1.SetName_5thVariable("Phi");
		Helicity0.SetUnevenBins(bin_q2,bin_x,bin_z,bin_pt,bin_phi);
		Helicity1.SetUnevenBins(bin_q2,bin_x,bin_z,bin_pt,bin_phi);
		Helicity0.GenerateHistograms("Helicity0"); Helicity1.GenerateHistograms("Helicity1");
		if(PolygonalBinClas==true) {Helicity0.InizializeClasPolygons(); Helicity1.InizializeClasPolygons();}

		// Counting particle with MissinM as function of Phi bins. 
		//it was 200 bins now it is just 3 because not used 
		double []binMissinM= new double[300];
		for(int i=0; i<300; i++) {
			double binsize =3.0/299 ;
			double minvalue=0.0;
			double conto=minvalue+binsize*i;
			binMissinM[i]=conto;
		}
		MissingM = new MultiBins(Q2_Bins,Xb_Bins,Z_Bins,PT_Bins,300 );
		MissingM.SetName_1stVariable("Q2");Helicity0.SetName_2ndVariable("Xb");
		MissingM.SetName_3rdVariable("z");Helicity0.SetName_4thVariable("Pt");
		MissingM.SetName_5thVariable("XMass");Helicity1.SetName_5thVariable("XMass");
		MissingM.SetUnevenBins(bin_q2,bin_x,bin_z,bin_pt,  binMissinM);
		MissingM.GenerateHistograms("MissingM");
		if(PolygonalBinClas==true) {MissingM.InizializeClasPolygons();}
		// Creating MultiDimensional Histograms for MC Analysis of Given Particle 

		MCParticles = new MultiBins(Q2_Bins,Xb_Bins,Z_Bins,PT_Bins,0);
		MCParticles.SetName_1stVariable("Q2");MCParticles.SetName_2ndVariable("Xb");
		MCParticles.SetName_3rdVariable("z");MCParticles.SetName_4thVariable("Pt");
		MCParticles.SetUnevenBins(bin_q2,bin_x,bin_z,bin_pt,bin_phi );
		MCParticles.GenerateHistograms("MCParticles");
		if(PolygonalBinClas==true) {MCParticles.InizializeClasPolygons();}
		// Creating an H2F for the information of the electron 

		Electron_counts = new MultiBins(Q2_Bins,Xb_Bins,0,0,0);
		Electron_counts.SetName_1stVariable("Q2");MCParticles.SetName_2ndVariable("Xb");
		Electron_counts.SetUnevenBins(bin_q2,bin_x, new double[]{0,0}, new double[]{0,0}, new double[]{0,0} );
		Electron_counts.GenerateHistograms("Electron_counts");
		System.out.println(Electron_counts.FstBin+" max "+ Electron_counts.GetUnevenBins(1)[0]+" "+ Electron_counts.GetUnevenBins(1)[1]);
		if(PolygonalBinClas==true) { Electron_counts.SetPolygonal(PolygonalBinClas); // it will compute properly the borders}
		}

	}
}
