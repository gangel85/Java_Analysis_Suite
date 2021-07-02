package PHD_Clean;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;

import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.data.H3F;
import org.jlab.groot.data.TDirectory;
import org.jlab.groot.fitter.DataFitter;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.math.F1D;
import org.jlab.groot.ui.TGCanvas;
import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.data.Schema;
import org.jlab.jnp.hipo4.data.Schema.SchemaBuilder;
import org.jlab.jnp.hipo4.io.HipoWriter;
import org.jlab.jnp.physics.LorentzVector;


public class AnalysisSIDIS {
	private double pi_peak=0.135; // this quantity is used to find the center of the peak 
	private double EC_res=0.012; //sigma of the pi0 peak in the calorimeter
private double Polarization = 0.863;
	public static H1F Phi_Dist_Overall,Z_Dist_Overall,XB_Dist_Overall,Pt_Dist_Overall, BSA_Phi;
	static List<H1F> InvMassPlots = new ArrayList<H1F>();
	static List<H1F> N0_Phi = new ArrayList<H1F>();
	static List<H1F> N1_Phi = new ArrayList<H1F>();
	static List<H1F> Pt_Dist = new ArrayList<H1F>();
	static List<H1F> Z_Dist = new ArrayList<H1F>();
	static List<H1F> XB_Dist = new ArrayList<H1F>();


	private H1F H_MissingMass = new H1F("H_MissingMass",200,0,3);

	private int zbins =0;
	private int PTbins=0;
	private int PT_Bins=0;
	private int Z_Bins=0;
	private int PhiBins=10;
	private double Phi_Min=-180, Phi_Max=180;
	private H3F InvMass;
	private H3F N0_Count;
	private H3F N1_Count;
	private boolean SingleHadron;


	private MultiBins Dist_Q2,Dist_xB,Dist_W,Dist_y,Dist_eps,Dist_Pt_ALL,Dist_z_ALL;
	private MultiBins Dist_zPT_Q2_ALL,Dist_zPT_xB_ALL,Dist_zPT_W_ALL,Dist_zPT_y_ALL,Dist_zPT_eps_ALL,Dist_zPT_PTZQ_ALL,Dist_zPT_EtaB_ALL,Dist_zPT_XF_ALL;
	private MultiBins Dist_ZPT;
	private  Schema schema;
	private HipoWriter writeroutput;
	private Bank tuple;
	
	public AnalysisSIDIS(String workoutdir,boolean PolygonalBinClas, int Z_bins, double Z_min, double Z_max, int PT_bins, double PT_min, double PT_max, boolean HadronBoolean)
	{
		
		

		this.SingleHadron=HadronBoolean;
		this.Z_Bins=Z_bins;
		this.PT_Bins=PT_bins;
		//InvMass = new H3F(Z_bins, Z_min, Z_max, PT_bins, PT_min, PT_max,Mass_Bins,min_pion,max_pion);
		N0_Count = new H3F(Z_bins, Z_min,Z_max, PT_bins, PT_min, PT_max,12,0,360);
		N1_Count = new H3F(Z_bins, Z_min, Z_max, PT_bins, PT_min, PT_max,12,0,360);	
		//	Plots();
		Z_Dist_Overall = new H1F("Z_Dist_Overall","Z_Dist_Overall",100,0,1);
		Z_Dist_Overall.setTitle("Z distribution integrated");
		Z_Dist_Overall.setTitleX("Z");

		XB_Dist_Overall = new H1F("Xb_Dist_Overall","Xb_Dist_Overall",100,0,1);
		XB_Dist_Overall.setTitle("Xb distribution integrated in Z");
		XB_Dist_Overall.setTitleX("Xb");

		Pt_Dist_Overall = new H1F("Pt_Dist_Overall","Pt_Dist_Overall",100,0,2);
		Pt_Dist_Overall.setTitle("Pt distribution integrated in Z");
		Pt_Dist_Overall.setTitleX("Pt ");

		Phi_Dist_Overall = new H1F("Phi_Dist_Overall","Phi_Dist_Overall",12,0,360);
		Phi_Dist_Overall.setTitle("Phi distribution integrated in Z");
		Phi_Dist_Overall.setTitleX("Phi");

		H_MissingMass.setTitle("H_MissingMass");
		H_MissingMass.setTitleX("Missing Mass [GeV]");

	}

	public void Analyze(boolean PolygonalBinClas, LorentzVector beam, LorentzVector target, ParticleREC electronRec, int helicity, H1F xMass, MultiBins invariantMassBins, MultiBins missingM, int sig, double xB, double Q2, Pi0_Particle new_Pi0s, MultiBins PionCounts, MultiBins Counts_Phi, MultiBins H0_Counts, MultiBins H1_Counts , MultiBins Electron_counts, double W, double y, double eps, long evento) {
		
			

		if(this.SingleHadron==true) {
			int k= new_Pi0s.GetIndex_MostEnergy();

			Z_Dist_Overall.fill(new_Pi0s.getZ().get(k));
			XB_Dist_Overall.fill(xB);
			Pt_Dist_Overall.fill(new_Pi0s.getPt().get(k));

			
			LorentzVector MissM= new LorentzVector();
			LorentzVector electronLV = new LorentzVector();
			electronLV.setPxPyPzE(electronRec.px(), electronRec.py(),electronRec.pz(), electronRec.e());
			MissM.add(beam).add(target).sub(electronLV).sub(new_Pi0s.getLorentzVector().get(k));
			//double MissingMass= MissM.mass();
			//xMass.fill(MissingMass);
			int Q2bin,XBbin;
			if(PolygonalBinClas==true) {
				 Q2bin= PionCounts.ComputePolygonalBins(Q2, xB).get(0);
				 XBbin = PionCounts.ComputePolygonalBins(Q2, xB).get(1);
				}
				else {
				 Q2bin = PionCounts.ComputeBin(1, Q2);
				 XBbin = PionCounts.ComputeBin(2, xB);
				}
			

			//missingM.GetH3F(Q2bin, XBbin).fill(new_Pi0s.getZ().get(k), new_Pi0s.getPt().get(k), MissM.mass());
	
				if(Q2bin>=0 && XBbin>=0 ) {
					invariantMassBins.GetH3F(Q2bin,XBbin).fill(new_Pi0s.getZ().get(k), new_Pi0s.getPt().get(k), new_Pi0s.getMass().get(k));
					if(helicity==-1)	{
						N0_Count.fill(new_Pi0s.getZ().get(k), new_Pi0s.getPt().get(k), new_Pi0s.getPhi().get(k)); 
						if (new_Pi0s.getMass().get(k) <(pi_peak+sig*EC_res) && new_Pi0s.getMass().get(k)> (pi_peak-sig*EC_res)) {
							H0_Counts.GetH3F(Q2bin, XBbin).fill(new_Pi0s.getZ().get(k), new_Pi0s.getPt().get(k), new_Pi0s.getPhi().get(k));
							
						}
					}
					else if(helicity==1) 
					{
						N1_Count.fill(new_Pi0s.getZ().get(k), new_Pi0s.getPt().get(k), new_Pi0s.getPhi().get(k));
						if (new_Pi0s.getMass().get(k) <(pi_peak+sig*EC_res) && new_Pi0s.getMass().get(k)> (pi_peak-sig*EC_res)) {
							H1_Counts.GetH3F(Q2bin, XBbin).fill(new_Pi0s.getZ().get(k), new_Pi0s.getPt().get(k), new_Pi0s.getPhi().get(k));
						}
					}
				}
		}
		else if (this.SingleHadron==false){	
			
			for(int kk=0;kk<new_Pi0s.getNrPi0s(); kk++) {			
			if (new_Pi0s.getMass().get(kk) <(pi_peak+sig*EC_res) && new_Pi0s.getMass().get(kk)> (pi_peak-sig*EC_res)) {
//System.out.println(" I am here ");
			Z_Dist_Overall.fill(new_Pi0s.getZ().get(kk));
			XB_Dist_Overall.fill(xB);
			Pt_Dist_Overall.fill(new_Pi0s.getPt().get(kk));			
			
		int Q2bin,XBbin;
	if(PolygonalBinClas==true) {}
			else {}
			

	
			}//Within X Sigmas from peak			
		}//kk index of pions
			
			


		}
	}

		

		
		





	public void Plots()
	{


		BSA_Phi = new H1F("BSA_Phi","BSA_Phi",12,0,360);
		BSA_Phi.setTitle("BSA Phi Dist");
		BSA_Phi.setTitleX("Phi [rad]");

		Z_Dist_Overall = new H1F("Z_Dist_Overall","Z_Dist_Overall",100,0,1);
		Z_Dist_Overall.setTitle("Z distribution integrated");
		Z_Dist_Overall.setTitleX("Z");

		XB_Dist_Overall = new H1F("Xb_Dist_Overall","Xb_Dist_Overall",100,0,1);
		XB_Dist_Overall.setTitle("Xb distribution integrated in Z");
		XB_Dist_Overall.setTitleX("Xb");

		Pt_Dist_Overall = new H1F("Pt_Dist_Overall","Pt_Dist_Overall",100,0,2);
		Pt_Dist_Overall.setTitle("Pt^2 distribution integrated in Z");
		Pt_Dist_Overall.setTitleX("Pt^2 ");

		Phi_Dist_Overall = new H1F("Phi_Dist_Overall","Phi_Dist_Overall",12,0,360);
		Phi_Dist_Overall.setTitle("Phi distribution integrated in Z");
		Phi_Dist_Overall.setTitleX("Phi");


	}

	
	public void CountPi0s(int fit, int sig, H1F xMass,MultiBins invariantMassBins, MultiBins MissingM, MultiBins PionCounts,MultiBins Counts_Phi, MultiBins H0_Counts, MultiBins H1_Counts, MultiBins Electron_counts, String time, boolean mC, String workdirout) throws FileNotFoundException {
		System.out.println("*************************");
		System.out.println(" Counting the Neutral Pion ");
		System.out.println("*************************");
		List<F1D> functions = new ArrayList<F1D>(); 
		double totEntry=0;

		int Mass_Bins=invariantMassBins.FthBin;
		double min_pion=invariantMassBins.Fth_Min;
		double max_pion= invariantMassBins.Fth_Max;


		int Phi_Bins= H0_Counts.FthBin;
		double min_phi = 0.0;
		double max_phi = 360.0;

		System.out.println(" Phi_Bins " + Phi_Bins);
		System.out.println(" Min Phi " + min_phi + " Max Phi "+ max_phi);

	

		String s1=workdirout;
		new File(s1+"/Plots").mkdir();
		TDirectory directoryMain = new TDirectory();	
		String PathFolder = s1+"/Plots/"+time;
		File dir = new File(PathFolder);dir.mkdir();
		PathFolder=s1+"/Plots/"+time+"/MultiBins";
		new File (PathFolder).mkdir();
		System.out.println(" (Pi) - > I am Plotting Multi Bins Histograms");
		// Looping over all bins and fill the 3Dimensional Invariant Mass bins (z,pt,mass) by using the 5dimensional information of InvariantMassBins
		for(int x=0 ; x <(invariantMassBins.FstBin);x++) { //Q2bins
			String number = Integer.toString(x);
			PathFolder=s1+"/Plots/"+time+"/MultiBins/"+invariantMassBins.GetName_1stVariable()+ "_"+number;
			File dirbin = new File(PathFolder);dirbin.mkdir();
			//	directory.mkdir(PathFolder); directory.cd("PathFolder");
			for (int y=0 ; y<(invariantMassBins.SndBin); y++) { //Xb bins

				
			

				
				
			
				
				
				
				
				
				
				
				
				
				

				String number2=Integer.toString(y);
				PathFolder=s1+"/Plots/"+time+"/MultiBins/"+invariantMassBins.GetName_1stVariable()+ "_"+number+"/"+invariantMassBins.GetName_2ndVariable()+ "_"+number2;
				File dirbin2 = new File(PathFolder);dirbin2.mkdir();
				File File_txt = new File(PathFolder+"/results.txt");
				directoryMain.mkdir("/main/"); directoryMain.cd("/main/");

				//directoryMain.addDataSet(GaussianCounts.GetH2F(, k));

				try {
					File_txt.createNewFile();
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				if(File_txt.exists()) {
		
				
					int myindex=0;
					int badZ = 0;
					int badPT= 0; 
					boolean badFit=false;
			
					System.out.println("----- ");
					directoryMain.addDataSet(Electron_counts.GetH2());
					//System.out.println("The pion Energy entries are " +Pion_Energy.GetH1F(x, y).getEntries());
					directoryMain.writeFile(PathFolder+"/HistogramsCounts.hipo");
				} //File exitst

			} //Xb bins 
		} //Q2 bins


		return;
	}


	public void CountPions( H1F xMass, MultiBins Electron_counts, String time, boolean mC, String workdirout) throws FileNotFoundException {

	
      
		
		
		System.out.println(" %%%%%%%%%%%%%%%%%%%%%%%%%%%%");
		System.out.println(" %%%%%% " + Electron_counts.Fst_Min+"%%%%%%%");
		System.out.println(" %%%%%%%%%%%%%%%%%%%%%%%%%%%%");
		List<F1D> functions = new ArrayList<F1D>(); 


	

		String s1=workdirout;
		new File(s1+"/Plots").mkdir();
		// TDirectory was here 
		
		String PathFolder = s1+"/Plots/"+time;
		File dir = new File(PathFolder);dir.mkdir();
		PathFolder=s1+"/Plots/"+time+"/MultiBins";
		new File (PathFolder).mkdir();
		System.out.println(" (Pi) - > I am Plotting Multi Bins Histograms");
		// Looping over all bins and fill the 3Dimensional Invariant Mass bins (z,pt,mass) by using the 5dimensional information of InvariantMassBins
		for(int x=0 ; x <(Electron_counts.FstBin);x++) { //Q2bins
			String number = Integer.toString(x);
			PathFolder=s1+"/Plots/"+time+"/MultiBins/"+Electron_counts.GetName_1stVariable()+ "_"+number;	
			File dirbin = new File(PathFolder);dirbin.mkdir();
			for (int y=0 ; y<(Electron_counts.SndBin); y++) { //Xb bins
				
				TDirectory directoryMain = new TDirectory();
				


	
				

				
				
				
				
				String number2=Integer.toString(y);
				PathFolder=s1+"/Plots/"+time+"/MultiBins/"+Electron_counts.GetName_1stVariable()+ "_"+number+"/"+Electron_counts.GetName_2ndVariable()+ "_"+number2;
				File dirbin2 = new File(PathFolder);dirbin2.mkdir();
				File File_txt = new File(PathFolder+"/results.txt");
				directoryMain.mkdir("/main/"); directoryMain.cd("/main/");
				try {
					File_txt.createNewFile();
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				if(File_txt.exists()) {
					System.out.println(" Output File Found" );
					PrintWriter outP = new PrintWriter(File_txt);
					int myindex=0;

			
					//It was here
					System.out.println("----- ");
				
					directoryMain.addDataSet(Electron_counts.GetH2());
		
					directoryMain.writeFile(PathFolder+"/HistogramsCounts.hipo");
					outP.close();						    // Printing stuff: 	
				} //File exitst



			} //Xb bins 
		} //Q2 bins


		return;
	}





	/**
	 * Save multi dimensional BSA
	 */
	public void Save_BSA() {


	}
	public void Histo (String time)
	{

		String s1="/Users/gangelini/Desktop/";
		new File(s1+"/Plots").mkdir();
		String PathFolder = s1+"/Plots/"+time;
		File dir = new File(PathFolder);
		dir.mkdir();
		PathFolder=s1+"/Plots/"+time+"/MultiBins";
		new File (PathFolder).mkdir();
		//new File("/Users/gangelini/Desktop/Plots/Inv_Mass").mkdir();


		System.out.println(" (Pi) - > I am Plotting Invariant Mass Histograms");

		EmbeddedCanvas photons = new EmbeddedCanvas();
		photons.setSize(1600,1000);
		photons.divide(PT_Bins,Z_Bins);
		photons.setAxisTitleSize(20);
		photons.setAxisFontSize(20);
		photons.setTitleSize(20);
		// I plot z from  from 0.1 to 0.9 because first and last bin are useless
		for(int i=0 ; i<(Z_Bins*PT_Bins); i++) {
			photons.cd(i);photons.draw(InvMassPlots.get(i));
		}
		for (int ii=0; ii<Z_Bins;ii++) {
			System.out.println(" Z bins:");
			System.out.println(this.InvMass.getXAxis().getBinCenter(ii));
		}
		for (int jj=0; jj<PT_Bins;jj++) {
			System.out.println(" PT bins:");
			System.out.println(this.InvMass.getYAxis().getBinCenter(jj));
		}


		String strg0 = String.format("%s/Pi0s.png",PathFolder);
		System.out.println("Saving plots in "+ PathFolder);
		photons.save(strg0);	

		PathFolder=s1+"/Plots/"+time+"/Distribution";
		new File (PathFolder).mkdir();
		//new File("/Users/gangelini/Desktop/Plots/Inv_Mass").mkdir();


		System.out.println(" (Pi) - > I am Plotting Invariant Mass Histograms");

		EmbeddedCanvas distribution = new EmbeddedCanvas();
		distribution.setSize(1600,1000);
		distribution.divide(3,1);
		distribution.setAxisTitleSize(20);
		distribution.setAxisFontSize(20);
		distribution.setTitleSize(20);
		// I plot z from  from 0.1 to 0.9 because first and last bin are useless

		distribution.cd(0);distribution.draw(Z_Dist_Overall);
		distribution.cd(1);distribution.draw(Pt_Dist_Overall);
		distribution.cd(2);distribution.draw(XB_Dist_Overall);
		String strg1= String.format("%s/Distributions.png",PathFolder);
		System.out.println("Saving plots in "+ PathFolder);
		distribution.save(strg1);	



		
}

}
