package PHD_Clean;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;

import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.data.H3F;
import org.jlab.groot.ui.TGCanvas;
import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.data.SchemaFactory;
import org.jlab.jnp.hipo4.io.HipoReader;



public class An_ReadMCALL {
	private static MultiBins PionCountsData,PionCounts,PionsGenerated,PionsGeneratedReco, EtaDist,PtDist,zDist,yDist,XDist,Q2Dist,WDist;

	private static H2F Electron_counts,Electron_countsData;

	
//	private static double correctionfactor= 0.86; //For negative
	
	private static double correctionfactor= 0.87; 

	
	public static void main(String[] args) throws FileNotFoundException {
		 Electron_counts = new H2F("Electron_counts","Electron_counts",6,0.5,6.5,7,0.5,7.5);
		 Electron_countsData = new H2F("Electron_countsData","Electron_countsData",6,0.5,6.5,7,0.5,7.5);
		// TODO Auto-generated method stub
		 DecimalFormat df = new DecimalFormat();
		df.setMaximumFractionDigits(3);
		
		double MinP	=Double.parseDouble(args[0]);
		double MaxP=Double.parseDouble(args[1]);
		double MX =Double.parseDouble(args[2]);
		double MinPGen=Double.parseDouble(args[3]);
		double MaxPGen=Double.parseDouble(args[4]);
		double MXGen = Double.parseDouble(args[5]);
		String workdir1 = args[6];
		String workdir2 = args[7];
		String workdirD =	args[8];
		String Outputdir= args[9];
		
				
		/*
		double MinP = 1.25;
		double MaxP = 3;
		double MX = 1.6;
		
		double MinPGen= 0;
		double MaxPGen = 8;
		double MXGen = 1.6;
	
		
		 
		String workdir1 = "/Volumes/My Passport for Mac/MCProduction/";
		String workdir2 = "/Volumes/My Passport for Mac/GeneratedMC/";
		String workdirD = "/Volumes/My Passport for Mac/DataProduction/";

		String Outputdir= "/Users/gangelini/work/ResultsMultiplicity";
		*/

		boolean debug =false;


		int binmass=300;
		double []bin_mass= new double[300];
		for(int i=0; i<300; i++) {
			double binsize =0.25/299 ;
			double minvalue=0.05;
			double conto=minvalue+binsize*i;
			bin_mass[i]=conto;
		}
		int bincount_phi = 12;
		double bin_phi[] = {0.0, 30.0, 60.0, 90.0, 120.0, 150.0, 180.0, 210.0, 240.0, 270.0, 300.0, 330.0, 360.0};
		int bincount_q2 = 6;
		double bin_q2[] = {1,12};
		int bincount_x = 7;
		double bin_x[] = {0,1};
		int bincount_z = 12;
		double bin_z[] = {0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.75,0.80};
		int bincount_pt =45; 
		double bin_pt[] = {0.0,0.03,0.06,0.09,0.12,0.15,0.18,0.21,0.24,0.27,0.30,0.33,0.36,0.39,0.42,0.45,0.48,0.51,0.54,0.57,0.60,
							0.63,0.66,0.69,0.72,0.75,0.78,0.81,0.84,0.87,0.90,0.93,0.96,0.99,1.02,1.05,1.08,1.11,1.14,1.17,1.20,
							1.23,1.26,1.29,1.32,1.35};
				
		
		InitializeMultiBins( true,  bincount_q2, bincount_x,  bincount_z,  bincount_pt,  binmass,  bincount_phi,
				bin_q2,bin_x,  bin_z,  bin_pt, bin_mass,  bin_phi);
		
		
		ArrayList<H2F> B3_PT = new ArrayList<H2F>();
		ArrayList<H2F> B3_z = new ArrayList<H2F>();




		
		File folder1 = new File(workdir1); File[] listOfFiles1 = folder1.listFiles();
		File folder2 = new File(workdir2); File[] listOfFiles2 = folder2.listFiles();
		File folderD = new File(workdirD); File[] listOfFilesD = folderD.listFiles();


		for (File file : listOfFilesD) {
			HipoReader readerD = new HipoReader(); 

			System.out.println(" Reading  " + file.getName());
			readerD.open(workdirD+file.getName());
			Event event=new Event();
			SchemaFactory schema = readerD.getSchemaFactory();
			Bank mia = new Bank(schema.getSchema("SIDIS::ePiX"));
			
			long eventonrD = 0;
			while(readerD.hasNext()==true){
				eventonrD++;
				if(eventonrD%200000 == 0) System.out.println(eventonrD + " events");

			readerD.nextEvent(event);
			event.read(mia); 
			//mia.show();

			for(int i=0;i<mia.getRows(); i++) {

				double e_p = mia.getDouble("e_p", i);
				double e_theta = Math.toDegrees(mia.getDouble("e_theta", i));	
				double e_phi = Math.toDegrees(mia.getDouble("e_phi", i));
				double Q2 = mia.getDouble("Q2", i);
				double x = mia.getDouble("x", i);
				double W = mia.getDouble("W", i);	
				double eta = mia.getDouble("eta", i);
				double y = mia.getDouble("y", i);

				
				long evento =	mia.getLong("event", i);
				double pion_p = mia.getDouble("pi_p", i);
			//	if(pion_p==0) mia.show();
				double pion_theta = mia.getDouble("pi_theta", i);
				double pion_phi = mia.getDouble("pi_phi", i);
				double pion_z = mia.getDouble("z", i);
				double pion_pt = mia.getDouble("pt", i);
				double pion_phiT = mia.getDouble("phi", i);
				double pion_xF = mia.getDouble("xf", i);
				double pion_mx = mia.getDouble("mX", i);
				
				
				int BinQ2= PionCountsData.ComputePolygonalBins(Q2, x).get(0);
	    		int BinXB = PionCountsData.ComputePolygonalBins(Q2, x).get(1);
	    		Electron_countsData.fill(BinQ2+1, BinXB+1);

	    		if(pion_xF>0.0 &&pion_p>MinP&& pion_p<MaxP && pion_mx>MX ) {	
	    		PionCountsData.GetH2F(BinQ2,BinXB).fill(pion_z, pion_pt);	
	    	 	}
				}
			}
			if (debug==true &&eventonrD==1000 )break;
		}
		
	
		
	
		for (File file : listOfFiles2) {
			HipoReader reader2 = new HipoReader(); 

			System.out.println(" Reading  " + file.getName());
			reader2.open(workdir2+file.getName());
			Event event2=new Event();
			SchemaFactory schema2 = reader2.getSchemaFactory();
			Bank mia2 = new Bank(schema2.getSchema("LUND::ePiX"));

		 
		long eventonr2 = 0; 
		 int conta2=0;
		while(reader2.hasNext()==true) {
				eventonr2++;
				if(eventonr2%200000 == 0) System.out.println(eventonr2 + " events");

			reader2.nextEvent(event2);
			event2.read(mia2);
			for(int i=0;i<mia2.getRows(); i++) {
				
				long evento=	mia2.getLong("event", i);
				int pid= mia2.getInt("pid", i);
				double e_p = mia2.getDouble("e_p", i);
				double e_theta = mia2.getDouble("e_theta", i);	
				double e_phi = mia2.getDouble("e_phi", i);
				double eta = mia2.getDouble("eta", i);
				double y = mia2.getDouble("y", i);

				double Q2 = mia2.getDouble("Q2", i);
				double x = mia2.getDouble("x", i);
				double W = mia2.getDouble("W", i);	
				double pion_p = mia2.getDouble("pi_p", i);
				double pion_theta = mia2.getDouble("pi_theta", i);
				double pion_phi = mia2.getDouble("pi_phi", i);
				double pion_z = mia2.getDouble("z", i);
				double pion_pt = mia2.getDouble("pt", i);
				double pion_phiT = mia2.getDouble("phi", i);
				double pion_xF = mia2.getDouble("xf", i);
				double pion_mx = mia2.getDouble("mX", i);
			
				
				//if (pion_mx < 1.2) System.out.println(" pion missing mass "+ pion_mx);
				
				int BinQ2= PionsGenerated.ComputePolygonalBins(Q2, x).get(0);
	    		int BinXB = PionsGenerated.ComputePolygonalBins(Q2, x).get(1);
				if(pion_xF>0.0 &&pion_p>MinPGen&& pion_p<MaxPGen && pion_mx>MXGen ) {
					
				PionsGenerated.GetH2F(BinQ2,BinXB).fill(pion_z, pion_pt);
				}
				 if (evento>eventonr2) eventonr2=evento;
				 conta2++;
			}
			

		//	mia2.show();
		
			
		}
		if (debug==true &&eventonr2==1000 )break;
		}
		System.out.println(" NR 2  FINITO" );

		
		for (File file : listOfFiles1) {
			HipoReader reader1 = new HipoReader(); 
			System.out.println(" Reading  " + file.getName());
			reader1.open(workdir1+file.getName());
			Event event=new Event();
			SchemaFactory schema = reader1.getSchemaFactory();
			Bank mia = new Bank(schema.getSchema("SIDIS::ePiX"));
			long eventonr1 = 0;int conta1=0;
			while(reader1.hasNext()==true){
				eventonr1++;
				if(eventonr1%200000 == 0) System.out.println(eventonr1 + " events");
			reader1.nextEvent(event);
			event.read(mia); 
			//mia.show();

			for(int i=0;i<mia.getRows(); i++) {
				double e_p = mia.getDouble("e_p", i);
				double e_theta = Math.toDegrees(mia.getDouble("e_theta", i));	
				double e_phi = Math.toDegrees(mia.getDouble("e_phi", i));
				double Q2 = mia.getDouble("Q2", i);
				double x = mia.getDouble("x", i);
				double W = mia.getDouble("W", i);	
				long evento =	mia.getLong("event", i);
				double pion_p = mia.getDouble("pi_p", i);
				double pion_theta = mia.getDouble("pi_theta", i);
				double pion_phi = mia.getDouble("pi_phi", i);
				double pion_z = mia.getDouble("z", i);
				double pion_pt = mia.getDouble("pt", i);
				double pion_phiT = mia.getDouble("phi", i);
				double pion_xF = mia.getDouble("xf", i);
				double pion_mx = mia.getDouble("mX", i);
				double genpion_z = mia.getDouble("generated_z", i);
				double genpion_pt = mia.getDouble("generated_pt", i);
				double eta = mia.getDouble("eta", i);
				double y = mia.getDouble("y", i);
				
				H2F originalval = new H2F("originalval",bincount_z,0.2,0.8,bincount_pt,0,1.35);
				H2F generval = new H2F("generval",bincount_z,0.2,0.8,bincount_pt,0,1.35);
				
				int recbin=	originalval.findBin(pion_z, pion_pt);
				int genbin=	generval.findBin(genpion_z, genpion_pt);
				
				int BinQ2= PionCounts.ComputePolygonalBins(Q2, x).get(0);
	    		int BinXB = PionCounts.ComputePolygonalBins(Q2, x).get(1);
	    		Electron_counts.fill(BinQ2+1, BinXB+1);

	    		if(pion_xF>0.0 &&pion_p>MinP&& pion_p<MaxP && pion_mx>MX ) {
	    			EtaDist.GetH3F(BinQ2, BinXB).fill(pion_z, pion_pt,eta );
	    			PtDist.GetH3F(BinQ2, BinXB).fill(pion_z, pion_pt,genpion_pt );
	    			zDist.GetH3F(BinQ2, BinXB).fill(pion_z, pion_pt,genpion_z );
	    			yDist.GetH3F(BinQ2, BinXB).fill(pion_z, pion_pt,y );
	    			XDist.GetH3F(BinQ2, BinXB).fill(pion_z, pion_pt,x );
	    			Q2Dist.GetH3F(BinQ2, BinXB).fill(pion_z, pion_pt,Q2 );
	    			WDist.GetH3F(BinQ2, BinXB).fill(pion_z, pion_pt,W);
	    			
	    		PionCounts.GetH2F(BinQ2,BinXB).fill(pion_z, pion_pt);
	    		if(recbin==genbin)PionsGeneratedReco.GetH2F(BinQ2,BinXB).fill(pion_z, pion_pt);

			
	    		}
				}
			}
			if (debug==true &&eventonr1==1000 )break;

			
		}

		
	System.out.println( " Elettroni in 0-0 " + Electron_counts.getBinContent(0, 0) + " elettroni in 1,1 "+ Electron_counts.getBinContent(1, 1));
	AnalyzeBins(0,0,Outputdir,bincount_z,bincount_pt);
	AnalyzeBins(1,0,Outputdir,bincount_z,bincount_pt);
	AnalyzeBins(1,1,Outputdir,bincount_z,bincount_pt);
	AnalyzeBins(1,2,Outputdir,bincount_z,bincount_pt);
	AnalyzeBins(2,1,Outputdir,bincount_z,bincount_pt);
	AnalyzeBins(2,2,Outputdir,bincount_z,bincount_pt);
	AnalyzeBins(2,3,Outputdir,bincount_z,bincount_pt);
	AnalyzeBins(2,4,Outputdir,bincount_z,bincount_pt);
	AnalyzeBins(3,2,Outputdir,bincount_z,bincount_pt);
	AnalyzeBins(3,3,Outputdir,bincount_z,bincount_pt);
	AnalyzeBins(3,4,Outputdir,bincount_z,bincount_pt);
	AnalyzeBins(4,3,Outputdir,bincount_z,bincount_pt);
	AnalyzeBins(4,4,Outputdir,bincount_z,bincount_pt);
	AnalyzeBins(5,4,Outputdir,bincount_z,bincount_pt);
	//AnalyzeBinsZ(0,0,bincount_z,bincount_pt);
	
	}	
	



	private static ArrayList<GraphErrors> GetMultiplicityGenerated(H2F GenH2F, double electrons) {
		ArrayList<GraphErrors> Multiplicity = new ArrayList<GraphErrors>();

		for(int i=0; i<GenH2F.getSlicesX().size();i++) {
			double BinSize= GenH2F.getXAxis().getBinWidth(i)*GenH2F.getYAxis().getBinWidth(0);

			H1F SlicePionGen = GenH2F.getSlicesX().get(i);
			GraphErrors Graph= new GraphErrors();
			for (int k =0 ; k<SlicePionGen.getXaxis().getNBins(); k++) {
				
				double multiplicity = SlicePionGen.getBinContent(k)/(electrons*BinSize);
				double error = Math.pow(Math.sqrt(SlicePionGen.getBinContent(k))/(electrons*BinSize),2)+Math.pow(Math.sqrt(electrons)*SlicePionGen.getBinContent(k)/(electrons*electrons*BinSize), 2);
				Graph.addPoint(SlicePionGen.getXaxis().getBinCenter(k),multiplicity ,0,0);
			}
			Multiplicity.add(Graph);
		}
		return Multiplicity;
	}



	private static ArrayList<GraphErrors> GetMultiplicityGeneratedY(H2F GenH2F, double electrons) {
		ArrayList<GraphErrors> Multiplicity = new ArrayList<GraphErrors>();

		for(int i=0; i<GenH2F.getSlicesY().size();i++) {
			double BinSize= GenH2F.getXAxis().getBinWidth(i)*GenH2F.getYAxis().getBinWidth(0);

			H1F SlicePionGen = GenH2F.getSlicesY().get(i);
			GraphErrors Graph= new GraphErrors();
			for (int k =0 ; k<SlicePionGen.getYaxis().getNBins(); k++) {
				
				double multiplicity = SlicePionGen.getBinContent(k)/(electrons*BinSize);
				double error = Math.pow(Math.sqrt(SlicePionGen.getBinContent(k))/(electrons*BinSize),2)+Math.pow(Math.sqrt(electrons)*SlicePionGen.getBinContent(k)/(electrons*electrons*BinSize), 2);
				Graph.addPoint(SlicePionGen.getYaxis().getBinCenter(k),multiplicity ,0,0);
			}
			Multiplicity.add(Graph);
		}
		return Multiplicity;
	}





	private static ArrayList<GraphErrors> GetMultiplicity(String NameF,int BinID,String Outputdir, H2F H2Pion,H2F H2PionMC, H2F H2Gen,H2F H2GenREC, double electrons,H3F Q2Multi,H3F XMulti, H3F YMulti,H3F WMulti, H3F ZMulti, H3F PTMulti) throws FileNotFoundException {
		ArrayList<GraphErrors> Multiplicity = new ArrayList<GraphErrors>();
		//Hi
		     File folder = new File(Outputdir+"");
			folder.mkdir();
			File File_txt = new File(Outputdir+"/"+BinID+".txt");
			try {
				File_txt.createNewFile();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			if(File_txt.exists()) {
				System.out.println(" Output File Found" );
				PrintWriter outP = new PrintWriter(File_txt);
				NumberFormat nf = NumberFormat.getInstance();
		        nf.setMinimumFractionDigits(4);
				outP.println("======================================================================================================================================");
			outP.println("====================================     Multiplicity Analysis Bin nr  " + BinID+" from Giovanni Angelini code     ====================================");
			outP.println("===================================================================================================================================");
			outP.println("| Bin nr z | Bin nr PT^2 | <Q2> | <xB> | <W> | <y> | <z> | <PT^2> | purity | Multiplicity | δMulti stat |");

			outP.println("===================================================================================================================================");
		
			
			
		//	H2F BinPurity = new H2F("BinPurity",12 , 0.2,0.8,45 ,0,1.35);
		
			
		for(int i=0; i<H2Pion.getSlicesX().size();i++) {//Z
			double BinSize= H2Pion.getXAxis().getBinWidth(i)*H2Pion.getYAxis().getBinWidth(0);
			H1F SlicePionRec = H2Pion.getSlicesX().get(i);
			H1F SlicePionRecMC = H2PionMC.getSlicesX().get(i);
			H1F SlicePionGen = H2Gen.getSlicesX().get(i);
			H1F SlicePionGenandRec =H2GenREC.getSlicesX().get(i);
			
			GraphErrors Graph= new GraphErrors();
			
			
			
			for (int k =0 ; k<SlicePionRec.getXaxis().getNBins(); k++) {//PT
				double acceptance = SlicePionRecMC.getBinContent(k)/SlicePionGen.getBinContent(k);
				double purity=0;
				if(SlicePionRecMC.getBinContent(k)>20) {
				 purity = SlicePionGenandRec.getBinContent(k)/SlicePionRecMC.getBinContent(k);
				 if(purity<0.9 && purity>0.8) purity=purity+0.05;
				 else if (purity<0.8&&purity>=0.7)purity=purity+0.07;
				 else if (purity<0.7)purity=purity+0.07;
				}
		//		BinPurity.setBinContent(i, k, purity);
				if (acceptance>0.01 && acceptance<100) {
				double error2= (Math.pow(Math.sqrt(SlicePionRecMC.getBinContent(k))/SlicePionGen.getBinContent(k),2)+Math.pow((SlicePionRecMC.getBinContent(k)*Math.sqrt(SlicePionGen.getBinContent(k))/Math.pow(SlicePionRecMC.getBinContent(k),2)),2));
				double errorA=Math.sqrt(error2);
				double rescaledPions=	SlicePionRec.getBinContent(k)/acceptance;
				double multiplicity = rescaledPions/(electrons*(BinSize));
				double errorM1 = Math.sqrt(SlicePionRec.getBinContent(k))*(SlicePionGen.getBinContent(k)/SlicePionRecMC.getBinContent(k))*(1/(electrons*BinSize));
				double errorM2 = Math.sqrt(SlicePionGen.getBinContent(k))*(SlicePionRec.getBinContent(k)/SlicePionRecMC.getBinContent(k))*(1/(electrons*BinSize));
				double errorM3 = Math.sqrt(SlicePionRecMC.getBinContent(k))*((SlicePionRec.getBinContent(k)*SlicePionGen.getBinContent(k))/(SlicePionRecMC.getBinContent(k)*SlicePionRecMC.getBinContent(k)))*(1/(electrons*BinSize));
				double errorM4 = Math.sqrt(electrons)*((SlicePionRec.getBinContent(k)*SlicePionGen.getBinContent(k))/SlicePionRec.getBinContent(k))*(1/(electrons*electrons*BinSize));
				
				double errorMultiplicity = Math.sqrt(errorM4*errorM4+errorM3*errorM3+errorM2*errorM2+errorM1*errorM1);
				Graph.addPoint(SlicePionRec.getXaxis().getBinCenter(k),multiplicity ,0,errorMultiplicity);
				//if(BinID>12) System.out.println(" ========");

				//if(BinID>12) System.out.println(" Pion " + SlicePionRec.getBinContent(k) +" acceptance" + acceptance+ " generated"+SlicePionGen.getBinContent(k)+" rec mc "+  SlicePionRecMC.getBinContent(k));
				//if(BinID>12) System.out.println(" Electron  " + electrons +" bin size " + BinSize+ " denominator "+(electrons*(BinSize))+" rescaled pion"+  rescaledPions+ "Multi "+multiplicity );

				H1F ValoreQ2 = new H1F("ValoreQ2",200,1,12);
				H1F ValoreX = new H1F("ValoreX",200,0,1);
				H1F ValoreY = new H1F("ValoreY",200,0,1);
				H1F ValoreZ = new H1F("ValoreZ",200,0,1);
				H1F ValorePT = new H1F("ValorePT",200,0,1.5);
				H1F ValoreW = new H1F("ValoerW",200,1,12);
				/*
				System.out.println(" ====> Total ammount bins i "+ Q2Multi.getXAxis().getNBins() + " J "+ Q2Multi.getYAxis().getNBins()+ " K " +Q2Multi.getZAxis().getNBins() );
				for(int i=0 ; i < Q2Multi.getXAxis().getNBins(); i++) {
					for(int j=0 ; j < Q2Multi.getYAxis().getNBins(); j++) {
						for(int k=0 ; k < Q2Multi.getZAxis().getNBins(); k++) {
							if( Q2Multi.getBinContent(i, j, k)!=0) {
							System.out.println(" ==> i " + i+ " j " + j + " k "+ k + " value " + Q2Multi.getBinContent(i, j, k) );
							}
							}
						} 
				}
*/				
				double binsizey =1.0/199 ;
				double minvaluey=0.0;
				double binsizex =1.0/199 ;
				double minvaluex=0.0;
				double binsizeq2 =11.0/199 ;
				double minvalueq2=1.0;
				double binsizez =1.0/199 ;
				double minvaluez=0.0;
				double binsizept =1.5/199 ;
				double minvaluept=0.0;
				for (int qq=0; qq<200 ; qq ++) {
				
					double contoq2=minvalueq2+binsizeq2*qq;		
					double contox=minvaluex+binsizex*qq;
					double contoy=minvaluey+binsizey*qq;
					double contoz=minvaluez+binsizez*qq;
					double contopt=minvaluept+binsizept*qq;
					ValoreQ2.fill(contoq2,Q2Multi.getBinContent(i, k, qq));
					ValoreX.fill(contox,XMulti.getBinContent(i, k, qq));
					ValoreY.fill(contoy,YMulti.getBinContent(i, k, qq));
					ValoreZ.fill(contoz,ZMulti.getBinContent(i, k, qq));
					ValorePT.fill(contopt,PTMulti.getBinContent(i, k, qq));
					ValoreW.fill(contoq2,WMulti.getBinContent(i, k, qq));
				}
				
				//if(BinID>12) System.out.println(" Multi" + nf.format(multiplicity) + " multu non nf" +multiplicity + " error multi "  + nf.format(errorMultiplicity) );
				if (ValoreZ.getMean()!=0) {
					outP.println(i+" "+k+" " +nf.format(ValoreQ2.getMean())+" "+nf.format(ValoreX.getMean())+" "+nf.format(ValoreW.getMean())+" "+nf.format(ValoreY.getMean())+" "+nf.format(ValoreZ.getMean())+" "+nf.format(ValorePT.getMean())+" "+nf.format(purity) + " " + nf.format(multiplicity) + " " + nf.format(errorMultiplicity));
				}
				

				}
				
				}
			
			Multiplicity.add(Graph);
			  
		}
		String Title= "Purity "+BinID;
		//TGCanvas Purity = new TGCanvas(Title,Title,1600,1400);
		//BinPurity.setTitle(Title);
		//Purity.draw(BinPurity);
		//Purity.getPad().getAxisX().setRange(0, 0.6);
	//	Purity.getPad().getAxisX().setTitle(" z " );
	//	Purity.getPad().getAxisY().setTitle(" PT^2 " );
	//	Purity.getPad().setTitleFontSize(42);
	//	Purity.getPad().setAxisTitleFontSize(42);
	//	Purity.getPad().setAxisLabelFontSize(42);
		 outP.close();
		return Multiplicity;
		
				}
			else {
				System.out.println( " I couldnt create the file and therefore compute the Multiplicity");
				return null;
			}
	}


	private static ArrayList<GraphErrors> GetMultiplicityY(H2F H2Pion,H2F H2PionMC, H2F H2Gen, double electrons) {
		ArrayList<GraphErrors> Multiplicity = new ArrayList<GraphErrors>();
		for(int i=0; i<H2Pion.getSlicesY().size();i++) {
			double BinSize= H2Pion.getXAxis().getBinWidth(i)*H2Pion.getYAxis().getBinWidth(0);

			//System.out.println(" Slice H2Pion " + i+ " bin widht " + BinSize);
			H1F SlicePionRec = H2Pion.getSlicesY().get(i);
			H1F SlicePionRecMC = H2PionMC.getSlicesY().get(i);
			H1F SlicePionGen = H2Gen.getSlicesY().get(i);
			GraphErrors Graph= new GraphErrors();
			for (int k =0 ; k<SlicePionRec.getYaxis().getNBins(); k++) {
				double acceptance = SlicePionRecMC.getBinContent(k)/SlicePionGen.getBinContent(k);
				if (acceptance >0.01 && acceptance<100) {
				double error2= (Math.pow(Math.sqrt(SlicePionRecMC.getBinContent(k))/SlicePionGen.getBinContent(k),2)+Math.pow((SlicePionRecMC.getBinContent(k)*Math.sqrt(SlicePionGen.getBinContent(k))/Math.pow(SlicePionRecMC.getBinContent(k),2)),2));
				double errorA=Math.sqrt(error2);
				double rescaledPions=	SlicePionRec.getBinContent(k)/acceptance;
				double multiplicity = rescaledPions/(electrons*(BinSize));
				
				double errorMultp1 = Math.sqrt(SlicePionRec.getBinContent(k))/(acceptance*electrons*BinSize);
				double errorMultp2 = errorA*SlicePionRec.getBinContent(k)/(acceptance*acceptance*electrons*BinSize);
				double errorMultp3 = Math.sqrt(electrons)*SlicePionRec.getBinContent(k)/(acceptance*electrons*electrons*BinSize);
				double errorMultiplicity = Math.sqrt(errorMultp1*errorMultp1+errorMultp2*errorMultp2+errorMultp3*errorMultp3);
				Graph.addPoint(SlicePionRec.getXaxis().getBinCenter(k),multiplicity ,0,errorMultiplicity);
				}
				
				}
			
			Multiplicity.add(Graph);
		}
		return Multiplicity;
	}
	
	
	
	public static void AnalyzeBins(int binq, int binx , String Outputdir,int bincount_z, int bincount_pt ) throws FileNotFoundException {
		int BinNr = 0;
		if (binq==0 && binx==0)BinNr=1;
		else if (binq==1 && binx==0)BinNr=2;
		else if (binq==1 && binx==1)BinNr=3;
		else if (binq==1 && binx==2)BinNr=4;
		else if (binq==2 && binx==1)BinNr=5;
		else if (binq==2 && binx==2)BinNr=6;
		else if (binq==2 && binx==3)BinNr=7;
		else if (binq==2 && binx==4)BinNr=8;
		else if (binq==3 && binx==2)BinNr=9;
		else if (binq==3 && binx==3)BinNr=10;
		else if (binq==3 && binx==4)BinNr=11;
		else if (binq==4 && binx==3)BinNr=12;
		else if (binq==4 && binx==4)BinNr=13;
		else if (binq==5 && binx==4)BinNr=14;
		H2F GPions = PionsGenerated.GetH2F(binq,binx);
		H2F RPions = PionCounts.GetH2F(binq, binq);
		H2F DPions= PionCountsData.GetH2F(binq,binx);
		
		H3F Q2Multi = Q2Dist.GetH3F(binq, binx);
		H3F XMulti = XDist.GetH3F(binq, binx);
		H3F YMulti = yDist.GetH3F(binq, binx);
		H3F ZMulti = zDist.GetH3F(binq, binx);
		H3F PTMulti = PtDist.GetH3F(binq, binx);
		H3F WMulti= WDist.GetH3F(binq, binx);
		
		
		ArrayList<GraphErrors> GenPions = SliceHistogramX(PionsGenerated.GetH2F(binq, binx),bincount_z);
		ArrayList<GraphErrors> RecPions = SliceHistogramX(PionCounts.GetH2F(binq, binx),bincount_z);
		ArrayList<GraphErrors> DataPions = SliceHistogramX(PionCountsData.GetH2F(binq, binx),bincount_z);
		
		ArrayList<GraphErrors> MultiplicityGenerated = GetMultiplicityGenerated(PionsGenerated.GetH2F(binq, binx),Electron_counts.getBinContent(binq, binx));
		ArrayList<GraphErrors> MultiplicityData = GetMultiplicity("Data",BinNr,Outputdir,PionCountsData.GetH2F(binq, binx),PionCounts.GetH2F(binq, binx),PionsGenerated.GetH2F(binq, binx),PionsGeneratedReco.GetH2F(binq, binx),Electron_countsData.getBinContent(binq, binx),Q2Multi,XMulti,YMulti,WMulti,ZMulti,PTMulti);
		ArrayList<GraphErrors> Acceptance = GetAcceptance(PionCounts.GetH2F(binq, binx),PionsGenerated.GetH2F(binq, binx));

		
	for(int i=0; i <bincount_z ; i++) {	
		String Title = "Bin " +BinNr+ "Zbin "+i;
	//TGCanvas Z1 = new TGCanvas(Title,Title,2400,1400);Z1.divide(4, 2);

	
	
	RPions.setTitle("Reconstructed Pions from MC");	RPions.setTitleX(" z "); RPions.setTitleY(" PT^2 [GeV^2] "); 
	GPions.setTitle("Generated Pions in MC");GPions.setTitleX(" z "); GPions.setTitleY(" PT^2 [GeV^2] "); 
	DPions.setTitle("Data Rec Pions from MC");DPions.setTitleX(" z "); DPions.setTitleY(" PT^2 [GeV^2] "); 
	
	
	//Z1.cd(0);Z1.draw(GPions);Z1.getPad().getAxisZ().setLog(true);
	//Z1.cd(1);Z1.draw(RPions);Z1.getPad().getAxisZ().setLog(true);
	//Z1.cd(2);Z1.draw(DPions);Z1.getPad().getAxisZ().setLog(true);


	
	
	GraphErrors SlicedGenerated= GenPions.get(i);
	GraphErrors SlicedReconstruction= RecPions.get(i);
	GraphErrors SlicedData= DataPions.get(i);
	
	SlicedGenerated.setTitle(" Generated pion slice z "+i );SlicedGenerated.setTitleX("PT^2 [GeV^2]");
	SlicedReconstruction.setTitle(" Reconstruction pion slice z "+i );SlicedReconstruction.setTitleX("PT^2 [GeV^2]");
	SlicedData.setTitle(" Data Rec pion slice z "+i );SlicedData.setTitleX("PT^2 [GeV^2]");
	SlicedData.setMarkerColor(2);

	
	//Z1.cd(3);Z1.draw(SlicedGenerated);
	//Z1.cd(4);Z1.draw(SlicedReconstruction);Z1.draw(SlicedData,"same");Z1.getPad().getAxisY().setLog(true);
	
	GraphErrors SlicedAcceptance= Acceptance.get(i);
	SlicedAcceptance.setTitle(" Acceptance pion slice z "+i );
	SlicedAcceptance.setTitleX("PT^2 [GeV^2]");SlicedAcceptance.setTitleY("Acceptance+Efficiency");
	//Z1.cd(5);Z1.draw(SlicedAcceptance);
	//Z1.getPad().getAxisY().setRange(0.0,1);
	
	

	
	
	GraphErrors SlicedMultiplicityGenerated =  MultiplicityGenerated.get(i);
	GraphErrors SlicedMultiplicityData = MultiplicityData.get(i);
	SlicedMultiplicityGenerated.setMarkerColor(2);SlicedMultiplicityGenerated.setMarkerSize(4);SlicedMultiplicityGenerated.setLineThickness(3);
	//Z1.cd(6); Z1.draw(SlicedMultiplicityData);  Z1.draw(SlicedMultiplicityGenerated,"same,L");
	//Z1.getPad().getAxisY().setRange(0.01,30);
	//Z1.getPad().getAxisY().setLog(true);
	//Z1.getPad().getAxisY().setTitle(" Multiplicity ");
	//Z1.getPad().getAxisX().setTitle(" PT^2 [GeV^2]");
	
	
		}
	}
	
	public static void AnalyzeBinsZ(int binq, int binx , int bincount_z, int bincount_pt ) {
		int BinNr = 0;
		if (binq==0 && binx==0)BinNr=1;
		else if (binq==1 && binx==0)BinNr=2;
		else if (binq==1 && binx==1)BinNr=3;
		else if (binq==1 && binx==2)BinNr=4;
		else if (binq==2 && binx==1)BinNr=5;
		else if (binq==2 && binx==2)BinNr=6;
		else if (binq==2 && binx==3)BinNr=7;
		else if (binq==2 && binx==4)BinNr=8;
		else if (binq==3 && binx==2)BinNr=9;
		else if (binq==3 && binx==3)BinNr=10;
		else if (binq==3 && binx==4)BinNr=11;
		else if (binq==4 && binx==3)BinNr=12;
		else if (binq==4 && binx==4)BinNr=13;
		else if (binq==5 && binx==4)BinNr=14;
				
		//TGCanvas ZDep = new TGCanvas("ZDep","ZDep",2400,1400);
		//ZDep.divide(6, 6);

		ArrayList<GraphErrors> GenPions = SliceHistogramY(PionsGenerated.GetH2F(binq, binx));
		ArrayList<GraphErrors> RecPions = SliceHistogramY(PionCounts.GetH2F(binq, binx));
		ArrayList<GraphErrors> DataPions = SliceHistogramY(PionCountsData.GetH2F(binx, binx));
		ArrayList<GraphErrors> Acceptance = GetAcceptanceY(PionCounts.GetH2F(binq, binx),PionsGenerated.GetH2F(binq, binx));
		ArrayList<GraphErrors> MultiplicityGenerated = GetMultiplicityGeneratedY(PionsGenerated.GetH2F(binq, binx),Electron_counts.getBinContent(binq, binx));
		ArrayList<GraphErrors> MultiplicityData = GetMultiplicityY(PionCountsData.GetH2F(binq, binx),PionCounts.GetH2F(binq, binx),PionsGenerated.GetH2F(binq, binx),Electron_countsData.getBinContent(binq, binx));

		for(int i=0; i <bincount_pt ; i++) {	
			GraphErrors SlicedGenerated= GenPions.get(i);
			GraphErrors SlicedReconstruction= RecPions.get(i);
			GraphErrors SlicedData= DataPions.get(i);
			
			GraphErrors SlicedAcceptance= Acceptance.get(i);
			GraphErrors SlicedMultiGen = MultiplicityGenerated.get(i);
			GraphErrors SlicedMultiData= MultiplicityData.get(i);
			
			//ZDep.cd(i);ZDep.draw(SlicedMultiData);ZDep.draw(SlicedMultiGen,"same,L");
		}


	}
	
	
	
	private static ArrayList<GraphErrors> GetAcceptance(H2F H2Pion,
			H2F H2Gen) {
		ArrayList<GraphErrors> Accettanze = new ArrayList<GraphErrors>();

		for(int i=0; i<H2Pion.getSlicesX().size();i++) {
			H1F SlicePionRec = H2Pion.getSlicesX().get(i);
			H1F SlicePionGen = H2Gen.getSlicesX().get(i);
			GraphErrors Graph= new GraphErrors();
			for (int k =0 ; k<SlicePionRec.getXaxis().getNBins(); k++) {
				double ratio = SlicePionRec.getBinContent(k)/SlicePionGen.getBinContent(k);
				double error2= (Math.pow(Math.sqrt(SlicePionRec.getBinContent(k))/SlicePionGen.getBinContent(k),2)+Math.pow((SlicePionRec.getBinContent(k)*Math.sqrt(SlicePionGen.getBinContent(k))/Math.pow(SlicePionGen.getBinContent(k),2)),2));
				Graph.addPoint(SlicePionRec.getXaxis().getBinCenter(k),ratio ,0, Math.sqrt(error2));
			}
			Accettanze.add(Graph);
		}
	
		return Accettanze;
	}

	private static ArrayList<GraphErrors> GetAcceptanceY(H2F H2Pion,
			H2F H2Gen) {
		ArrayList<GraphErrors> Accettanze = new ArrayList<GraphErrors>();

		for(int i=0; i<H2Pion.getSlicesY().size();i++) {
			H1F SlicePionRec = H2Pion.getSlicesY().get(i);
			H1F SlicePionGen = H2Gen.getSlicesY().get(i);
			GraphErrors Graph= new GraphErrors();
			for (int k =0 ; k<SlicePionRec.getYaxis().getNBins(); k++) {
				double ratio = SlicePionRec.getBinContent(k)/SlicePionGen.getBinContent(k);
				double error2= (Math.pow(Math.sqrt(SlicePionRec.getBinContent(k))/SlicePionGen.getBinContent(k),2)+Math.pow((SlicePionRec.getBinContent(k)*Math.sqrt(SlicePionGen.getBinContent(k))/Math.pow(SlicePionGen.getBinContent(k),2)),2));
				Graph.addPoint(SlicePionRec.getYaxis().getBinCenter(k),ratio ,0, Math.sqrt(error2));
			}
			Accettanze.add(Graph);
		}
	
		return Accettanze;
	}
	




	private static ArrayList<GraphErrors> SliceHistogramX(H2F getH2F, int bincount_z) {
	ArrayList<GraphErrors> MyGraphs = new ArrayList<GraphErrors>();
	for(int i=0; i<getH2F.getSlicesX().size();i++) {
		H1F Slice = getH2F.getSlicesX().get(i);
		GraphErrors Graph= new GraphErrors();
		for (int k =0 ; k<Slice.getXaxis().getNBins(); k++) {
			Graph.addPoint(Slice.getXaxis().getBinCenter(k),Slice.getBinContent(k) ,0, Math.sqrt(Slice.getBinContent(k)));
		}
		MyGraphs.add(Graph);
	}
	return MyGraphs;
	}

	private static ArrayList<GraphErrors> SliceHistogramY(H2F getH2F) {
		ArrayList<GraphErrors> MyGraphs = new ArrayList<GraphErrors>();
		for(int i=0; i<getH2F.getSlicesY().size();i++) {
			H1F Slice = getH2F.getSlicesY().get(i);
			GraphErrors Graph= new GraphErrors();
			for (int k =0 ; k<Slice.getYaxis().getNBins(); k++) {
				Graph.addPoint(Slice.getYaxis().getBinCenter(k),Slice.getBinContent(k) ,0, Math.sqrt(Slice.getBinContent(k)));
			}
			MyGraphs.add(Graph);
		}
		return MyGraphs;
		}



	void InitializeHistos() {
		
	}
	
	private static void InitializeMultiBins(boolean PolygonalBinClas, int bincount_q2,int bincount_x, int bincount_z, int bincount_pt, int binmass, int bincount_phi,
			double[] bin_q2,double[] bin_x, double[] bin_z, double[] bin_pt, double[] bin_mass, double[] bin_phi) {
		
		PionCounts = new MultiBins(6,7,bincount_z,bincount_pt,0);
		PionCounts.SetName_1stVariable("Q2");PionCounts.SetName_2ndVariable("Xb");
		PionCounts.SetName_3rdVariable("z");PionCounts.SetName_4thVariable("Pt");	 System.out.println( " Multidimensional Analysis code ");
		System.out.println(" Bin1: " + PionCounts.GetName_1stVariable() + " Bin2: " + PionCounts.GetName_2ndVariable());
		System.out.println(" Bin3: " + PionCounts.GetName_3rdVariable() + " Bin4: " + PionCounts.GetName_4thVariable());
		PionCounts.SetUnevenBins(bin_q2,bin_x,bin_z,bin_pt,bin_mass);
		PionCounts.GenerateHistograms("PionCounts");
		if(PolygonalBinClas==true) { PionCounts.InizializeClasPolygons();}
		
		
		PionCountsData = new MultiBins(6,7,bincount_z,bincount_pt,0);
		PionCountsData.SetName_1stVariable("Q2");PionCountsData.SetName_2ndVariable("Xb");
		PionCountsData.SetName_3rdVariable("z");PionCountsData.SetName_4thVariable("Pt");	 System.out.println( " Multidimensional Analysis code ");
		System.out.println(" Bin1: " + PionCountsData.GetName_1stVariable() + " Bin2: " + PionCountsData.GetName_2ndVariable());
		System.out.println(" Bin3: " + PionCountsData.GetName_3rdVariable() + " Bin4: " + PionCountsData.GetName_4thVariable());
		PionCountsData.SetUnevenBins(bin_q2,bin_x,bin_z,bin_pt,bin_mass);
		PionCountsData.GenerateHistograms("PionCountsData");
		if(PolygonalBinClas==true) { PionCountsData.InizializeClasPolygons();}
		
		PionsGenerated = new MultiBins(6,7,bincount_z,bincount_pt,0);
		PionsGenerated.SetName_1stVariable("Q2");PionsGenerated.SetName_2ndVariable("Xb");
		PionsGenerated.SetName_3rdVariable("z");PionsGenerated.SetName_4thVariable("Pt");	 
		PionsGenerated.SetUnevenBins(bin_q2,bin_x,bin_z,bin_pt,bin_mass);
		PionsGenerated.GenerateHistograms("PionsGenerated");
		if(PolygonalBinClas==true) { PionsGenerated.InizializeClasPolygons();}
				
		
		PionsGeneratedReco = new MultiBins(6,7,bincount_z,bincount_pt,0);;	 
		PionsGeneratedReco.SetUnevenBins(bin_q2,bin_x,bin_z,bin_pt,bin_mass);
		PionsGeneratedReco.GenerateHistograms("PionsGeneratedReco");
		if(PolygonalBinClas==true) { PionsGeneratedReco.InizializeClasPolygons();}
		/*
		Electron_counts = new MultiBins(6,7,0,0,0);
		Electron_counts.SetName_1stVariable("Q2");Electron_counts.SetName_2ndVariable("Xb");
		Electron_counts.SetUnevenBins(bin_q2,bin_x, new double[]{0,0}, new double[]{0,0}, new double[]{0,0} );
		Electron_counts.GenerateHistograms("Electron_counts");
		System.out.println(Electron_counts.FstBin+" max "+ Electron_counts.GetUnevenBins(1)[0]+" "+ Electron_counts.GetUnevenBins(1)[1]);
		if(PolygonalBinClas==true) { Electron_counts.SetPolygonal(PolygonalBinClas); }
		
		Electron_countsData = new MultiBins(6,7,0,0,0);
		Electron_countsData.SetName_1stVariable("Q2");Electron_countsData.SetName_2ndVariable("Xb");
		Electron_countsData.SetUnevenBins(bin_q2,bin_x, new double[]{0,0}, new double[]{0,0}, new double[]{0,0} );
		Electron_countsData.GenerateHistograms("Electron_counts");
		System.out.println(Electron_countsData.FstBin+" max "+ Electron_countsData.GetUnevenBins(1)[0]+" "+ Electron_countsData.GetUnevenBins(1)[1]);
		if(PolygonalBinClas==true) { Electron_countsData.SetPolygonal(PolygonalBinClas); }
		
		*/
		double []bin200Q2 =new double[200];
		for(int i=0; i<200; i++) {
			double binsize =11.0/199 ;
			double minvalue=1.0;
			double conto=minvalue+binsize*i;
			bin200Q2[i]=conto;
		}
		Q2Dist = new MultiBins(6,7,bincount_z,bincount_pt,200);
		Q2Dist.SetUnevenBins(bin_q2,bin_x,bin_z,bin_pt,bin200Q2);
		Q2Dist.GenerateHistograms("Q2Dist");
		if(PolygonalBinClas==true) { Q2Dist.InizializeClasPolygons();}
		
		WDist = new MultiBins(6,7,bincount_z,bincount_pt,200);
		WDist.SetUnevenBins(bin_q2,bin_x,bin_z,bin_pt,bin200Q2);
		WDist.GenerateHistograms("WDist");
		if(PolygonalBinClas==true) { WDist.InizializeClasPolygons();}
		
		
		
		double []bin200X =new double[200];
		for(int i=0; i<200; i++) {
			double binsize =1.0/199 ;
			double minvalue=0.0;
			double conto=minvalue+binsize*i;
			bin200X[i]=conto;
		}
		XDist = new MultiBins(6,7,bincount_z,bincount_pt,200);
		XDist.SetUnevenBins(bin_q2,bin_x,bin_z,bin_pt,bin200X);
		XDist.GenerateHistograms("XDist");
		if(PolygonalBinClas==true) { XDist.InizializeClasPolygons();}
		
		double []bin200y =new double[200];
		for(int i=0; i<200; i++) {
			double binsize =1.0/199 ;
			double minvalue=0.0;
			double conto=minvalue+binsize*i;
			bin200y[i]=conto;
		}
		yDist = new MultiBins(6,7,bincount_z,bincount_pt,200);
		yDist.SetUnevenBins(bin_q2,bin_x,bin_z,bin_pt,bin200y);
		yDist.GenerateHistograms("yDist");
		if(PolygonalBinClas==true) { yDist.InizializeClasPolygons();}
		
		
		double []bin200Z= new double[200];
		for(int i=0; i<200; i++) {
			double binsize =1.0/199 ;
			double minvalue=0.0;
			double conto=minvalue+binsize*i;
			bin200Z[i]=conto;
		}
		zDist = new MultiBins(6,7,bincount_z,bincount_pt,200);
		zDist.SetUnevenBins(bin_q2,bin_x,bin_z,bin_pt,bin200Z);
		zDist.GenerateHistograms("zDist");
		if(PolygonalBinClas==true) { zDist.InizializeClasPolygons();}
	
		double []bin200PT= new double[200];
		for(int i=0; i<200; i++) {
			double binsize =1.5/199 ;
			double minvalue=0.0;
			double conto=minvalue+binsize*i;
			bin200PT[i]=conto;
		}
		PtDist = new MultiBins(6,7,bincount_z,bincount_pt,200);
		PtDist.SetUnevenBins(bin_q2,bin_x,bin_z,bin_pt,bin200PT);
		PtDist.GenerateHistograms("PtDist");
		if(PolygonalBinClas==true) { PtDist.InizializeClasPolygons();}
	
		
		double []bin200Eta=new double[200];
		for(int i=0; i<200; i++) {
			double binsize =10/199 ;
			double minvalue=-10.0;
			double conto=minvalue+binsize*i;
			bin200Eta[i]=conto;
		}
		EtaDist = new MultiBins(6,7,bincount_z,bincount_pt,200);
	
		EtaDist.SetUnevenBins(bin_q2,bin_x,bin_z,bin_pt,bin200Eta);
		EtaDist.GenerateHistograms("EtaDist");
		if(PolygonalBinClas==true) { EtaDist.InizializeClasPolygons();}

		
	}

	
}
