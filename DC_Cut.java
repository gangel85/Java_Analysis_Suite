package PHD_Clean;

import java.util.ArrayList;

import org.jlab.groot.data.H2F;
import org.jlab.groot.math.Axis;
import org.jlab.groot.math.F1D;
import org.jlab.jnp.physics.Vector3;

import eu.mihosoft.vrl.v3d.Polygon;
import eu.mihosoft.vrl.v3d.Vector3d;

public class DC_Cut implements DetectorCut {
	//Paticle excluded by the cuts
	H2F DC1Histo_Pre = new H2F("DC1Pre",200,-200,200,200,-200,200);
	H2F DC2Histo_Pre = new H2F("DC2Pre",200,-300,300,200,-300,300);
	H2F DC3Histo_Pre = new H2F("DC3Pre",200,-500,500,200,-500,500);
	//Particle that surves the cuts 
	H2F DC1Histo_Aft = new H2F("DC1Aft",200,-200,200,200,-200,200);
	H2F DC2Histo_Aft = new H2F("DC2Aft",200,-300,300,200,-300,300);
	H2F DC3Histo_Aft = new H2F("DC3Aft",200,-500,500,200,-500,500);
    ArrayList<H2F> Histograms_Cut_Pre = new ArrayList<H2F>();
    ArrayList<H2F> Histograms_Cut_Aft = new ArrayList<H2F>();
	private String name;
	
	// THESE ARE INTERNAL VALUES TO MACH THE PARTICLE REC DEFINIOTNI
	private int DC1=61; //The index used to identify the first DC Region
	private int DC2=62 ;//The index useto identify the second DC region 
	private int DC3=63; //The index used to identify the third DC region 	 
	// A list of boolean used in the code
	private boolean DC1_cut=false;
	private boolean DC2_cut=false;
	private boolean DC3_cut=false;
	// Variables used in the cuts 
	private int elDCRoutine;
	private int elDCCutValue;
	private int hadronDCRoutine;
	private int hadronDCCutValue;
	private boolean MC;
	
	//These numbers are form Stefan's Diehl studies and plot on RG-A common analysis 
	double[][] e_min_rad= {{50,30,30},{50,30,30}};
    double[][] e_max_rad={{115,300,300},{115,300,300}};
    
    
	double[][] pim_min_rad= {{78,110,30},{78,110,30}};
    double[][] pim_max_rad={{140,300,300},{140,300,300}};
     
	double[][][] e_x_y_min_const = {{{19,15.3,24.2},{14.2,11.2,23.15},{16.9,13.9,31.5},{15,11.5,24.00},{17.5,20.4,24.11},{17.0,12.5,23.3}},{{16,19,24},{16,19,24},{16,19,24},{16,19,24},{16,19,24},{16,19,24}}};
	double[][][] e_x_y_min_slope = {{{-0.679,-0.563,-0.61},{-0.640,-0.552,-0.60},{-0.642,-0.559,-0.65},{-0.640,-0.570,-0.61},{-0.662,-0.610,-0.611},{-0.680,-0.645,-0.636}},{{-0.64,-0.61,-0.605},{-0.64,-0.61,-0.605},{-0.64,-0.61,-0.605},{-0.64,-0.61,-0.605},{-0.64,-0.61,-0.605},{-0.64,-0.61,-0.605}}};

	double[][][] e_x_y_max_const = {{{-20,-18.5,-30.1},{-15.1,-16.9,-29.7},{-22.6,-23.3,-29.5},{-16.2,-19.5,-29.3},{-21.1,-28.1,-30.5},{-19.3,-20.7,-23.2}},{{-16,-19,-24},{-16,-19,-24},{-16,-19,-24},{-16,-19,-24},{-16,-19,-24},{-16,-19,-24}}};
	double[][][] e_x_y_max_slope = {{{0.677,0.563,0.61},{0.657,0.594,0.65},{0.715,0.622,0.649},{0.660,0.60,0.63},{0.691,0.687,0.598},{0.670,0.605,0.657}},{{0.64,0.61,0.605},{0.64,0.61,0.605},{0.64,0.61,0.605},{0.64,0.61,0.605},{0.64,0.61,0.605},{0.64,0.61,0.605}}};

		
	double[][] pip_min_rad= {{60,90,100},{60,90,100}};
    double[][] pip_max_rad={{118,178,280},{118,178,280}};


    
	double[][][] pip_x_y_min_const = {{{4,8,8},{4,8,8},{4,8,8},{4,8,8},{4,8,8},{4,8,8}},{{4,8,8},{4,8,8},{4,8,8},{4,8,8},{4,8,8},{4,8,8}}};
	double[][][] pip_x_y_min_slope = {{{-0.55,-0.55,-0.55},{-0.55,-0.55,-0.55},{-0.55,-0.55,-0.55},{-0.55,-0.55,-0.55},{-0.55,-0.55,-0.55},{-0.55,-0.55,-0.55}},{{-0.55,-0.55,-0.55},{-0.55,-0.55,-0.55},{-0.55,-0.55,-0.55},{-0.55,-0.55,-0.55},{-0.55,-0.55,-0.55},{-0.55,-0.55,-0.55}}};

	double[][][] pip_x_y_max_const = {{{-5,-9,-9},{-5,-9,-9},{-5,-9,-9},{-5,-9,-9},{-5,-9,-9},{-5,-9,-9}},{{-5,-9,-9},{-5,-9,-9},{-5,-9,-9},{-5,-9,-9},{-5,-9,-9},{-5,-9,-9}}};
	double[][][] pip_x_y_max_slope = {{{0.55,0.55,0.55},{0.55,0.55,0.55},{0.55,0.55,0.55},{0.55,0.55,0.55},{0.55,0.55,0.55},{0.55,0.55,0.55}},{{0.55,0.55,0.55},{0.55,0.55,0.55},{0.55,0.55,0.55},{0.55,0.55,0.55},{0.55,0.55,0.55},{0.55,0.55,0.55}}};
	
	
		 			 // PID needed for fiducial cuits
		 			 int pid =0; 
		 			 int sector =0 ; 
	
	/* ---- Constructors
       ----         	 */

	DC_Cut(ReadParameters parameters){
		elDCRoutine = parameters.Get_el_DCRoutine();
		elDCCutValue = parameters.Get_el_DCCUT();
		hadronDCRoutine = parameters.Get_Hadron_DCRoutine();
		hadronDCCutValue= parameters.Get_Hadron_DCCUT();
		MC = parameters.GetMCBoolean();
		//This constructor initizializes standard parameter for the cuts
	
	}

	public void setName(String nameDetector) {
		this.name=nameDetector;
	}
	public String toString() {
		return this.name;
	}
	
	
/**
 * This function return the sector of a given hit .
 * @param Hit is a Vector3 representing x,y,z of a detector reconstructed cluster
 **/

	public int getSector(Vector3 Hit) {
		double myphi = Math.toDegrees(Hit.phi());
		int sec =0;
		if (myphi>0 && myphi<=30) sec= 1;
		else if (myphi>30 && myphi <=90) sec= 2;
		else if (myphi>90 && myphi <=150 ) sec= 3;
		else if (myphi>150 && myphi<= 180 ) sec=4;
		else if ( myphi>-30 && myphi<=0) sec=1;
		else if (myphi>-90 && myphi<=-30) sec=6;
		else if (myphi>-150 && myphi <= -90) sec=5;
		else if ( myphi>=-180 && myphi <=-150) sec= 4;
		//if(sec==0) System.out.println( " AO 0 " + Hit.x()+" "+Hit.y()+" "+Hit.z());
		return sec;  
		}

	/**
	 * This function check if the particle passes the cuts
	 **/
	public boolean Status(ParticleREC particle) {
		if(particle.pid()==11 && elDCCutValue!=0)return DCCut(particle);
		else if (particle.pid()!=11 && elDCCutValue!=0)return DCCut(particle);
		else return true; 
	}

	//--- Internal functions
	private boolean DCCut(ParticleREC particle) {
		
		if ( particle.pid()==11) {
			this.pid=0; 
		//	System.out.println(" I am looking at electrons");
		}
		else if ( particle.pid() == 2212) {
			this.pid=1; 
			System.out.println(" I am looking at 2212");
		}
		else if ( particle.pid() == 211) {
			this.pid=1 	; 
			//System.out.println(" I am looking at pi+ 211");
		}
		else if ( particle.pid() == -211) {
			this.pid=2; 
				//We set pi-  == to electron so to define the fiducial cuts 
 			//System.out.println(" I am looking at pi+ 211");
		}
		else if ( particle.pid() == 321) {
			this.pid=4; System.out.println(" I am looking at K id 321");
		}
		else if ( particle.pid() == -321) {
			this.pid=5; System.out.println(" I am looking at K id -321");
		}
		//Vector3 Hit = new Vector3();
			// cuts only on region 2 
			int settore=0;
			if(particle.hasDetector(DC2)==true){
		   Vector3	Hit2=particle.getDetectorHits(DC2);
			 settore = this.getSector(Hit2);
			//System.out.println(" Ho un settore"); 
			}
		 else {
			 return false;}

			if(particle.hasDetector(DC1)==true  && settore !=0){
				Vector3 Hit=particle.getDetectorHits(DC1);
			
				if(this.pid==0) {		
					//System.out.println(" electron - DC ROUTINE "+ elDCRoutine);
					//Use XYZ Cuts
					 DC1_cut=DC_fiducial_cut_XY(Hit,1,this.pid,settore);
					//Use Polar cuts 
						
				}
				else {	
					 DC1_cut=DC_fiducial_cut_XY(Hit,1,this.pid,settore);
						
					}
				}
			if(particle.hasDetector(DC2)==true && settore !=0){
				Vector3 Hit=particle.getDetectorHits(DC2);
				if(this.pid==0) {		
					//Use XYZ Cuts
					 DC2_cut=DC_fiducial_cut_XY(Hit,2,this.pid,settore);
					//Use Polar cuts 
						
				}
				else {	
			DC2_cut=DC_fiducial_cut_XY(Hit,2,this.pid,settore);
						
					}
				
			}
			if(particle.hasDetector(DC3)==true && settore !=0){				
				Vector3 Hit=particle.getDetectorHits(DC3);
				if(this.pid==0) {		
					//Use XYZ Cuts
					 DC3_cut=DC_fiducial_cut_XY(Hit,3,this.pid,settore);
				
					
				}
				else {	
					 DC3_cut=DC_fiducial_cut_XY(Hit,3,this.pid,settore);
					
					}
			}
			if(DC1_cut==true && DC2_cut==true && DC3_cut==true)  { 
			return true; } // all thre cuts passed  }
			return false;	
	}

	//pid is converted pid from the Stefan definition 

	boolean DC_fiducial_cut_XY(Vector3 Hit,int region, int pid, int settore){
				double X= Hit.x();
				double Y= Hit.y();

				if(settore == 2)
				{
				 double X_new = X * Math.cos(-60.0 * Math.PI / 180.0) - Y * Math.sin(-60 * Math.PI / 180);
				Y = X * Math.sin(-60.0 * Math.PI / 180.0) + Y * Math.cos(-60.0 * Math.PI / 180.0);
				X = X_new;
				}

				if(settore == 3)
				{
				 double X_new = X * Math.cos(-120.0 * Math.PI / 180) - Y * Math.sin(-120 * Math.PI / 180);
				Y = X * Math.sin(-120.0 * Math.PI / 180.0) + Y * Math.cos(-120 * Math.PI / 180.0);
				X = X_new;
				}

				if(settore == 4)
				{
					
				X = -X ;
				Y= - Y ; 
				
				}

				if(settore == 5)
				{
				 double X_new = X * Math.cos(120 * Math.PI / 180) - Y * Math.sin(120 * Math.PI / 180);
				Y = X * Math.sin(120 * Math.PI / 180) + Y * Math.cos(120 * Math.PI / 180);
				X = X_new;
				}

				if(settore == 6)
				{
				 double X_new = X * Math.cos(60.0 * Math.PI / 180.0) - Y * Math.sin(60.0 * Math.PI / 180.0);
				Y = X * Math.sin(60 * Math.PI / 180.0) + Y * Math.cos(60 * Math.PI / 180);
				X = X_new;
				//System.out.println(" I am in sector 6 " + X + " y " + Y);
				}

				--region;
				double Constant_Min=0, Constant_Max=0, Slope_Min=0, Slope_Max=0,Radius_Min=0,Radius_Max=0;
				// electron
				if(pid==0) {
					if (MC==true) {
						if(elDCCutValue==2) {
					Constant_Min=e_x_y_min_const[1][settore-1][region];
					Constant_Max=e_x_y_max_const[1][settore-1][region];
					Slope_Min=e_x_y_min_slope[1][settore-1][region];
					Slope_Max=e_x_y_max_slope[1][settore-1][region];
					Radius_Min=e_min_rad[1][region];
					Radius_Max=e_max_rad[1][region];
						}
						else if (elDCCutValue==1) {
							Constant_Min=e_x_y_min_const[1][settore-1][region];
							Constant_Max=e_x_y_max_const[1][settore-1][region];
							Constant_Min=Constant_Min-2;
							Constant_Max=Constant_Max+2;
							Slope_Min=e_x_y_min_slope[1][settore-1][region];
							Slope_Max=e_x_y_max_slope[1][settore-1][region];
							Radius_Min=e_min_rad[1][region];
							Radius_Max=e_max_rad[1][region];
						}
						else if (elDCCutValue==3) {
							Constant_Min=e_x_y_min_const[1][settore-1][region];
							Constant_Max=e_x_y_max_const[1][settore-1][region];
							Constant_Min=Constant_Min+2;
							Constant_Max=Constant_Max-2;
							Slope_Min=e_x_y_min_slope[1][settore-1][region];
							Slope_Max=e_x_y_max_slope[1][settore-1][region];
							Radius_Min=e_min_rad[1][region];
							Radius_Max=e_max_rad[1][region];
						}
					}
					else {
						Constant_Min=e_x_y_min_const[0][settore-1][region];
						Constant_Max=e_x_y_max_const[0][settore-1][region];
						Slope_Min=e_x_y_min_slope[0][settore-1][region];
						Slope_Max=e_x_y_max_slope[0][settore-1][region];
						Radius_Min=e_min_rad[0][region];
						Radius_Max=e_max_rad[0][region];
					}
				}
						//pi p 
				else if (pid==1) {
					if (MC==true) {
						if(hadronDCCutValue==2) {
						Constant_Min=pip_x_y_min_const[1][settore-1][region];
						Constant_Max=pip_x_y_max_const[1][settore-1][region];
						Slope_Min=pip_x_y_min_slope[1][settore-1][region];
						Slope_Max=pip_x_y_max_slope[1][settore-1][region];
						Radius_Min=pip_min_rad[1][region];
						Radius_Max=pip_max_rad[1][region];
						}
						else if (hadronDCCutValue==1) {
							Constant_Min=pip_x_y_min_const[1][settore-1][region];
							Constant_Max=pip_x_y_max_const[1][settore-1][region];
							Constant_Min=Constant_Min-2;
							Constant_Max=Constant_Max+2;
							Slope_Min=pip_x_y_min_slope[1][settore-1][region];
							Slope_Max=pip_x_y_max_slope[1][settore-1][region];
							Radius_Min=pip_min_rad[1][region];
							Radius_Max=pip_max_rad[1][region];
						}
						else if (hadronDCCutValue==3) {
							Constant_Min=pip_x_y_min_const[1][settore-1][region];
							Constant_Max=pip_x_y_max_const[1][settore-1][region];
							Constant_Min=Constant_Min+2;
							Constant_Max=Constant_Max-2;
							Slope_Min=pip_x_y_min_slope[1][settore-1][region];
							Slope_Max=pip_x_y_max_slope[1][settore-1][region];
							Radius_Min=pip_min_rad[1][region];
							Radius_Max=pip_max_rad[1][region];
						}
						}
						else {
							Constant_Min=pip_x_y_min_const[0][settore-1][region];
							Constant_Max=pip_x_y_max_const[0][settore-1][region];
							Slope_Min=pip_x_y_min_slope[0][settore-1][region];
							Slope_Max=pip_x_y_max_slope[0][settore-1][region];
							Radius_Min=pip_min_rad[0][region];
							Radius_Max=pip_max_rad[0][region];
						}
				}
				
				
				
				else if (pid==2) {
					if (MC==true) {
						Constant_Min=e_x_y_min_const[1][settore-1][region];
						Constant_Max=e_x_y_max_const[1][settore-1][region];
						Slope_Min=e_x_y_min_slope[1][settore-1][region];
						Slope_Max=e_x_y_max_slope[1][settore-1][region];
						Radius_Min=pim_min_rad[1][region];
						Radius_Max=pim_max_rad[1][region];
						}
						else {
							Constant_Min=e_x_y_min_const[0][settore-1][region];
							Constant_Max=e_x_y_max_const[0][settore-1][region];
							Slope_Min=e_x_y_min_slope[0][settore-1][region];
							Slope_Max=e_x_y_max_slope[0][settore-1][region];
							Radius_Min=pim_min_rad[0][region];
							Radius_Max=pim_max_rad[0][region];
						}
				}
				
				
				double calc_min = Constant_Min + Slope_Min * X;
				double calc_max = Constant_Max + Slope_Max * X;
				//System.out.println(" Cal Min " + calc_min + " Cal Max" + calc_max + " Equation 1 " +Constant_Min +" + "+Slope_Min+"*x"+" Equation 2 " +Constant_Max +" + "+Slope_Max+"*x"  );
				
			//	double calc_min = x_y_minparams_in[pid][settore - 1][region][0] + x_y_minparams_in[pid][settore - 1][region][1] * X;
			//	double calc_max = x_y_maxparams_in[pid][settore - 1][region][0] + x_y_maxparams_in[pid][settore - 1][region][1] * X;
				
				if ((Y > calc_min) && (Y < calc_max) && Math.sqrt(Hit.x()*Hit.x()+Hit.y()*Hit.y())>Radius_Min &&  Math.sqrt(Hit.x()*Hit.x()+Hit.y()*Hit.y())<Radius_Max ) 
					{
					//if (region ==0 && pid==1)System.out.println(" I am here babe! Angle   " + Math.toDegrees(Hit.theta())+ " radius " +Math.sqrt(Hit.x()*Hit.x()+Hit.y()*Hit.y())  + " PID "+ pid); 
					if (region ==0 )DC1Histo_Aft.fill(Hit.x(), Hit.y());
					else if (region ==1 )DC2Histo_Aft.fill(Hit.x(), Hit.y());
					else if (region ==2 )DC3Histo_Aft.fill(Hit.x(), Hit.y());
						
						return true;
					}
					else {
						if (region ==0 )DC1Histo_Pre.fill(Hit.x(), Hit.y());
						else if (region ==1 )DC2Histo_Pre.fill(Hit.x(), Hit.y());
						else if (region ==2 )DC3Histo_Pre.fill(Hit.x(), Hit.y());

						return false;
					}
				}
	
	

	
	public ArrayList<H2F> Histograms_Pre(){
		DC1Histo_Pre.setTitle(" DC1 " );DC1Histo_Pre.setTitleX(" X [cm]" );DC1Histo_Pre.setTitleX(" Y [cm]" );
		DC2Histo_Pre.setTitle(" DC2 " );DC2Histo_Pre.setTitleX(" X [cm]" );DC2Histo_Pre.setTitleX(" Y [cm]" );
		DC3Histo_Pre.setTitle(" DC3 " );DC3Histo_Pre.setTitleX(" X [cm]" );DC3Histo_Pre.setTitleX(" Y [cm]" );

		this.Histograms_Cut_Pre.add(DC1Histo_Pre);
		this.Histograms_Cut_Pre.add(DC2Histo_Pre);
		this.Histograms_Cut_Pre.add(DC3Histo_Pre);
		return this.Histograms_Cut_Pre;
	}
	
	public ArrayList<H2F> Histograms_Aft(){
		DC1Histo_Aft.setTitle(" DC1 " );DC1Histo_Aft.setTitleX(" X [cm]" );DC1Histo_Aft.setTitleX(" Y [cm]" );
		DC2Histo_Aft.setTitle(" DC2 " );DC2Histo_Aft.setTitleX(" X [cm]" );DC2Histo_Aft.setTitleX(" Y [cm]" );
		DC3Histo_Aft.setTitle(" DC3 " );DC3Histo_Aft.setTitleX(" X [cm]" );DC3Histo_Aft.setTitleX(" Y [cm]" );
		
		this.Histograms_Cut_Aft.add(DC1Histo_Aft);
		this.Histograms_Cut_Aft.add(DC2Histo_Aft);
		this.Histograms_Cut_Aft.add(DC3Histo_Aft);
		return this.Histograms_Cut_Aft;
	}
	
}
