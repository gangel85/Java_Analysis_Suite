package PHD_Clean;

import java.security.Policy.Parameters;
import java.util.ArrayList;
import java.util.List;

import org.jlab.groot.data.H2F;
import org.jlab.jnp.physics.Vector3;

public class EC_Cut implements DetectorCut{
	//Particle that doesn't survive the cuts
	H2F EC1Histo_Pre = new H2F("EC1_Pre",450,0,450,450,0,450);
	H2F EC2Histo_Pre = new H2F("EC2_Pre",450,0,450,450,0,450);
	H2F EC3Histo_Pre = new H2F("EC3_Pre",450,0,450,450,0,450);
	// Particle that does survive the cuts 
	H2F EC1Histo_Aft = new H2F("EC1_Aft",450,0,450,450,0,450);
	H2F EC2Histo_Aft = new H2F("EC2_Aft",450,0,450,450,0,450);
	H2F EC3Histo_Aft = new H2F("EC3_Aft",450,0,450,450,0,450);
    ArrayList<H2F> Histograms_Cut_Pre = new ArrayList<H2F>();
    ArrayList<H2F> Histograms_Cut_Aft = new ArrayList<H2F>();
     private String name;
	// THIS IS MY INTERNAL NOTATION USED IN PARTICLE REC 
    private int PCALHIT=70; // I am taking the x,y,z,on PCAL 
	private int PCALID=71; //The index used to idenfity the PCAL
	private int ECAL1ID=74; //The index used to idenfity the PCAL
	private int ECAL2ID=77; //The index used to idenfity the PCAL
	// ------------------------------------------------
	private double lu_cutMin,lv_cutMin,lw_cutMin,lu_cutMax,lv_cutMax,lw_cutMax;
	// Variables use internally
	private boolean PCAL_cut,EC1_cut,EC2_cut ;
	private Vector3 Hits;

	
	  // Cut definitions
	  double data_min_v_tight_inb[] = {24.0, 24.0, 24.0, 24.0, 24.0, 24.0};
	  double data_min_v_med_inb[]   = {19.0, 19.0, 19.0, 19.0, 19.0, 19.0};
	  double data_min_v_loose_inb[] = {15.0,  15.0,  15.0,  15.0,  15.0,  15.0 };
	  //
	  double data_max_v_tight_inb[] = {400, 400, 400, 400, 400, 400};
	  double data_max_v_med_inb[]   = {400, 400, 400, 400, 400, 400};
	  double data_max_v_loose_inb[] = {400, 400, 400, 400, 400, 400};
	  //
	  double data_min_w_tight_inb[] = {24.0, 24.0, 24.0, 24.0, 24.0, 24.0};
	  double data_min_w_med_inb[]   = {19.0, 19.0, 19.0, 19.0, 19.0, 19.0};
	  double data_min_w_loose_inb[] = {15.0,  15.0,  15.0,  15.0,  15.0,  15.0 };
	  // 
	  double data_max_w_tight_inb[] = {400, 400, 400, 400, 400, 400};
	  double data_max_w_med_inb[]   = {400, 400, 400, 400, 400, 400};
	  double data_max_w_loose_inb[] = {400, 400, 400, 400, 400, 400};
	  
	  // Cut definitions MC
	  double mc_min_v_tight_inb[] = {23.0, 23.0, 23.0, 23.0, 23.0, 23.0};
	  double mc_min_v_med_inb[]   = {18.0, 18.0, 18.0, 18.0, 18.0, 18.0};
	  double mc_min_v_loose_inb[] = {14.0,  14.0,  14.0,  14.0,  14.0,  14.0 };
	  //
	  double mc_max_v_tight_inb[] = {400, 400, 400, 400, 400, 400};
	  double mc_max_v_med_inb[]   = {400, 400, 400, 400, 400, 400};
	  double mc_max_v_loose_inb[] = {400, 400, 400, 400, 400, 400};
	  //
	  double mc_min_w_tight_inb[] = {23.0, 23.0, 23.0, 23.0, 23.0, 23.0};
	  double mc_min_w_med_inb[]   = {18.0, 18.0, 18.0, 18.0, 18.0, 18.0};
	  double mc_min_w_loose_inb[] = {14.0,  15.0,  15.0,  15.0,  15.0,  15.0 };
	  // 
	  double mc_max_w_tight_inb[] = {400, 400, 400, 400, 400, 400};
	  double mc_max_w_med_inb[]   = {400, 400, 400, 400, 400, 400};
	  double mc_max_w_loose_inb[] = {400, 400, 400, 400, 400, 400};
	  
	  
	  int CutIndex=0;
	  boolean MC;
	EC_Cut(ReadParameters parameters){
		//Implment to some ok values as inizialization
		this.lu_cutMin=19;this.lv_cutMin=19;this.lw_cutMin=19;this.lu_cutMax=400;this.lv_cutMax=400;this.lw_cutMax=400;
		 CutIndex = parameters.Get_el_PCALCut();
		 MC=parameters.GetMCBoolean();
	}
	
	//Constructor with costum limits to apply to all calorimeters
	EC_Cut(double lum, double lvm, double lwm ,double luM, double lvM, double lwM){
		this.lu_cutMin=lum;this.lv_cutMin=lvm;this.lw_cutMin=lwm;this.lu_cutMax=luM;this.lv_cutMax=lvM;this.lw_cutMax=lwM;	
	}



	public void setName(String nameDetector) {
		this.name=nameDetector;
	}
	public String toString() {
		return this.name;
	}
	
	/**
	 * With this function the user can set new cuts for the variable lu,lv,lw
	 * @author gangelini
	 * @param new cuts value
	 * lu min, lv min , lw min, lu max, lv max ,lw max 
	 */
	public void SetNewLimits(double lu_min, double lv_min, double lw_min, double lu_max, double lv_max, double lw_max) {
		this.lu_cutMin=lu_min;this.lv_cutMin=lv_min;this.lw_cutMin=lw_min;this.lu_cutMax=lu_max;this.lv_cutMax=lv_max;this.lw_cutMax=lw_max;
	}
	/**
	 *  This function return the result of the cut 
	 *  case implemented: Electron and Photon
	 *  No other particles are implemented 
	 */
	public boolean Status(ParticleREC particle) {
		if (particle.pid()==11) { 
			return this.ElectronCut(particle) ;
		}
		else if (particle.pid()==22) return this.PhotonCut(particle);
		else return false;
		
	//	}
	}

	// --- Internal operation 
	private boolean PhotonCut(ParticleREC particle) {
		if(particle.hasDetector(PCALID)==true){ 
			Hits=particle.getDetectorHits(PCALID);
			return EC_cut(Hits.x(),Hits.y(),Hits.z());
		}
		else return false;
	}

	private boolean ElectronCut(ParticleREC particle){
		if(particle.hasDetector(PCALID)==true) {
			PCAL_cut=false;
			Hits=particle.getDetectorHits(PCALID);
			int sector = getSector(Hits);
			// Limits
			if(MC==true) {
			if(CutIndex==0 ) return true;
			else if (CutIndex==1)SetNewLimits(0,mc_min_v_loose_inb[sector-1],mc_min_w_loose_inb[sector-1],450,mc_max_v_loose_inb[sector-1],mc_max_w_loose_inb[sector-1]);
			else if (CutIndex==2)SetNewLimits(0,mc_min_v_med_inb[sector-1],mc_min_w_med_inb[sector-1],450,mc_max_v_med_inb[sector-1],mc_max_w_med_inb[sector-1]);
			else if (CutIndex==3)SetNewLimits(0,mc_min_v_tight_inb[sector-1],mc_min_w_tight_inb[sector-1],450,mc_max_v_med_inb[sector-1],mc_max_w_tight_inb[sector-1]);
			}
			else if (MC==false) {
				if(CutIndex==0 ) return true;
				else if (CutIndex==1)SetNewLimits(0,data_min_v_loose_inb[sector-1],data_min_w_loose_inb[sector-1],450,data_max_v_loose_inb[sector-1],data_max_w_loose_inb[sector-1]);
				else if (CutIndex==2)SetNewLimits(0,data_min_v_med_inb[sector-1],data_min_w_med_inb[sector-1],450,data_max_v_med_inb[sector-1],data_max_w_med_inb[sector-1]);
				else if (CutIndex==3)SetNewLimits(0,data_min_v_tight_inb[sector-1],data_min_w_tight_inb[sector-1],450,data_max_v_med_inb[sector-1],data_max_w_tight_inb[sector-1]);
			}
		//	System.out.println(" Hit "+particle.getDetectorHits(71).y() +" "+ particle.getDetectorHits(71).z() + " values "+ lv_cutMin+ " " + lv_cutMax);
			if(particle.getDetectorHits(71).y()>lv_cutMin && lv_cutMax>particle.getDetectorHits(71).y()&& lw_cutMax>particle.getDetectorHits(71).z() && lw_cutMin<particle.getDetectorHits(71).z()	)
			{	
			//System.out.println(" FC EC Sector is " + particle.getSector(PCALID));
				if (particle.getSector(PCALID) ==1  && particle.getDetectorHits(71).z() >70 &&particle.getDetectorHits(71).z() <98 ) return false;
				else if (particle.getSector(PCALID) ==2 && particle.getDetectorHits(71).y() >96 &&particle.getDetectorHits(71).y()<118) return false ;
				else return true;
        	
		}
			/*
			if(particle.getDetectorHits(71).y()>lv_cutMin && lv_cutMax>particle.getDetectorHits(71).y()&& lw_cutMax>particle.getDetectorHits(71).z() && lw_cutMin<particle.getDetectorHits(71).z()	)
			{	
			//	System.out.println(" This is good! ");
            return true;
		}
		*/
		
		}
		return PCAL_cut;
		
	}

	/**
	 * The class takes the lu,lv,lw coordinates from EC banks and apply the cut.
	 * @param d 
	 * @param e
	 * @param f
	 * @return boolean T if pass, F if not passed 
	 */
	private  boolean EC_cut (double d, double e, double f) {
		if (d>=this.lu_cutMin && e >= this.lv_cutMin && f >= this.lw_cutMin && d<= this.lu_cutMax && e<=this.lv_cutMax & f<=this.lw_cutMax) {
			return true;
		}
		else {
			//System.out.println("Che schifo ! the velues of d, e and f are" + d + " " + e  + " " + f );
			return false ; 
		}
	}
	@Override
	/**
	 * This function return the sector of a given hit .
	 * @param Hit is a Vector3 representing x,y,z of a detector reconstructed cluster
	 **/
	
	
		public int getSector(Vector3 Hit) {
		//	System.out.println (" Calcolo dell' hit angolo" + Math.toDegrees(Hit.phi()));
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
			return sec; 
			
			}
	@Override
	public List<H2F> Histograms_Pre() {
		this.Histograms_Cut_Pre.add(EC1Histo_Pre);
		this.Histograms_Cut_Pre.add(EC2Histo_Pre);
		this.Histograms_Cut_Pre.add(EC3Histo_Pre);
		return this.Histograms_Cut_Pre;
		// TODO Auto-generated method stub
	}
	public List<H2F> Histograms_Aft() {
		this.Histograms_Cut_Aft.add(EC1Histo_Aft);
		this.Histograms_Cut_Aft.add(EC2Histo_Aft);
		this.Histograms_Cut_Aft.add(EC3Histo_Aft);
		return this.Histograms_Cut_Aft;
		// TODO Auto-generated method stub
	  }
	}

