package PHD_Clean;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;

public class ReadParameters {

	private static String inputdir,outputdir;//in and out locations
	private static int MCB,PionID,filterBoolean; // basic settings for the analysis 
    private static int Current;// Defining the current to analyze 
	private static int valueBinStatus; // 0: Differential in Q2 , 1 : Differential in xB, 2: Differential in z ...
	private static double W2Cut,Q2Cut,YCut; //  Basic Kinematic cuts
	private static double LowQ2,HighQ2,LowxB,HighxB,LowZ,HighZ; // Limits of Kinematics
	private static int PhiCut,PhiCutBoolean; // Define eventual shaving Phi cuts
	private static double E_Mom_Min,El_min_Angle,El_max_Angle,Vertex_Max,Vertex_Min;//El basic cuts 
	private static double Min_H_momentum,Max_H_momentum,Min_H_theta,Max_H_theta,H_vertexcut;//Hadron basic cuts
	private static double Cut_MissingMass_H; //ehX missing mass cuts 
	private static int ElStatus; //0 no selection - 1 select Status Forward
	private static int HadronStatus; //0 no selection - 1 select Status Forward
    private static int el_FCON,Hadron_FCON; // 1 means apply Fiducial cuts on electron - 0 don't apply fiducial cuts
    private static int el_DCCUT,Hadron_DCCUT;//Dirft Chambers
	private static double el_SF_sigma,el_MinEnCut;//El calorimeter
    private static int el_PCAL_Border,el_ECAL_Border;//lu,lv,lw Cuts
    private static int ElDCRoutine,HadronDCRoutine; // 0 for loading polar cuts, 1 for loading XYZ cuts in local cooridnate system. 
    private static int Chi2PIDBoolean;//0: No chi2Pid Cut | 1: Yes Chi2PID cut
    private static double Chi2PIDSigma;// How many sigma to cut on Chi2PID
    private static int PhotonMethod;//The methods for selecting pi0 from photon
    private static double PhotonSigma;//The number of simga for Pi0 peak selection
    private static int PhotonBackground,Pi0CombBoolean; // Bkg shape and Type of Pi0Analysis from photons 
    private static double PhotonEnergyCut;//Min energy to select photon 
    private static double Ph_el_angle; // Min angle from electron to select photonb 
    private static int Debug; // Debug parameter 
    private static int Resolution; // To set as boolean in case I want to compute resolutions of MC simulation
    private static int Hadr_TrackNDF_min,Hadr_TrackNDF_max ; // Min and Max value of the track NDF
    private static int El_Vert_Sigma,Hadr_Vert_Sigma; // Number of sigma for the vertex resolution cut 
    
	//  public static void main(String args[]) throws FileNotFoundException {
    ReadParameters(String Location) throws FileNotFoundException{
    	  //  File text = new File("/Users/gangelini/testParameters.txt");	
    	  File text = new File(Location);	     
          Scanner scnr = new Scanner(text);	     
	        while(scnr.hasNextLine()){
	            String line = scnr.nextLine();
	            String[] Splitted= line.split(":");	            
	            if(Splitted[0].contentEquals("Inputdir")) inputdir=Splitted[1];
	            else if(Splitted[0].contentEquals("Outputdir")) outputdir=Splitted[1];
	            else if(Splitted[0].contentEquals("Current")) Current=Integer.parseInt(Splitted[1]);
	            else if(Splitted[0].contentEquals("Debug")) Debug=Integer.parseInt(Splitted[1]);
	            else if(Splitted[0].contentEquals("MC_Boolean")) MCB=Integer.parseInt(Splitted[1]);
	            else if(Splitted[0].contentEquals("MC_Resolution")) Resolution=Integer.parseInt(Splitted[1]);
	            else if(Splitted[0].contentEquals("PID")) PionID=Integer.parseInt(Splitted[1]);
	            else if(Splitted[0].contentEquals("Filter")) filterBoolean=Integer.parseInt(Splitted[1]);
	            else if(Splitted[0].contentEquals("Dimension")) valueBinStatus=Integer.parseInt(Splitted[1]);
	            else if(Splitted[0].contentEquals("W2Cut")) W2Cut=Double.parseDouble(Splitted[1]);
	            else if(Splitted[0].contentEquals("Q2Cut")) Q2Cut=Double.parseDouble(Splitted[1]);
	            else if(Splitted[0].contentEquals("YCut")) YCut=Double.parseDouble(Splitted[1]);
	            else if(Splitted[0].contentEquals("Q2_low_limit")) LowQ2=Double.parseDouble(Splitted[1]);
	            else if(Splitted[0].contentEquals("Q2_max_limit")) HighQ2=Double.parseDouble(Splitted[1]);
	            else if(Splitted[0].contentEquals("xB_low_limit")) LowxB=Double.parseDouble(Splitted[1]);
	            else if(Splitted[0].contentEquals("xB_max_limit")) HighxB=Double.parseDouble(Splitted[1]);
	            else if(Splitted[0].contentEquals("Z_low_limit"))  LowZ=Double.parseDouble(Splitted[1]);
	            else if(Splitted[0].contentEquals("Z_max_limit"))  HighZ=Double.parseDouble(Splitted[1]);
	            else if(Splitted[0].contentEquals("Phi_Shaved"))  PhiCutBoolean=Integer.parseInt(Splitted[1]);
	            else if(Splitted[0].contentEquals("Phi_Shaved_cut"))  PhiCut=Integer.parseInt(Splitted[1]);
	            else if(Splitted[0].contentEquals("El_min_Theta"))  El_min_Angle=Double.parseDouble(Splitted[1]);
	            else if(Splitted[0].contentEquals("El_max_Theta"))  El_max_Angle=Double.parseDouble(Splitted[1]);
	            else if(Splitted[0].contentEquals("El_min_vertex"))  Vertex_Min=Double.parseDouble(Splitted[1]);
	            else if(Splitted[0].contentEquals("El_max_vertex"))  Vertex_Max=Double.parseDouble(Splitted[1]);
	            else if(Splitted[0].contentEquals("El_status"))  { 
	            						if(Splitted[1].contentEquals("Forward")) ElStatus=1;
	            						else ElStatus=0;
	            		}
	            else if(Splitted[0].contentEquals("El_FC_Boolean"))  el_FCON=Integer.parseInt(Splitted[1]);
	            else if(Splitted[0].contentEquals("SF_Cut_sigma"))  el_SF_sigma=Double.parseDouble(Splitted[1]);
	            else if(Splitted[0].contentEquals("CAL_MinEn_Cut"))  el_MinEnCut=Double.parseDouble(Splitted[1]);
	            else if(Splitted[0].contentEquals("PCAL_Border_Cut"))  el_PCAL_Border=Integer.parseInt(Splitted[1]);
	            else if(Splitted[0].contentEquals("ECAL_Border_Cut"))  el_ECAL_Border=Integer.parseInt(Splitted[1]);
	            else if(Splitted[0].contentEquals("DC_Cut_Routine_electron")) {
	            if(Splitted[1].contentEquals("Polar")) ElDCRoutine=0;
	            			else if(Splitted[1].contentEquals("XYZ")) ElDCRoutine=1;
	            				else System.out.println("[X] Attention error with DC_Cut_Routine_electron: unrecognized setting. Chose 'Polar' or 'XYZ' a");
	            		}
	            else if(Splitted[0].contentEquals("DC_Cut_Value_electron"))   el_DCCUT=Integer.parseInt(Splitted[1]);
	            
	            else if(Splitted[0].contentEquals("ElVertexSigma"))   El_Vert_Sigma=(int)Double.parseDouble(Splitted[1]);
	            else if(Splitted[0].contentEquals("HadrVertexSigma"))   Hadr_Vert_Sigma=(int)Double.parseDouble(Splitted[1]);
	            
	            else if(Splitted[0].contentEquals("Hadr_TrackNDFMin"))   Hadr_TrackNDF_min=Integer.parseInt(Splitted[1]);
	            else if(Splitted[0].contentEquals("Hadr_TrackNDFMax"))   Hadr_TrackNDF_max=Integer.parseInt(Splitted[1]);

	       
	            else if(Splitted[0].contentEquals("Hadr_mim_mom"))   Min_H_momentum=Double.parseDouble(Splitted[1]);
	            else if(Splitted[0].contentEquals("Hadr_max_mom"))   Max_H_momentum=Double.parseDouble(Splitted[1]);
	            else if(Splitted[0].contentEquals("Hadr_min_Theta")) Min_H_theta=Double.parseDouble(Splitted[1]);
	            else if(Splitted[0].contentEquals("Hadr_max_Theta")) Max_H_theta=Double.parseDouble(Splitted[1]);
	            else if(Splitted[0].contentEquals("Hadr_vertex"))    H_vertexcut=Double.parseDouble(Splitted[1]);
	            else if(Splitted[0].contentEquals("Hadr_status"))  { 
					if(Splitted[1].contentEquals("Forward")) HadronStatus=1;
					else HadronStatus=0;
	}
	            else if(Splitted[0].contentEquals("el_Had_MX_cut"))    Cut_MissingMass_H=Double.parseDouble(Splitted[1]);
	            else if(Splitted[0].contentEquals("Hadron_FC_Boolean"))  Hadron_FCON=Integer.parseInt(Splitted[1]);
	            else if(Splitted[0].contentEquals("DC_Cut_Routine_hadron")) {
	            		if(Splitted[1].contentEquals("Polar")) HadronDCRoutine=0;
		            			else if(Splitted[1].contentEquals("XYZ")) HadronDCRoutine=1;
		            				else System.out.println("[X] Attention error with DC_Cut_Routine_electron: unrecognized setting. Chose 'Polar' or 'XYZ' a");
		            		}
	            else if(Splitted[0].contentEquals("DC_Cut_Value_hadron"))   Hadron_DCCUT=Integer.parseInt(Splitted[1]);
	            else if(Splitted[0].contentEquals("Chi2PID_Boolean"))   Chi2PIDBoolean=Integer.parseInt(Splitted[1]);
	            else if(Splitted[0].contentEquals("Chi2PID_Sigma"))   Chi2PIDSigma=Double.parseDouble(Splitted[1]);
	            else if(Splitted[0].contentEquals("Photon_method"))  PhotonMethod=Integer.parseInt(Splitted[1]);
	            else if(Splitted[0].contentEquals("Photon_sigma_fit"))  PhotonSigma=Double.parseDouble(Splitted[1]);
	            else if(Splitted[0].contentEquals("Photon_Background"))  PhotonBackground=Integer.parseInt(Splitted[1]);
	            else if(Splitted[0].contentEquals("Photon_Energy_min"))  PhotonEnergyCut=Double.parseDouble(Splitted[1]);
	            else if(Splitted[0].contentEquals("Pi0_Combinatory_Boolean"))  Pi0CombBoolean=Integer.parseInt(Splitted[1]);
	            else if(Splitted[0].contentEquals("Photon_Electron_Angle"))  Ph_el_angle=Double.parseDouble(Splitted[1]);

	        }      
	       System.out.println(" ==== Printing the parameters for the Data Analysis ===== ");
	       System.out.println(" -> Input Dir: "+ inputdir); 
	       System.out.println(" -> Output Dir: "+ outputdir);
	       System.out.println(" -> Current analyzed: "+ Current);
	       System.out.println(" -> MC?: "+ MCB);
	       System.out.println(" -> PID: "+ PionID);
	       System.out.println(" -> Index Analysis: "+ valueBinStatus); 
	       System.out.println(" -> eHX filter? : "+ filterBoolean); 
	       System.out.println(" ==== Kinematic cuts ===== ");
	       System.out.println(" -> Q2 >: "+ Q2Cut);
	       System.out.println(" -> W >: "+ W2Cut);
	       System.out.println(" -> y <: "+ YCut);
	       System.out.println(" ==== Electron cuts ===== ");
	       System.out.println(" -> El min theta: "+ El_min_Angle); 
	       System.out.println(" -> El max theta: "+ El_max_Angle);
	       System.out.println(" -> El vertex between "+ Vertex_Min+ " and "+Vertex_Max);

	       System.out.println(" -> El vertex sigma cut: "+El_Vert_Sigma);
	       System.out.println(" -> El only in Forward? "+ ElStatus);
	       System.out.println(" -> El Fiducial Cuts? "+ el_FCON);
	       System.out.println(" -> Sampling Fraction Sigma "+ el_SF_sigma);
	       System.out.println(" -> Minimum Energy PCAL  "+ el_MinEnCut);
	       System.out.println(" -> PCAL border cuts (none/v.lose/lose/medium/tight/v.tight)  "+ el_PCAL_Border);
	       System.out.println(" -> ECAL border cuts (none/v.lose/lose/medium/tight/v.tight)  "+ el_ECAL_Border);
	       System.out.println(" -> DC Routine (Polar/XYZ)  "+ ElDCRoutine);
	       System.out.println(" -> DC cuts (none/v.lose/lose/medium/tight/v.tight)  "+ el_DCCUT);
	       System.out.println(" -> Hadron Track NDF min "+Hadr_TrackNDF_min);
	       System.out.println(" -> Hadron Track NDF max "+Hadr_TrackNDF_max);
	       System.out.println(" ==== Hadron cuts ===== ");
	       System.out.println(" -> Hadron min momentum: "+ Min_H_momentum); 
	       System.out.println(" -> Hadron max momentumn: " + Max_H_momentum );
	       System.out.println(" -> Hadron min theta: " + Min_H_theta ); 
	       System.out.println(" -> Hadron max theta: "+ Max_H_theta);
	       System.out.println(" -> |Hadron-electron| vertex  range: " + H_vertexcut);
	       System.out.println(" -> Hadron vertex sigma cut: "+Hadr_Vert_Sigma);

	       System.out.println(" -> Missing Mass > : " + Cut_MissingMass_H);
	       System.out.println(" -> Hadron Fiducial Cuts? "+ Hadron_FCON);
	       System.out.println(" -> DC Routine (Polar/XYZ)  "+ HadronDCRoutine);
	       System.out.println(" -> DC cuts (none/v.lose/lose/medium/tight/v.tight)  "+ Hadron_DCCUT);
	       System.out.println(" -> Chi2PID?  "+ Chi2PIDBoolean);
	       System.out.println(" -> Chi2PID sigma:  "+ Chi2PIDSigma);
	       System.out.println(" ==== Photon + Pi0 cuts ===== ");
	       System.out.println(" -> Photon Selection Method : "+ PhotonMethod);
	       System.out.println(" -> Photon Background Function : "+ PhotonBackground);
	       System.out.println(" -> Photon min Energy Cut :  "+ PhotonEnergyCut);
	       System.out.println(" -> Pion  Selection Sigmas :  "+ PhotonSigma);
	       System.out.println(" -> Pion  Combinatory of photons?   :  "+ Pi0CombBoolean);
	       System.out.println(" -> Photons to electtron min angle cut   :  "+ Ph_el_angle);
	      
	       
	       System.out.println(" ==================================== ");
	  }  
	  // ----------------------
	  //    RETURN FUNCTIONS 
	  //------------------------
    
    /**
  	  * Return the Input Dir   
  	  */
  	  public String Get_InputDir(){return inputdir;}  
  	 /**
   	  * Return the Output Dir   
   	  */
   	  public String Get_OutputDir(){return outputdir;}  
   	 /**
    	  * Return the Current we are analyizing   
    	  */
   	  public Integer Get_Current(){return Current;}
   	 /**
	   * Return the 0/1 if DEBUG variable is OFF/ON 

	   */
	  public boolean Get_DebugBoolean() {
		  if(Debug==1)return true;
		  else return false;
	  }
  	  /**
	  /**
	  * Return the Boolean  
	  * To identify if the MC is ON (true) or OFF (false)
	  */
	  public boolean GetMCBoolean(){if(MCB==1)return true; else return false;}  
	  /**
		  * Return the Boolean  
		  * To identify if MC resolution are on or OFF
		  */
		  public boolean GetMCResolutions(){if(MCB==1&&Resolution==1)return true; else return false;}  
	  /**
	   * Return the PID used in the Analysis 
	   */
	  public int GetPID(){return PionID;}
	  /**
	   * Return the Index for analysis (0-5)
	   * 0:1D Q2 | 1: 1D xB | 2: 1D xz | 3: 1D PT | 4 : 2D Q2-X | 5: 4D 
	   */
	  public int GetIndexAnalysis() {return valueBinStatus;}
	  /**
	   * Return 0/1 based if we want the SemiInclusive filter on or off
	   */
	  public int GetFilter() {return filterBoolean;}
	  /**
	   * Return The W Cut
	   */
	  public Double GetW2Cut() {return W2Cut;}
	  /**
	   * Return The Q2 Cut
	   */
	  public Double GetQ2Cut() {return Q2Cut;}
	  /**
	   * Return The Y Cut
	   */
	  public Double GetYCut() {return YCut;}
	  /**
	   * Return The Q2 low Limit 
	   */
	  public Double GetQ2LowLimit() {return LowQ2;}
	  /**
	   * Return The Q2 Upper Limit 
	   */
	  public Double GetQ2HighLimit() {return HighQ2;}
	  /**
	   * Return The xB low Limit 
	   */
	  public Double GetxBLowLimit() {return LowxB;}
	  /**
	   * Return The xB Upper Limit 
	   */
	  public Double GetxBHighLimit() {return HighxB;}
	  /**
	   * Return The Z low Limit 
	   */
	  public Double GetZLowLimit() {return LowZ;}
	  /**
	   * Return The Z Upper Limit 
	   */
	  public Double GetZHighLimit() {return HighZ;}
	  /**
	   * Return 1 if PHI Shaved Cuts are activated
	   */
	  public boolean GetPhiCutBoolean() {if(PhiCutBoolean==1) return true; else return false;}
	  /**
	   * Return the Phi range for phi shaved cuts
	   */
	  public Integer GetPhiCutRange() {return PhiCut;}
	  /**
	   * Return the minimum Theta angle for electron 
	   */
	  public Double Get_el_minTheta() {return El_min_Angle;}
	  /**
	   * Return the maximum Theta angle for electron 
	   */
	  public Double Get_el_maxTheta() {return El_max_Angle;}
	  /**
	   * Return the minimum vertex value electron 
	   */
	  public Double Get_el_minVertex() {return Vertex_Min;}
	  /**
	   * Return the maximum vertex value for electron 
	   */
	  public Double Get_el_maxVertex() {return Vertex_Max;}
	  /**
	   * Return a boolean for electron : TRUE if only Forward dector needs to be consider 
	   */
	  public boolean Get_el_Forward() {if (ElStatus==1) return true; else return false;} 
	  /**
	   * Return a boolean for the hadron: TRUE if only Forward dector needs to be consider 
	   */
	  public boolean Get_Hadron_Forward() {if (HadronStatus==1) return true; else return false;}
	  /**
	   * Return a boolean for Fiducial Cuts of electron : TRUE 
	   */
	  public boolean Get_el_FC() {if (el_FCON==1) return true; else return false;} 
	  /**
	   * Return a boolean for the Fiducial Cuts of  hadron: TRUE 
	   */
	  public boolean Get_Hadron_FC() {if (Hadron_FCON==1) return true; else return false;}
	  /**
	   * Return the min momenta for Hadron
	   */
	  public Double Get_Hadron_minP() {return Min_H_momentum;} 
	  /**
	   * Return the max momenta for Hadron
	   */
	  public Double Get_Hadron_maxP() {return Max_H_momentum; }
	  /**
	   * Return the min Theta for Hadron
	   */
	  public Double Get_Hadron_minTheta() {return Min_H_theta;} 
	  /**
	   * Return the max  Theta for Hadron
	   */
	  public Double Get_Hadron_maxTheta() {return Max_H_theta; }
	  /**
	   * Return the difference between electron and pion vertex 
	   */
	  public Double Get_Hadron_vertexCut() {return H_vertexcut; }
	  /**
	   * Return the Min Missing Mass to cut for eHX 
	   */
	  public Double Get_Hadron_MissingMassCut() {return Cut_MissingMass_H; }
	  /**
	   * Return the electron DC routine (0: Polar | 1 : XYZ) 
	   */
	  public Integer Get_el_DCRoutine() {return ElDCRoutine; }
	  /**
	   * Return the electron DC routine (0: Polar | 1 : XYZ) 
	   */
	  public Integer Get_Hadron_DCRoutine() {return HadronDCRoutine; }
	  /**
	   * Return the electron DC cut 
	   * 0:None | 1: Very Loose | 2: Loose
	   * 3: Medium | 4 : Tight | 5: Very Tight
	   */
	  public Integer Get_el_DCCUT() {return el_DCCUT; }
	  /**
	   * Return the electron DC cut 
	   * 0:None | 1: Very Loose | 2: Loose
	   * 3: Medium | 4 : Tight | 5: Very Tight
	   */
	  public Integer Get_Hadron_DCCUT() {return Hadron_DCCUT; }
	  
	  /**
	   * Return the min cut for track NDF 
	   */
	  public Integer Get_Hadron_TrackNDF_Min() {return Hadr_TrackNDF_min; }
	  /**
	   * Return the max cut for track NDF 
	   */
	  public Integer Get_Hadron_TrackNDF_Max() {return Hadr_TrackNDF_max; }

	  /**
	   * Return the ammount of sigma to cut the vertex resolution for electrons
	   */
	  public Integer Get_el_VertexSigma() {return El_Vert_Sigma;}

	  /**
	   * Return the ammount of sigma to cut the vertex resolution for hadrons
	   */
	  public Integer Get_Hadron_VertexSigma() {return El_Vert_Sigma;}
	  /**
	   * Return the electron SF sigma for cuts 
	   */
	  public Double Get_el_SFSigma() {return el_SF_sigma; }
	  /**
	   * Return the electron Minimum PCAL Energy cut 
	   */
	  public Double Get_el_MinEn() {return el_MinEnCut; }
	  /**
	   * Return the electron PCAL Border Cut
	   * * 0:None | 1: Very Loose | 2: Loose
	   * 3: Medium | 4 : Tight | 5: Very Tight
	   */
	  public Integer Get_el_PCALCut() {return el_PCAL_Border; }
	  /**
	   * Return the electron ECAL Border Cut
	   * * 0:None | 1: Very Loose | 2: Loose
	   * 3: Medium | 4 : Tight | 5: Very Tight
	   */
	  public Integer Get_el_ECALCut() {return el_ECAL_Border; }
	  /**
	   * Return TRUE if the Chi2PID cuts are active 
	   */
	  public boolean Get_Hadron_Ch2PIDBoolean() {if(Chi2PIDBoolean==1)return true;else return false; }
	  /**
	   * Return the Chi2PID sigma for Hadron cuts 
	   */
	  public Double Get_Hadron_Ch2PIDSigma() {return Chi2PIDSigma; }
	  /**
	   * Return The method for photon analysis 
	   */
	  public Integer Get_Gamma_Method() {return PhotonMethod; }
	  /**
	   * Return The number of sigma for Pi0 fits on 2gamma Invariant Mass 
	   */
	  public Double Get_Pi0_Sigma() {return PhotonSigma; }
	  /**
	   * Return an integer rappresenting the function used for Pion subtraction 
	   */
	  public Integer Get_Pi0_BackgroundID() {return PhotonBackground; }
	  /**
	   * Return TRUE if the the Pi0 is done with combinatory method
	   * Return FALSE if the Pi0 is done with the 2 most energetic photons 
	   */
	  public boolean Get_Pi0_CombinatoryBoolean() {if(Pi0CombBoolean==1)return true;else return false; }
	  /**
	   * Return The minimum energy in order to select a good photon 
	   */
	  public Double Get_Gamma_EnergyCut() {return PhotonEnergyCut; }
	  /**
	   * Return The min angle between electron and photon to cut on 
	   */
	  public Double Get_Gamma_elAngle() {return Ph_el_angle; }
	  
} 
