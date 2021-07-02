package PHD_Clean;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;

import org.jlab.jnp.physics.LorentzVector;
import org.jlab.jnp.physics.Vector3;
import org.jlab.jnp.physics.reaction.TransMatrix;
import org.jlab.groot.data.H1F;
import org.jlab.groot.fitter.DataFitter;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.math.F1D;
import org.jlab.groot.ui.TGCanvas;



public class Pi0_Particle {
	private	double Sector1_B=-30,Sector1_E=30,Sector2_B=30, Sector2_E=90, Sector3_B=90, Sector3_E=150, Sector4_B=150, Sector4_E=-150, Sector5_B=-150,Sector5_E=-90, Sector6_B=-90, Sector6_E=-30;

	private int PhiB0= -180,PhiB1= -150; private int PhiB2= -120;  private int PhiB3= -90;  private int PhiB4= -60; private int PhiB5= -30; private int PhiB6= 0;private int PhiB7= 30; private int PhiB8= 60; private int PhiB9= 90; private int PhiB10= 120;  private int PhiB11= 150; private int PhiB12= 180; 
private double minMass = 0.01;
private double maxMass= 0.015;
	private double Phi_LAB=99;
	private int nrPi0=0;
	
	private double Xf=0;
	
	
	
	static ArrayList<LorentzVector> Photons4Vectors = new ArrayList<>();
	static ArrayList<Integer> Pi0_sector = new ArrayList<Integer> ();
	static ArrayList<Double> Pi0_Z = new ArrayList<Double> ();
	static ArrayList<Double> Pi0_Pt = new ArrayList<Double> ();
	static ArrayList<Double> Pi0_Phi = new ArrayList<Double> ();
	static ArrayList<Double> Pi0_Theta = new ArrayList<Double> ();
	static ArrayList<Double> Pi0_Momentum = new ArrayList<Double> ();
	static ArrayList<Double> Pi0_Energy = new ArrayList<Double> ();
	static ArrayList<Double> Pi0_PhiCLAS12 = new ArrayList<Double> ();
	static ArrayList<Double> Pi0_Mass = new ArrayList<Double> ();
	static ArrayList<Double> Pi0_XMass = new ArrayList<Double> ();
	static ArrayList<Double> Pi0_XF = new ArrayList<Double> ();
	static ArrayList<Double> Pi0_Eta2 = new ArrayList<Double> ();
	static ArrayList<Double> Pi0_EtaBar = new ArrayList<Double> ();
	

	
	
	static ArrayList<LorentzVector> Pi0_4Vector = new ArrayList<LorentzVector> ();
	static ArrayList<LorentzVector> Pi0_X4Vector = new ArrayList<LorentzVector> ();
	
private double minPhi=999;
private double maxPhi=-999;

	int Method =0 ; 
	double EnergyCut=0;
	 public int N0_B1=0,N0_B2=0,N0_B3=0,N0_B4=0,N0_B5=0,N0_B6=0,N0_B7=0,N0_B8=0,N0_B9=0,N0_B10=0,N0_B11=0,N0_B0=0;
	 public int N1_B1=0,N1_B2=0,N1_B3=0,N1_B4=0,N1_B5=0,N1_B6=0,N1_B7=0,N1_B8=0,N1_B9=0,N1_B10=0,N1_B11=0,N1_B0=0;

	public void setMethod(int k) {
		if (k<=2) {
			this.Method= k;}
		else this.Method=0;
	}
	public void setEnergyCut(double En) {
		this.EnergyCut= En;
	}

	public  Pi0_Particle() {
		this.Method=0;
		this.EnergyCut=0;
	}

	//	public void Neutral_Pion_Analysis(ElectronSelection electron, double Degrees) {
	    public void getPi0s(LorentzVector beam, LorentzVector target, LorentzVector vecE,LorentzVector vecQ2, LorentzVector vecW2, double xB, double Q2, double Degrees) {
		this.Phi_LAB=99;
		this.nrPi0=0;
		this.Pi0_Z.clear(); this.Pi0_Pt.clear(); this.Pi0_Phi.clear(); this.Pi0_Mass.clear(); this.Pi0_4Vector.clear();
this.Pi0_sector.clear();
	this.Pi0_Theta.clear();
	this.Pi0_PhiCLAS12.clear();
		this.Pi0_XMass.clear();
		this.Pi0_XF.clear();
	this.Pi0_Eta2.clear();
	this.Pi0_EtaBar.clear();

		double Angl_Cut = Math.toRadians(Degrees);
		boolean AngleCheck=true; 
		double Pq = -target.px()*vecQ2.px()-target.py()*vecQ2.py()-target.pz()*vecQ2.pz()+target.e()*vecQ2.e();
		LorentzVector q = new LorentzVector(0,0,0,0);
		double modulus =  (Math.sqrt(vecQ2.px()*vecQ2.px()+vecQ2.py()*vecQ2.py()+vecQ2.pz()*vecQ2.pz()));
		q.setPxPyPzE(vecQ2.px()/modulus, vecQ2.py()/modulus,vecQ2.pz()/modulus, vecQ2.e());
		if(Method==0 || EnergyCut==0 ) { System.out.println(" Attention the Pi0 extrapolation method should be set and the Energy min cut should be set");return;}
		else if( Method==1) {
			double sector1 =0, sector2=0;
			double PPh=0, Z_var=0,Angolo=0, AngleAA=0,AngleSin=0;
			double Pi0M=0,Pls = 0, Pt=0, Plod = 0;
			double amp2=0.88035;
			//System.out.println( "   -> Photon4Vector size "+Photons4Vectors.size() );
			for( int kk=0 ; kk<Photons4Vectors.size()-1; kk ++)
			{
				sector1=getSector(Photons4Vectors.get(kk).vect()); 
				for(int qq=kk+1 ; qq< Photons4Vectors.size(); qq++) 
				{   
					//System.out.println(" I fotoni li ho ");
					sector2= getSector(Photons4Vectors.get(qq).vect()); // check the two photons are in the same sector
					double angolo = (Photons4Vectors.get(kk).vect().dot(Photons4Vectors.get(qq).vect()))/(Photons4Vectors.get(kk).vect().mag()*Photons4Vectors.get(qq).vect().mag()) ;
					//System.out.println("Angolo Cut " + Math.cos(Angl_Cut) + " Angolo between photon case 2 : " + angolo );
					
					//Angle smalle than 5 is cos close to 1, so I want cosine to be smaller than the maximum cosine 
					if(angolo >=Math.cos(Angl_Cut) ) {
						double angolo2 = (Photons4Vectors.get(kk).vect().dot(Photons4Vectors.get(qq).vect()))/(Photons4Vectors.get(kk).vect().mag()*Photons4Vectors.get(qq).vect().mag()) ;
						//System.out.println("Angolo Cut " + Math.cos(Angl_Cut) + " Angolo between photon case 2 : " + angolo2 );
						AngleCheck=false;}
					else if(angolo>= Math.cos(Angl_Cut)){
						double angolo2 = (Photons4Vectors.get(qq).vect().dot(Photons4Vectors.get(kk).vect()))/(Photons4Vectors.get(kk).vect().mag()*Photons4Vectors.get(qq).vect().mag());
						//System.out.println("Angolo Cut " + Angl_Cut + " Angolo between photon case 2 : " + angolo2 );
						AngleCheck=false;
						}
					if((vecE.vect().dot(Photons4Vectors.get(qq).vect()))/(vecE.vect().mag()*Photons4Vectors.get(qq).vect().mag())>=Math.cos(Angl_Cut) ||(vecE.vect().dot(Photons4Vectors.get(kk).vect()))/(vecE.vect().mag()*Photons4Vectors.get(kk).vect().mag())>=Math.cos(Angl_Cut) ) {
						AngleCheck=false;
			//		System.out.println("Radiated photon "+(vecE.vect().dot(Photons4Vectors.get(qq).vect()))/(vecE.vect().mag()*Photons4Vectors.get(qq).vect().mag()));
					}
				//	System.out.println(AngleCheck);
					if(Photons4Vectors.get(kk).e()>EnergyCut  && Photons4Vectors.get(qq).e()>EnergyCut && AngleCheck==true && sector1== sector2 ) {
						
						this.nrPi0++;

						//System.out.println("I am here !!! Hiii! ");
						LorentzVector PHM = new LorentzVector(0,0,0,0);
						PHM.add(Photons4Vectors.get(kk));
						PHM.add(Photons4Vectors.get(qq));
						PPh = -target.px()*PHM.px()-target.py()*PHM.py()-target.pz()*PHM.pz()+target.e()*PHM.e();
					//	double Z_new=PHM.e()/(beam.e()-vecE.e());
						Pi0_Momentum.add(Math.sqrt(PHM.px()*PHM.px()+PHM.py()*PHM.py()+PHM.pz()*PHM.pz()));
						Pi0_Energy.add(PHM.e());
						
						Z_var= PPh/Pq;
						LorentzVector Vphotons = new LorentzVector(0,0,0,0);
						LorentzVector Gammacm = new LorentzVector(0,0,0,0);
						LorentzVector Protoncm = new LorentzVector(0,0,0,0);
						LorentzVector PIcm = new LorentzVector(0,0,0,0); 
						LorentzVector Beamcm = new LorentzVector(0,0,0,0);
						LorentzVector X = new LorentzVector(0,0,0,0);X.add(target).add(beam).sub(vecE).sub(PHM);
						LorentzVector VGS = new LorentzVector(0,0,0,0); VGS.add(beam).sub(vecE); //THe proper virtual photon
						Vector3 VectorPh = VGS.vect(); // Virtual photon is beam-el scattered but adding the target doesent change tri momenta
						Vector3 VectorPi = PHM.vect();
						Vector3 VectorP = target.vect();
						Vector3 VectorB = beam.vect();
						Vector3 VectorPhNorm = new Vector3(VectorPh.x()/VectorPh.mag(),VectorPh.y()/VectorPh.mag(),VectorPh.z()/VectorPh.mag());				
						Vector3 Angl1 = VectorPhNorm.cross(VectorB);
						Vector3 Angl2 = VectorPhNorm.cross(VectorPi);
						Angl1.setXYZ(Angl1.x()/Angl1.mag(), Angl1.y()/Angl1.mag(),Angl1.z()/Angl1.mag());
						Angl2.setXYZ(Angl2.x()/Angl2.mag(), Angl2.y()/Angl2.mag(),Angl2.z()/Angl2.mag());
						double AngoloCos=Angl1.dot(Angl2);
						Angolo= Math.toDegrees(Math.acos(AngoloCos));
						//this.Phi_LAB = Math.toDegrees(Math.acos(Angolo/(Angl1.mag()*Angl2.mag())));
						Vector3 AngleA = VectorB.cross(VectorPi);
						AngleAA = AngleA.dot(VectorPhNorm);
						Vector3 AngleB=VectorPhNorm.cross(VectorB);
						Vector3 AngleBB =VectorPhNorm.cross(VectorPi);
						AngleSin= AngleAA/(AngleB.mag()*AngleBB.mag());
						double AlternativePhi= Math.toDegrees(Math.asin(AngleSin));
								
								this.Phi_LAB = Angolo;
								
								if(Angolo>0 && AlternativePhi>0) this.Phi_LAB=Angolo;
								else if(Angolo<0 && AlternativePhi>0) this.Phi_LAB=Angolo;
								else if (Angolo<0  && AlternativePhi<0) this.Phi_LAB=Angolo+90;
								else if (Angolo>0 && AlternativePhi<0 ) this.Phi_LAB=360-Angolo;
								
								
								//System.out.println("Angolo Ph e' " + this.Phi_LAB );
								// Nuovo Angolo ORLANDO
								//double AnglePHI_Orlando1= (VectorPh.cross(VectorB)).dot(VectorPi);
								//AnglePHI_Orlando1 = AnglePHI_Orlando1/Math.abs(AnglePHI_Orlando1);
								//System.out.println( "Angolo Orlando " + (AnglePHI_Orlando1*Angolo));
								Vector3 rotationVector = new Vector3(0,0,0);
								TransMatrix rotationMatrix = new TransMatrix();
								rotationMatrix.compose(VectorPh); // Rotating in the direction of virtual photon 
								Vector3 PhCM = rotationMatrix.mult(VectorPh); // Photon rotated
								Vector3 PiCM = rotationMatrix.mult(VectorPi);
								Vector3 PCM = rotationMatrix.mult(VectorP);
								Vector3 BCM= rotationMatrix.mult(VectorB);
								Gammacm.setPxPyPzE(PhCM.x(), PhCM.y(), PhCM.z(), VGS.e());
								//Check these energie! 
								Protoncm.setPxPyPzE(PCM.x(), PCM.y(), PCM.z(), target.e());
								PIcm.setPxPyPzE(PiCM.x(), PiCM.y(), PiCM.z(), PHM.e());
								Beamcm.setPxPyPzE(BCM.x(), BCM.y(), BCM.z(), beam.e());

								Vphotons.add(Gammacm).add(Protoncm);
								Vector3 boost = Vphotons.boostVector();
								boost.setXYZ(-boost.x(), -boost.y(), -boost.z());
								Vphotons.boost(boost);
								Gammacm.boost(boost);
								Protoncm.boost(boost);
								PIcm.boost(boost);
								Beamcm.boost(boost);
								
								 Pi0M= PHM.mass();
								 Pls = PIcm.vect().dot(Gammacm.vect())/Gammacm.vect().mag();
								 Pt = PIcm.vect().mag2()-Pls*Pls; // Questo e' PT^2 
								 Plod = PHM.vect().dot(VGS.vect())/VGS.vect().mag();
								
								Pi0_sector.add(	getSector(VectorPi));
								 Pi0_Z.add(Z_var);
								 Pi0_Pt.add(Pt);
								 Pi0_Phi.add(this.Phi_LAB);
								 Pi0_4Vector.add(PHM);
								 Pi0_Mass.add(PHM.mass());
								 Pi0_X4Vector.add(X);
								 Pi0_Theta.add(Math.toDegrees(VectorPi.theta()));
								Pi0_XMass.add(X.mass());
								Pi0_PhiCLAS12.add(Math.toDegrees(VectorPi.phi()));

		
								 double myXf= 2*Pls/vecW2.mass();
								
								/* Cosa di Harut 
								 * */
													
								 double ehs, pcm2;
								 double pts;
								 double mpro2 = target.e()*target.e();
								 double mhad2 = PHM.mass2();
								 double anu = beam.e() - vecE.e();
								 double ggenzi = PHM.e()/anu;

								 double qq2 = -VGS.mass2();
								 double ww2 = vecW2.mass2();
								 double mx2 = X.mass2();
								 double gxb = -VGS.mass2()/(2*0.938272*VGS.e());
								 double gq2 = qq2;

								 if(ww2 > 0.0){
								   ehs = (ww2+mhad2-mx2)/Math.sqrt(ww2)/2.0;
								   pcm2=(0.25*Math.pow((ww2-mpro2+qq2),2) + mpro2*qq2)/ww2;
								   
								 }
								 else{
								   ehs  = 0.0;
								   pcm2 = 0.0;
								 }

								 double pcm = Math.sqrt(pcm2);
								 double pls = (PHM.e()*0.93827-ehs*Math.sqrt(pcm2+mpro2))/pcm;

								 double pts2 = ehs*ehs - mhad2 - pls*pls;
								 if(pts2 > 0){
								   pts = Math.sqrt(pts2);
								 }
								 else{
								   pts = 0.0;
								 }

								 double xfhad = 2*pls/Math.sqrt(ww2);
								 double et2 = 0.5 * Math.log((ehs+pls)/(ehs-pls));

								 // define Breit frame eta
								 double g2 = 4*gxb*gxb*0.88035/gq2;
								 double xn = 2*gxb/(1.0 + Math.sqrt(1 + g2));
								 double g24 = xn*xn*0.88035/gq2;
								 double yhc1 = xn*xn*0.88035 + xn*gq2;
								 double yhc2 = (1.0 - xn)*gq2;
								 double yhc = Math.log(Math.sqrt(yhc1/yhc2));
								 double etabr = -et2 - yhc;
								
									Pi0_XF.add(xfhad);
									Pi0_Eta2.add(et2);
									Pi0_EtaBar.add(etabr);
								
						
					}					
				}
			}
		}
		
		else if( Method==2) {
			//TO BE FIXED TO BE IMPLEMENTED
			double Energy_max=0;
			LorentzVector Reference = new LorentzVector();
			LorentzVector MaxPh= new LorentzVector();
			LorentzVector SndPh= new LorentzVector();
			int Indice=-1;	
			int F1_sect=0;
			int F2_sect=0;
			// No vertex correction 
			for (int ll=0;ll < Photons4Vectors.size();ll++)
			{
				if(Photons4Vectors.get(ll).e() > Energy_max) { 
					Energy_max=Photons4Vectors.get(ll).e();
					MaxPh.copy(Photons4Vectors.get(ll));
					Indice=ll;
				}
			}
			Energy_max=0;
			for (int ii=0;ii<Photons4Vectors.size();ii++)
			{
				if(ii!=Indice)
				{	
					if( Photons4Vectors.get(ii).vect().dot(MaxPh.vect()) <= Angl_Cut ) AngleCheck=false;
					else if( MaxPh.vect().dot(Photons4Vectors.get(ii).vect()) <= Angl_Cut ) AngleCheck=false;
					if(Photons4Vectors.get(ii).e() > Energy_max &&Photons4Vectors.get(ii)!=Photons4Vectors.get(Indice) && AngleCheck==true) { 
						Energy_max=Photons4Vectors.get(ii).e();
						SndPh.copy(Photons4Vectors.get(ii));
					}
				}
			}
			Vector3 fotone1 = new Vector3();
			Vector3 fotone2 = new Vector3();
			// fotone1.setXYZ(MaxPh.vect().x(), MaxPh.vect().y(), MaxPh.vect().z()-e_vz);
			//  fotone2.setXYZ(SndPh.vect().x(), SndPh.vect().y(), SndPh.vect().z()-e_vz);
			fotone1.setXYZ(MaxPh.vect().x(), MaxPh.vect().y(), MaxPh.vect().z());
			fotone2.setXYZ(SndPh.vect().x(), SndPh.vect().y(), SndPh.vect().z());
			double Angle_ph=(fotone1.dot(fotone2))/(fotone1.mag()*fotone2.mag());
			//double Angle_ph=((MaxPh.vect()).dot(SndPh.vect()))/(MaxPh.vect().mag()*SndPh.vect().mag());
			LorentzVector PHM = new LorentzVector(0,0,0,0);
			PHM.add(MaxPh);
			PHM.add(SndPh);		
			double PPh = -target.px()*PHM.px()-target.py()*PHM.py()-target.pz()*PHM.pz()+target.e()*PHM.e();
			double Z_var= PPh/Pq;
			LorentzVector Vphotons = new LorentzVector(0,0,0,0);

			//check normalization+
			//Vphotons.setPxPyPzE(vecW2.px()/vecW2.mass(), vecW2.py()/vecW2.mass(),vecW2.pz()/vecW2.mass(), vecW2.e()/vecW2.mass());
			double Scalar_old = PHM.vect().dot(vecW2.vect())/vecW2.vect().mag();
			double Pt_old=Math.sqrt(PHM.vect().mag2()-Scalar_old*Scalar_old);

			//Z_Plot.fill(Z_var);

			if( MaxPh.e()> EnergyCut && SndPh.e()>EnergyCut ) {
				LorentzVector PHM2 = new LorentzVector(0,0,0,0);
				PHM2.add(MaxPh);
				PHM2.add(SndPh);		
				LorentzVector Vphotons2 = new LorentzVector(0,0,0,0);
				LorentzVector Gammacm = new LorentzVector(0,0,0,0);
				LorentzVector Protoncm = new LorentzVector(0,0,0,0);
				LorentzVector PIcm = new LorentzVector(0,0,0,0); 
				LorentzVector X= new LorentzVector(0,0,0,0);X.add(target).add(beam).sub(vecE).sub(PHM);
				LorentzVector VGS = new LorentzVector(0,0,0,0); VGS.add(beam).sub(vecE); //THe proper virtual photon
				//System.out.println("VecW2 x,y,z: " +vecW2.vect().x()+" . "+vecW2.vect().y()+" , "+  vecW2.vect().z());
				//	System.out.println("VecVGS x,y,z: " +VGS.vect().x()+" . "+VGS.vect().y()+" , "+  VGS.vect().z());
				Vector3 VectorPh = VGS.vect(); // Virtual photon is beam-el scattered but adding the target doesent change tri momenta
				Vector3 VectorPi=PHM2.vect();
				Vector3 VectorP = target.vect();
				Vector3 VectorB=beam.vect();
				Vector3 VectorPhNorm = new Vector3(VectorPh.x()/VectorPh.mag(),VectorPh.y()/VectorPh.mag(),VectorPh.z()/VectorPh.mag());
				Vector3 Angl1=  VectorPhNorm.cross(VectorB)	;
				Vector3 Angl2=  VectorPhNorm.cross(VectorPi)	;	
				double Angolo=  Angl1.dot(Angl2);
				this.Phi_LAB = Math.toDegrees(Math.acos(Angolo/(Angl1.mag()*Angl2.mag())));

				Vector3 rotationVector = new Vector3(0,0,0);
				TransMatrix rotationMatrix = new TransMatrix();
				rotationMatrix.compose(VectorPh);
				Vector3  PhCM= rotationMatrix.mult(VectorPh);
				Vector3 PiCM= rotationMatrix.mult(VectorPi);
				Vector3 PCM = rotationMatrix.mult(VectorP);
				Gammacm.setPxPyPzE(PhCM.x(), PhCM.y(), PhCM.z(), VGS.e());
				Protoncm.setPxPyPzE(PCM.x(), PCM.y(), PCM.z(), target.e());
				PIcm.setPxPyPzE(PiCM.x(), PiCM.y(), PiCM.z(), PHM.e());
				Vphotons.add(Gammacm).add(Protoncm);
				Vector3 boost = Vphotons.boostVector();
				boost.negative();
				Vphotons.boost(boost);
				Gammacm.boost(boost);
				Protoncm.boost(boost);
				PIcm.boost(boost);
				double Pi0M= PHM2.mass();
				double Pls = PIcm.vect().dot(Gammacm.vect())/Gammacm.vect().mag();
				double Pt=Math.sqrt(PIcm.vect().mag2()-Pls*Pls);
				double Plod = PHM2.vect().dot(VGS.vect())/VGS.vect().mag();
				Pi0_Z.add(Z_var);
				Pi0_Pt.add(Pt);
				Pi0_Phi.add(Angolo);
				Pi0_4Vector.add(PHM);
				Pi0_Mass.add(PHM.mass());
				
				this.nrPi0++;
				//Pi_0_M2_4_BK.fill(PHM2.mass());
				//Pi_0_M2_4.fill(PHM2.mass());
				if(Math.toDegrees(fotone1.phi())<Sector1_E && Math.toDegrees(fotone1.phi())>=Sector1_B) F1_sect=1;
				else if(Math.toDegrees(fotone1.phi())<Sector2_E && Math.toDegrees(fotone1.phi())>=Sector2_B) F1_sect=2;
				else if(Math.toDegrees(fotone1.phi())<Sector3_E && Math.toDegrees(fotone1.phi())>=Sector3_B) F1_sect=3;
				else if(Math.toDegrees(fotone1.phi())<Sector4_E || Math.toDegrees(fotone1.phi())>=Sector4_B) F1_sect=4;
				else if(Math.toDegrees(fotone1.phi())<Sector5_E && Math.toDegrees(fotone1.phi())>=Sector5_B) F1_sect=5;
				else if(Math.toDegrees(fotone1.phi())<Sector6_E && Math.toDegrees(fotone1.phi())>=Sector6_B) F1_sect=6;    
				if(Math.toDegrees(fotone2.phi())<Sector1_E && Math.toDegrees(fotone2.phi())>=Sector1_B) F2_sect=1;
				else if(Math.toDegrees(fotone2.phi())<Sector2_E && Math.toDegrees(fotone2.phi())>=Sector2_B) F2_sect=2;
				else if(Math.toDegrees(fotone2.phi())<Sector3_E && Math.toDegrees(fotone2.phi())>=Sector3_B) F2_sect=3;
				else if(Math.toDegrees(fotone2.phi())<Sector4_E || Math.toDegrees(fotone2.phi())>=Sector4_B) F2_sect=4;
				else if(Math.toDegrees(fotone2.phi())<Sector5_E && Math.toDegrees(fotone2.phi())>=Sector5_B) F2_sect=5;
				else if(Math.toDegrees(fotone2.phi())<Sector6_E && Math.toDegrees(fotone2.phi())>=Sector6_B) F2_sect=6;    
				//	if(F1_sect==1&F2_sect==1) Pi0M_S1.fill(PHM2.mass());
				//	if(F1_sect==2&F2_sect==2) Pi0M_S2.fill(PHM2.mass());
				//	if(F1_sect==3&F2_sect==3) Pi0M_S3.fill(PHM2.mass());
				//	if(F1_sect==4&F2_sect==4) Pi0M_S4.fill(PHM2.mass());
				//	if(F1_sect==5&F2_sect==5) Pi0M_S5.fill(PHM2.mass());
				//	if(F1_sect==6&F2_sect==6) Pi0M_S6.fill(PHM2.mass());
				if(X.mass()>1.5){
					
					double VPr= (PHM2.vect().dot(q.vect())/q.vect().mag());
					Vector3 P_l = new Vector3(VPr*q.vect().x(),VPr*q.vect().y(),VPr*q.vect().z() ) ;
					Vector3 P_t = new Vector3( PHM2.px()-P_l.x(),PHM2.py()-P_l.y(),PHM2.pz()-P_l.z());
					//  System.out.println(P_t.mag());	
					//	H_p_Pt.fill(P_t.mag());
					//H_p_Z.fill(Z_var);
				}
			}
		}
	}

	public double getPhiLab()
	{
		return this.Phi_LAB;
	}
	public void add(Photons photons)
	{
		for(int i =0  ; i< photons.getPhoton4Vects().size();i++) {
			Photons4Vectors.add(photons.getPhoton4Vects().get(i));
		}
	}
	public int getNrPi0s()
	{
		return this.nrPi0;
	}
	public int getNrPhotons()
	{
		return this.Photons4Vectors.size();
	}
public ArrayList<LorentzVector> getLorentzVector(){
	return this.Pi0_4Vector;
}

public ArrayList<LorentzVector> getInvMassLVector(){
	return this.Pi0_X4Vector;
}
public ArrayList<Double> getZ()
{
	return this.Pi0_Z;
}
public ArrayList<Double> getPt()
{
	return this.Pi0_Pt;
}
 
public ArrayList<Double> getTheta()
{
	return this.Pi0_Theta;
}

public ArrayList<Double> getMomentum()
{
	return this.Pi0_Momentum;
}

public ArrayList<Double> getEnergy()
{
	return this.Pi0_Energy;
}

public ArrayList<Double> getXF()
{
	return this.Pi0_XF;
}

public int GetIndex_MostEnergy()
{
	int indice = -9;
	double Energy =0;
	
	for (int i =0 ; i< this.Pi0_4Vector.size(); i++)
	{
		if( this.Pi0_4Vector.get(i).e()> Energy) {
			Energy=this.Pi0_4Vector.get(i).e();
			indice=i; 
		}
	}
	return indice;
}

public  ArrayList<Double> getMissMass(){
	return this.Pi0_XMass;
}


public  ArrayList<Double> getPhiClas12() {
	return this.Pi0_PhiCLAS12;
}

public  ArrayList<Double> getEtaBar() {
	return this.Pi0_EtaBar;
}


public ArrayList<Integer> getSector(){
	return this.Pi0_sector;

}
/**
 * 
 * @return List of pi0s angle in the event in radiant measured from lab frame
 * */
public ArrayList<Double> getPhi()
{
	return this.Pi0_Phi;
}

public ArrayList<Double> getMass()
{
	return this.Pi0_Mass;
}
	

	public void resetPhotons()
	{
		this.Photons4Vectors.clear();
	}
	public void printsHCounts() {
		System.out.println(" Angle min and max" + this.minPhi  +" " +this.maxPhi );
		System.out.println("H0 B1 " + N0_B1 + " H1 B1 " + N1_B1);
		System.out.println("H0 B2 " + N0_B2 + " H1 B2 " + N1_B2);
		System.out.println("H0 B3 " + N0_B3 + " H1 B3 " + N1_B3);
		System.out.println("H0 B4 " + N0_B4 + " H1 B4 " + N1_B4);
		System.out.println("H0 B5 " + N0_B5 + " H1 B5 " + N1_B5);
		System.out.println("H0 B6 " + N0_B6 + " H1 B6 " + N1_B6);
		System.out.println("H0 B7 " + N0_B7 + " H1 B7 " + N1_B7);
		System.out.println("H0 B8 " + N0_B8 + " H1 B8 " + N1_B8);
		System.out.println("H0 B9 " + N0_B9 + " H1 B9 " + N1_B9);
		System.out.println("H0 B10 " + N0_B10 + " H1 B10 " + N1_B10);
		System.out.println("H0 B11 " + N0_B11 + " H1 B11 " + N1_B11);
		// TODO Auto-generated method stub
		
	}
	
	public int getSector(Vector3 Hit) {
	//	System.out.println (" Calcolo dell' hit angolo" + Math.toDegrees(Hit.phi()));
		// remove 30 to put sector 1 on the O, add 180 to transform everything positively.
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

	
}
