package PHD_Clean;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.data.H3F;
import org.jlab.groot.data.TDirectory;
import org.jlab.groot.fitter.DataFitter;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.math.F1D;
import org.jlab.groot.ui.TGCanvas;
import org.jlab.jnp.hipo.data.HipoEvent;
import org.jlab.jnp.hipo.data.HipoGroup;
import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.data.Schema;
import org.jlab.jnp.hipo4.data.Schema.SchemaBuilder;
import org.jlab.jnp.hipo4.io.HipoWriter;
import org.jlab.jnp.physics.LorentzVector;
import org.jlab.jnp.physics.Vector3;
import org.jlab.jnp.physics.reaction.TransMatrix;


public class MCAnalysis {
	double z1,z2,z3,z4,z5,z6,z7,z8,z9,z10;
	double pt1,pt2,pt3,pt4,pt5,pt6,pt7,pt8,pt9,pt10;
	private int nr_particle;
	double beamEnergy = 10.6;
	private LorentzVector VB   = new LorentzVector(0.0,0.0,beamEnergy,beamEnergy); // Getting the Lorentz Vector associated to the beam
	private LorentzVector VT = new LorentzVector(0.0,0.0,0.0,0.93827); // Defining Lorentz Vector for the proton target

	private int parentID;
	private FileWriter fileWriter;
	private  double min_pion= 0.05; private  double max_pion=0.25;	private int Mass_Bins=100;
	private int Z_Bins=0; private double Z_min=0.2; private double Z_max=0.8; 
	private int PT_Bins=0;private double PT_min=0.1; private double PT_max=1; 
	private H2F InvMass,H_e_vx_vy;
	private double Q2_el_MC= 0;
	private double XB_el_MC=0; 
	private double y_el_MC=0;
	private double W_el_MC=0;
	private double MX=0;
	private boolean PolygonalBinClas;

	
	// Variables of the electron:
	
	public double e_xB, e_Q2, e_W, e_theta, e_phi, e_p;
	
   private int CountsID=0;
   private int PID=0;

	private SchemaBuilder schemaBuilder ; 
	private Schema schemaMC;
    private HipoWriter  writeroutputMC = new HipoWriter();

   
  

   
	
	
	
	public ArrayList <MultiBins>  Particles_LUND= new ArrayList<MultiBins>();
	
	private ArrayList<LorentzVector> Particle_4Vector = new ArrayList<>();	
	public ArrayList <Double>  Hadron_en= new ArrayList<Double>();
	public ArrayList <Double>  Hadron_mom= new ArrayList<Double>();
	public ArrayList <Double>  Hadron_theta= new ArrayList<Double>();
	public ArrayList <Double>  Hadron_phi= new ArrayList<Double>();
	public ArrayList <Double> Hadron_phiTrento= new ArrayList<Double>();
	public ArrayList <Double>  Hadron_z= new ArrayList<Double>();
	public ArrayList <Double>  Hadron_pt= new ArrayList<Double>();
	public ArrayList <Integer> Hadron_parent = new ArrayList<Integer>();
	public ArrayList <Integer> Hadron_pid = new ArrayList<Integer>();
	public ArrayList <Double> Hadron_xf = new ArrayList<Double>();
	public ArrayList <Double> Hadron_eta = new ArrayList<Double>();
	public ArrayList <Double> Hadron_mx = new ArrayList<Double>();
	public ArrayList <Double> Hadron_vz = new ArrayList<Double>();

	
	
//TGCanvas CountsPion= new TGCanvas("CountsPion","CountsPion",1600,1400);
	//H1F CountsParentID = new H1F("CountsID",799,-400,400);
	public MCAnalysis( boolean Poly, MultiBins invariantMassBins, int zbin, double zmin, double zmax, int ptbin, double ptmin, double ptmax, int phibins)
	{   

		
		schemaBuilder = new SchemaBuilder("LUND::ePiX",130,1);  
		schemaBuilder.addEntry("event", "L", "Event id");
		schemaBuilder.addEntry("pid", "I", "Pion id");
		schemaBuilder.addEntry("e_p", "D", "Trigger electron momentum");
		schemaBuilder.addEntry("e_theta", "D", "Trigger electron theta");
		schemaBuilder.addEntry("e_phi", "D", "Trigger electron phi ");
		schemaBuilder.addEntry("Q2", "D", "Q2 of the event");
		schemaBuilder.addEntry("W", "D", "W of the event");
		schemaBuilder.addEntry("x", "D", "x of the event");
		schemaBuilder.addEntry("y", "D", "y of the event");
		schemaBuilder.addEntry("pi_p", "D", "Pion momentum");
		schemaBuilder.addEntry("pi_theta", "D", "Pion theta");
		schemaBuilder.addEntry("pi_phi", "D", "Pion phi");
		schemaBuilder.addEntry("z", "D", "x z of the pion");
		schemaBuilder.addEntry("pt", "D", " PT  of the pion (with respect virtual photon)");
		schemaBuilder.addEntry("phi", "D", "Phi (trento) of the pion");
		schemaBuilder.addEntry("xf", "D", "x Feynman of the pion");
		schemaBuilder.addEntry("mX", "D", " Missing mass ehX");
		schemaBuilder.addEntry("eta", "D", "Rapidity in Breit frame");
		schemaBuilder.addEntry("father", "I", "Father of the particle ");
		schemaBuilder.addEntry("granfather", "I", "Granfather of the particle");
		schemaMC  = schemaBuilder.build();
		  writeroutputMC.getSchemaFactory().addSchema(schemaMC); // Should I use schema or factory? 
		  writeroutputMC.open("Generated_Skimmed_Hipo.hipo");
		

	}


	/**
	 * Process another event 
	 * @param event
	 */
	public void add(ReadParameters Parameters, Bank MCParticle, Bank bankMC, int Particle, long evento ) {
		
	
		ArrayList<Double> TuplHadron_momentum = new ArrayList<Double>();
		ArrayList<Double> TuplHadron_theta = new ArrayList<Double>();
		ArrayList<Double> TuplHadron_phiclas = new ArrayList<Double>();
		ArrayList<Double> TuplHadron_z = new ArrayList<Double>();
		ArrayList<Double> TuplHadron_pt = new ArrayList<Double>();
		ArrayList<Double> TuplHadron_phi = new ArrayList<Double>();
		ArrayList<Double> TuplHadron_xf = new ArrayList<Double>();
		ArrayList<Double> TuplHadron_MX = new ArrayList<Double>();
		ArrayList<Double> TuplHadron_eta = new ArrayList<Double>();
		ArrayList<Integer> TuplHadron_father = new ArrayList<Integer>();
		ArrayList<Integer> TuplHadron_granfather = new ArrayList<Integer>();

		
		ResetHadrons();
		Hadron_en.clear();Hadron_mom.clear(); Hadron_theta.clear();Hadron_phi.clear();
		 Hadron_z.clear(); Hadron_pt.clear();Hadron_phiTrento.clear();
		Particles_LUND.clear();Hadron_parent.clear();Hadron_pid.clear();
		Hadron_xf.clear();Hadron_eta.clear();
		Hadron_vz.clear();

		CountsID=0;
		// TODO Auto-generated method stub
		LorentzVector Ve = new LorentzVector ( 0,0,0,0); 
		// HipoGroup MCParticle = event.getGroup("MC::Particle");
		// Read the electrons and photons 
	//	 bankMC.show();
		//MCParticle.show();
		int nrows = MCParticle.getRows();
		float px=-100,py=-100,pz=-100,e_vx=-100,e_vy=-100,e_vz=-100;
		double e_mom=0, e_theta=0, e_phi=0;
		double Q2=0, xB=0;
		Vector3 ParticleVector= new Vector3(0,0,0);
		//MCParticle.show();
		//bankMC.show();
		LorentzVector VGS = new LorentzVector(0,0,0,0);
		LorentzVector vecW2= new LorentzVector(0,0,0,0);
		double Pq=0;
		double Pl=0;
		int nrows2 = bankMC.getRows();
		//MCParticle.show();
	//	System.out.println(" EVENTO MC ---------------------------------- ");
		for (int ipart =0 ; ipart<nrows; ipart++) {
			
			//System.out.println(MCParticle.getInt("pid",ipart));
			if(11 == MCParticle.getInt("pid",ipart)) {
				//MCParticle.show();
				px = MCParticle.getFloat("px",ipart);
				py = MCParticle.getFloat("py",ipart);
				pz = MCParticle.getFloat("pz",ipart);
				e_vx = MCParticle.getFloat("vx",ipart);
				e_vy = MCParticle.getFloat("vy",ipart);
				e_vz = MCParticle.getFloat("vz",ipart);
				e_mom = Math.sqrt(px*px+py*py+pz*pz);
				Ve.setPxPyPzE(px,py,pz,e_mom);
				
				this.e_p=e_mom;
	//			System.out.println(" Elettrone px " + px);
				//THIS IS WHAT MAKES IT CRASH: ParticleVector	
			   // ParticleVector.setXYZ(px,py,pz);
			    double myphi = Math.toDegrees(Math.atan(py/px));
			    double mytheta = Math.toDegrees(Math.acos(pz/e_mom));
			   
			    if (px<0 && py>0) myphi+=180;
			    else if (px<0 && py<0) myphi-=180;
			    this.e_phi=myphi;
			    this.e_theta=mytheta;
				//e_theta= Math.toDegrees(ParticleVector.theta());
				e_phi= Math.toDegrees(ParticleVector.phi());
			   // System.out.println("Theta from ParticleVector " + e_theta + " Phi " + e_phi + " Mine calculation is theta " + mytheta+" Phi "+myphi);
			
				break;
			}
		}
	
		VGS.add(this.VB);
		VGS.sub(Ve);
		Q2=-VGS.mass2(); 
		this.e_Q2=Q2;
		Pq = -this.VT.px()*VGS.px()-this.VT.py()*VGS.py()-this.VT.pz()*VGS.pz()+this.VT.e()*VGS.e();
		vecW2.add(this.VT);vecW2.add(this.VB);vecW2.sub(Ve);
		Pl= -this.VT.px()*this.VB.px()-this.VT.py()*this.VB.py()-this.VT.pz()*this.VB.pz()+this.VT.e()*this.VB.e();
		xB=Q2/(2*Pq);
		this.e_xB=xB;
		this.y_el_MC = Pq/Pl;
		XB_el_MC = xB; 
		Q2_el_MC=-VGS.mass2();
		W_el_MC= vecW2.mass();
		this.e_W=W_el_MC;
		//HipoGroup bankMC = event.getGroup("MC::Lund");
		//bankMC.show();
	

		double p_px=-10,p_py=-10,p_pz=-10, p_E=0;
		  
		
		    // Modify to look over the whole LUND structure
		    
		   // for(int k = 0; k < nrows; k++){
		for(int k = 0; k < nrows2; k++){
				//MCParticle.show();
				//bankMC.show();
				//System.out.println(Q2);
				//Hadron_pid.add(MCParticle.getInt("pid",k));
				//p_px = MCParticle.getFloat("px",k);
				//p_py = MCParticle.getFloat("py",k);
				//p_pz = MCParticle.getFloat("pz",k);
				//p_E = MCParticle.getFloat("energy",k);
			//	double mass=0.15;
				//p_E = (float) Math.sqrt(Math.pow(p_px*p_px+p_py*p_py+p_pz*p_pz,2)+Math.pow(mass,2));
		//		int parentID= 0;
				//System.out.println(" Ciao sono qui " + bankMC.getInt("pid",k));
			
			Hadron_pid.add(bankMC.getInt("pid",k));
			Hadron_vz.add((double)bankMC.getFloat("vz",k));
				p_px = (double) bankMC.getFloat("px",k);
				p_py =(double) bankMC.getFloat("py",k);
				p_pz =(double) bankMC.getFloat("pz",k);
				p_E =(double) bankMC.getFloat("energy",k);
				int parentID= bankMC.getInt("parent", k);
				//Considering only string produced pions
				//if((bankMC.getInt("pid",k)==113 || bankMC.getInt("pid",k)==213 || bankMC.getInt("pid",k)==-213) && Math.abs(bankMC.getInt("pid",parentID-1))<100 ) {
				//if((bankMC.getInt("pid",k)==223 ) && Math.abs(bankMC.getInt("pid",parentID-1))<100 ) {

			//if(bankMC.getInt("pid",k)==Particle && Math.abs(bankMC.getInt("pid",parentID-1))>100 &&Math.abs(bankMC.getInt("pid",parentID-1))<400 ) {
			//	if(bankMC.getInt("pid",k)==Particle && Math.abs(bankMC.getInt("pid",parentID-1))<100 ) {

					//if(bankMC.getInt("pid",k)==Particle && Math.abs(bankMC.getInt("pid",parentID-1))<100  ) {
			//if(bankMC.getInt("pid",k)==Particle ) {
				
				//	if(bankMC.getInt("pid",k)==Particle&& (Math.abs(bankMC.getInt("pid",parentID-1))==113 ||Math.abs(bankMC.getInt("pid",parentID-1))==213 || Math.abs(bankMC.getInt("pid",parentID-1))==-213 )) {
			int parentofparent = bankMC.getInt("parent",parentID-1);
	
				//if(bankMC.getInt("pid",k)==Particle&& bankMC.getInt("pid",parentID-1)==223 && bankMC.getInt("pid",parentofparent-1)<100  ) {

				 //CountsParentID.fill(bankMC.getInt("pid", parentID-1));
				//System.out.println(  Math.abs(bankMC.getInt("pid",parentID-1)));
				this.CountsID++;
			//	System.out.println(" The PARENT IS " + bankMC.getInt("pid",parentID-1));
			//	bankMC.show();
			    
			    LorentzVector PHM = new LorentzVector(0,0,0,0);
			    LorentzVector Gammacm = new LorentzVector(0,0,0,0);
			    LorentzVector Protoncm = new LorentzVector(0,0,0,0);
			    LorentzVector PIcm = new LorentzVector(0,0,0,0);
			    LorentzVector Vphotons = new LorentzVector(0,0,0,0);
			
				PHM.setPxPyPzE(p_px,p_py,p_pz,p_E);
				LorentzVector X = new LorentzVector(0,0,0,0);X.add(VT).add(VB).sub(Ve).sub(PHM);
                MX=X.mass();
        //        System.out.println(" Energy beam "+ VT.e() + " energy beam "+VB.e()+ " energy erlectron "+ Ve.e() +" Pion "+ PHM.e());
				//LorentzVector PHM = new LorentzVector(p_px,p_py,p_pz,p_E);
				double PPh = -this.VT.px()*PHM.px()-this.VT.py()*PHM.py()-this.VT.pz()*PHM.pz()+this.VT.e()*PHM.e();
				Vector3 VectorPh = VGS.vect(); // Virtual photon is beam-el scattered but adding the target doesent change tri momenta
				Vector3 VectorPi=PHM.vect();
				Vector3 VectorP = this.VT.vect();
				Vector3 rotationVector = new Vector3(0,0,0);
				TransMatrix rotationMatrix = new TransMatrix();
				rotationMatrix.compose(VectorPh);
				Vector3 PhCM = rotationMatrix.mult(VectorPh);
				Vector3 PiCM = rotationMatrix.mult(VectorPi);
				Vector3 PCM = rotationMatrix.mult(VectorP);
				Gammacm.setPxPyPzE(PhCM.x(), PhCM.y(), PhCM.z(), VGS.e());
				Protoncm.setPxPyPzE(PCM.x(), PCM.y(), PCM.z(), this.VT.e());
				PIcm.setPxPyPzE(PiCM.x(), PiCM.y(), PiCM.z(), PHM.e());
				Vphotons.add(Gammacm).add(Protoncm);
				Vector3 boost = Vphotons.boostVector();
				boost.negative();
				Vphotons.boost(boost);
				Gammacm.boost(boost);
				Protoncm.boost(boost);
				PIcm.boost(boost);
				double Pi0M= PHM.mass();
				double Pls = PIcm.vect().dot(Gammacm.vect())/Gammacm.vect().mag();
				double PT=PIcm.vect().mag2()-Pls*Pls;  
			//double PT=Math.sqrt(PIcm.vect().mag2()-Pls*Pls);

				double Z_var= PPh/Pq;
	           double Xf = (2 * Pls) / vecW2.mass();
				Vector3 VectorB = VB.vect();
				Vector3 VectorPhNorm = new Vector3(VectorPh.x()/VectorPh.mag(),VectorPh.y()/VectorPh.mag(),VectorPh.z()/VectorPh.mag());				
				Vector3 Angl1 = VectorPhNorm.cross(VectorB);
				Vector3 Angl2 = VectorPhNorm.cross(VectorPi);
				Angl1.setXYZ(Angl1.x()/Angl1.mag(), Angl1.y()/Angl1.mag(),Angl1.z()/Angl1.mag());
				Angl2.setXYZ(Angl2.x()/Angl2.mag(), Angl2.y()/Angl2.mag(),Angl2.z()/Angl2.mag());
				double AngoloCos=Angl1.dot(Angl2);
				double Angolo= Math.toDegrees(Math.acos(AngoloCos));
				//this.Phi_LAB = Math.toDegrees(Math.acos(Angolo/(Angl1.mag()*Angl2.mag())));
				Vector3 AngleA = VectorB.cross(VectorPi);
				double AngleAA = AngleA.dot(VectorPhNorm);
				Vector3 AngleB=VectorPhNorm.cross(VectorB);
				Vector3 AngleBB =VectorPhNorm.cross(VectorPi);
				double AngleSin= AngleAA/(AngleB.mag()*AngleBB.mag());
				double AlternativePhi= Math.toDegrees(Math.asin(AngleSin));
				//System.out.println("----------");
				//System.out.println(" Phi cos " +Angl1.dot(Angl2) + " sin " + AngleSin );
				double Phi_LAB = Angolo;
				//if(AngleSin<0) this.Phi_LAB=-this.Phi_LAB; 
				//System.out.println("Angolo LAB [degrees] " + Angolo);
				//if(AlternativePhi<0) AlternativePhi=-AlternativePhi; 
				//System.out.println(" Angolo LAB Alternativo " + AlternativePhi);
				if(Angolo>0 && AlternativePhi>0) Phi_LAB=Angolo;
				else if(Angolo<0 && AlternativePhi>0) Phi_LAB=Angolo;
				else if (Angolo<0  && AlternativePhi<0) Phi_LAB=Angolo+90;
				else if (Angolo>0 && AlternativePhi<0 ) Phi_LAB=360-Angolo;
				//System.out.println("PID " + Particle +  " My generted  Phi  is " + Phi_LAB + "  Particle Mass " + Pi0M + " The electron  momentum is (x,y,z):  " + px  +" " + py +" "+pz );
				
				// Cose di Harut.
				 double ehs, pcm2;
				 double pts;
				 double mpro2 = VT.e()*VT.e();
				 double mhad2 = PHM.mass2();
				 double anu = VB.e() - VB.e();
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
				
				//System.out.println( " XF HARUT "+xfhad + " Mine "+ Xf );
				
				
				
				double momhadron=  Math.sqrt(p_px*p_px+p_py*p_py+p_pz*p_pz);
				double thethadron =  Math.toDegrees(PHM.vect().theta());
				double phihadron =  Math.toDegrees(PHM.vect().phi());
		    	
				
				//if(MX>Parameters.Get_Hadron_MissingMassCut()) {
					//if(Z_var>0.6 && xB>0.4 && Q2>4 && PT>0.4)System.out.println(" Generated MX " + MX);

		    	Hadron_en.add(p_E);
		    	Hadron_mom.add( momhadron);
		    	Hadron_theta.add(thethadron);
		        Hadron_phi.add(phihadron);
		        Hadron_phiTrento.add( Phi_LAB);
		        Hadron_mx.add(MX);
		        
		    	Hadron_z.add(Z_var);
		    	Hadron_pt.add( PT);
		    	Hadron_parent.add(parentID);
		    	Hadron_xf.add(Xf);
		    	Hadron_eta.add(etabr);
				
		    	if(bankMC.getInt("pid",k)==Particle ) {
	
				TuplHadron_momentum.add((double)momhadron);
				TuplHadron_theta.add((double)thethadron);
				TuplHadron_phiclas.add((double)phihadron);
				TuplHadron_z.add(Z_var);
				TuplHadron_pt.add(PT);
				TuplHadron_phi.add(Phi_LAB) ;
				TuplHadron_xf.add(Xf);
				TuplHadron_MX.add(MX);
				TuplHadron_eta.add(etabr);
				TuplHadron_father.add(bankMC.getInt("pid",parentID-1));
				TuplHadron_granfather.add(bankMC.getInt("pid",parentofparent-1));
				
		//}// MX CUT
			}
			
		}
		
		Bank tuple = new Bank(schemaMC,TuplHadron_momentum.size());
		for (int kk=0; kk<TuplHadron_momentum.size();kk++) {
			tuple.putLong("event", kk, evento);
		tuple.putInt("pid", kk, Particle);
		tuple.putDouble("e_p", kk,   this.e_p);
		tuple.putDouble("e_theta", kk,  this.e_theta);
		tuple.putDouble("e_phi", kk,  this.e_phi);
		tuple.putDouble("Q2", kk,  this.e_Q2);
		tuple.putDouble("W", kk,  this.e_W);
		tuple.putDouble("x", kk,  this.e_xB);
		tuple.putDouble("y", kk,  this.y_el_MC);
		tuple.putDouble("pi_p", kk,  TuplHadron_momentum.get(kk));
		tuple.putDouble("pi_theta", kk, TuplHadron_theta.get(kk));
		tuple.putDouble("pi_phi", kk, TuplHadron_phiclas.get(kk) );
		tuple.putDouble("z", kk,  TuplHadron_z.get(kk));
		tuple.putDouble("pt", kk, TuplHadron_pt.get(kk));
		tuple.putDouble("phi", kk, TuplHadron_phi.get(kk));
		tuple.putDouble("xf", kk, TuplHadron_xf.get(kk));
		tuple.putDouble("mX", kk, TuplHadron_MX.get(kk));
		tuple.putDouble("eta", kk, TuplHadron_eta.get(kk) );
		tuple.putInt("father", kk, TuplHadron_father.get(kk) );
		tuple.putInt("granfather", kk, TuplHadron_granfather.get(kk) );
	//	System.out.println(" id "+ kk+ " mom " + TuplHadron_momentum.get(kk) +" e_theta "+ this.e_theta+" e_phi " +this.e_phi+ " h_theta " +TuplHadron_theta.get(kk)+ " h_phi "+TuplHadron_phiclas.get(kk) + " z " + TuplHadron_z.get(kk)+ " pt " + TuplHadron_pt.get(kk));
		}		
		Event  eventoMC = new Event();
		eventoMC.reset();
		eventoMC.write(tuple);
	    writeroutputMC.addEvent(eventoMC);	
	}

	
	
	
	/**
	 * Process another generated event 
	 * @param event
	 */
	public void addGenerated(Bank bankMC, int Particle ) {
		Hadron_en.clear();Hadron_mom.clear(); Hadron_theta.clear();Hadron_phi.clear();
		 Hadron_z.clear(); Hadron_pt.clear();Hadron_phiTrento.clear();
		Particles_LUND.clear();Hadron_parent.clear();Hadron_pid.clear();
		Hadron_xf.clear();	Hadron_eta.clear();
		CountsID=0;
		// TODO Auto-generated method stub
		LorentzVector Ve = new LorentzVector ( 0,0,0,0); 

		int nrows = bankMC.getRows();
		float px=-100,py=-100,pz=-100,e_vx=-100,e_vy=-100,e_vz=-100;
		double e_mom=0, e_theta=0, e_phi=0;
		double Q2=0, xB=0;
		Vector3 ParticleVector= new Vector3(0,0,0);
		LorentzVector VGS = new LorentzVector(0,0,0,0);
		LorentzVector vecW2= new LorentzVector(0,0,0,0);
		double Pq=0;
		double Pl=0;
		int nrows2 = bankMC.getRows();
		for (int ipart =0 ; ipart<nrows; ipart++) {	
			if(11 == bankMC.getInt("pid",ipart) && bankMC.getInt("parent", ipart)==1) {
				//bankMC.show();
				px = bankMC.getFloat("px",ipart);
				py = bankMC.getFloat("py",ipart);
				pz = bankMC.getFloat("pz",ipart);
				e_vx = bankMC.getFloat("vx",ipart);
				e_vy = bankMC.getFloat("vy",ipart);
				e_vz = bankMC.getFloat("vz",ipart);
				e_mom = Math.sqrt(px*px+py*py+pz*pz);
				Ve.setPxPyPzE(px,py,pz,e_mom);
				H_e_vx_vy.fill(e_vx, e_vy);;
			
				this.e_p=e_mom;
			    double myphi = Math.toDegrees(Math.atan(py/px));
			    double mytheta = Math.toDegrees(Math.acos(pz/e_mom));
			    if (px<0 && py>0) myphi+=180;
			    else if (px<0 && py<0) myphi-=180;
			    this.e_phi=myphi;
			    this.e_theta=mytheta;
				e_phi= Math.toDegrees(ParticleVector.phi());
			
				break;
			}
		}
		VGS.add(this.VB);
		VGS.sub(Ve);
		Q2=-VGS.mass2(); 
		this.e_Q2=Q2;
		Pq = -this.VT.px()*VGS.px()-this.VT.py()*VGS.py()-this.VT.pz()*VGS.pz()+this.VT.e()*VGS.e();
		vecW2.add(this.VT);vecW2.add(this.VB);vecW2.sub(Ve);
		Pl= -this.VT.px()*this.VB.px()-this.VT.py()*this.VB.py()-this.VT.pz()*this.VB.pz()+this.VT.e()*this.VB.e();
		xB=Q2/(2*Pq);
		this.e_xB=xB;
		y_el_MC = Pq/Pl;
		XB_el_MC = xB; 
		Q2_el_MC=-VGS.mass2();
		W_el_MC= vecW2.mass();
		this.e_W=W_el_MC;
//		bankMC.show();
		double p_px=-10,p_py=-10,p_pz=-10, p_E=0;
		for(int k = 0; k < nrows2; k++){
				Hadron_pid.add(bankMC.getInt("pid",k));
				p_px =(double) bankMC.getFloat("px",k);
				p_py =(double) bankMC.getFloat("py",k);
				p_pz =(double) bankMC.getFloat("pz",k);
				p_E =(double) bankMC.getFloat("E",k);
				int parentID= bankMC.getInt("parent", k);
				//Considering only string produced pions
				//if((bankMC.getInt("pid",k)==113 || bankMC.getInt("pid",k)==213 || bankMC.getInt("pid",k)==-213) && Math.abs(bankMC.getInt("pid",parentID-1))<100 ) {
				if((bankMC.getInt("pid",k)==223 ) && Math.abs(bankMC.getInt("pid",parentID-1))<100 ) {

			//if(bankMC.getInt("pid",k)==Particle && Math.abs(bankMC.getInt("pid",parentID-1))>100 &&Math.abs(bankMC.getInt("pid",parentID-1))<400 ) {
			//	if(bankMC.getInt("pid",k)==Particle && Math.abs(bankMC.getInt("pid",parentID-1))<100 ) {

				//	if(bankMC.getInt("pid",k)==Particle && Math.abs(bankMC.getInt("pid",parentID-1))<100  ) {
			//if(bankMC.getInt("pid",k)==Particle ) {
					//if(bankMC.getInt("pid",k)==Particle&& (Math.abs(bankMC.getInt("pid",parentID-1))==113 ||Math.abs(bankMC.getInt("pid",parentID-1))==213 || Math.abs(bankMC.getInt("pid",parentID-1))==-213 )) {
					int parentofparent = bankMC.getInt("parent",parentID-1);
	
				//if(bankMC.getInt("pid",k)==Particle&& bankMC.getInt("pid",parentID-1)==223 && bankMC.getInt("pid",parentofparent-1)<100  ) {
				 //CountsParentID.fill(bankMC.getInt("pid", parentID-1));
				//System.out.println(  Math.abs(bankMC.getInt("pid",parentID-1)));
				this.CountsID++;

			    LorentzVector PHM = new LorentzVector(0,0,0,0);
			    LorentzVector Gammacm = new LorentzVector(0,0,0,0);
			    LorentzVector Protoncm = new LorentzVector(0,0,0,0);
			    LorentzVector PIcm = new LorentzVector(0,0,0,0);
			    LorentzVector Vphotons = new LorentzVector(0,0,0,0);		
				PHM.setPxPyPzE(p_px,p_py,p_pz,p_E);
				LorentzVector X = new LorentzVector(0,0,0,0);X.add(VT).add(VB).sub(Ve).sub(PHM);
                MX=X.mass();
				double PPh = -this.VT.px()*PHM.px()-this.VT.py()*PHM.py()-this.VT.pz()*PHM.pz()+this.VT.e()*PHM.e();
				Vector3 VectorPh = VGS.vect(); // Virtual photon is beam-el scattered but adding the target doesent change tri momenta
				Vector3 VectorPi=PHM.vect();
				Vector3 VectorP = this.VT.vect();
				Vector3 rotationVector = new Vector3(0,0,0);
				TransMatrix rotationMatrix = new TransMatrix();
				rotationMatrix.compose(VectorPh);
				Vector3 PhCM = rotationMatrix.mult(VectorPh);
				Vector3 PiCM = rotationMatrix.mult(VectorPi);
				Vector3 PCM = rotationMatrix.mult(VectorP);
				Gammacm.setPxPyPzE(PhCM.x(), PhCM.y(), PhCM.z(), VGS.e());
				Protoncm.setPxPyPzE(PCM.x(), PCM.y(), PCM.z(), this.VT.e());
				PIcm.setPxPyPzE(PiCM.x(), PiCM.y(), PiCM.z(), PHM.e());
				Vphotons.add(Gammacm).add(Protoncm);
				Vector3 boost = Vphotons.boostVector();
				boost.negative();
				Vphotons.boost(boost);
				Gammacm.boost(boost);
				Protoncm.boost(boost);
				PIcm.boost(boost);
				double Pi0M= PHM.mass();
				double Pls = PIcm.vect().dot(Gammacm.vect())/Gammacm.vect().mag();
				double PT=PIcm.vect().mag2()-Pls*Pls;  
				double Z_var= PPh/Pq;
	            double Xf = (2 * Pls) / vecW2.mass();
				Vector3 VectorB = VB.vect();
				Vector3 VectorPhNorm = new Vector3(VectorPh.x()/VectorPh.mag(),VectorPh.y()/VectorPh.mag(),VectorPh.z()/VectorPh.mag());				
				Vector3 Angl1 = VectorPhNorm.cross(VectorB);
				Vector3 Angl2 = VectorPhNorm.cross(VectorPi);
				Angl1.setXYZ(Angl1.x()/Angl1.mag(), Angl1.y()/Angl1.mag(),Angl1.z()/Angl1.mag());
				Angl2.setXYZ(Angl2.x()/Angl2.mag(), Angl2.y()/Angl2.mag(),Angl2.z()/Angl2.mag());
				double AngoloCos=Angl1.dot(Angl2);
				double Angolo= Math.toDegrees(Math.acos(AngoloCos));
				Vector3 AngleA = VectorB.cross(VectorPi);
				double AngleAA = AngleA.dot(VectorPhNorm);
				Vector3 AngleB=VectorPhNorm.cross(VectorB);
				Vector3 AngleBB =VectorPhNorm.cross(VectorPi);
				double AngleSin= AngleAA/(AngleB.mag()*AngleBB.mag());
				double AlternativePhi= Math.toDegrees(Math.asin(AngleSin));
				double Phi_LAB = Angolo;
				if(Angolo>0 && AlternativePhi>0) Phi_LAB=Angolo;
				else if(Angolo<0 && AlternativePhi>0) Phi_LAB=Angolo;
				else if (Angolo<0  && AlternativePhi<0) Phi_LAB=Angolo+90;
				else if (Angolo>0 && AlternativePhi<0 ) Phi_LAB=360-Angolo;
		    	this.InvMass.fill(Z_var, PT);	
		    	Hadron_en.add(p_E);
		    	Hadron_mom.add( Math.sqrt(p_px*p_px+p_py*p_py+p_pz*p_pz));
		    	Hadron_theta.add( Math.toDegrees(PHM.vect().theta()));
		        Hadron_phi.add( Math.toDegrees(PHM.vect().phi()));
		        Hadron_phiTrento.add( Phi_LAB);
		    	Hadron_z.add( Z_var);
		    	Hadron_pt.add(PT);
		    	Hadron_parent.add(parentID);
		    	Hadron_xf.add(Xf);
				//System.out.println("Bins " + Hadron_Lund.FstBin + " | "+ Hadron_Lund.SndBin +" |  " + Hadron_Lund.TrdBin + "|  "+ Hadron_Lund.FrtBin);
				//System.out.println("First Min "  + Hadron_Lund.Fst_Min + " First max "+ Hadron_Lund.Fst_Max + " Second min "+ Hadron_Lund.Snd_Min + " Second Max" + Hadron_Lund.Snd_Max);
				//System.out.println(" Third Min " + Hadron_Lund.Trd_Min + " Third_Max " + Hadron_Lund.Trd_Max + " Forth Min " + Hadron_Lund.Frt_Min + " Forth Max " + Hadron_Lund.Frt_Max);
				int BinXB,BinQ2;
		  

		    	
		    	
		    	
			}
			
		}
		
	}
	/**
	 *  Set the Parent ID for searching PID
	 * @param i
	 */
	public void setParent(int i) {
		// Check if i has been already used in the class
		this.parentID=i;
	}

	/**
	 * Count how many particles in the events
	 * @param String of the bin : z or pt
	 * @return the nr of particles with the given id in the events analyzided
	 * @throws FileNotFoundException 
	 * @throws IOException 
	 */
	public void getNr(String s) throws FileNotFoundException {
		File File_txt = new File("/Users/gangelini/Desktop/Plots/Lund.txt");
		if(File_txt.exists()) {
			System.out.println("File found" );
			@SuppressWarnings("resource")
			PrintWriter outP = new PrintWriter(File_txt);
			for(int ii=0;ii<Z_Bins; ii++) {
				for (int jj=0;jj<PT_Bins; jj++) {
					System.out.println(" BIN Z: " +ii+ " BIN PT: " +jj);
					System.out.println(InvMass.getBinContent(ii, jj));
					outP.println(" BIN Z: " +ii+ " BIN PT: " +jj+" Counts Pi0s" + InvMass.getBinContent(ii, jj) );
				}
			} 
		}
		System.out.print(" Z COUNTS" );

		// TODO Auto-generated method stub
		if (s=="z") {
			System.out.println(" --- Z ---- " );
			System.out.println("MC counts Z<=0.1 " + this.z1);
			System.out.println("MC count 0.1<Z bin <=0.2 " + this.z2);
			System.out.println("MC count 0.2<Z bin <=0.3 " + this.z3);
			System.out.println("MC count 0.3<Z bin <=0.4 " + this.z4);
			System.out.println("MC count 0.4<Z bin <=0.5 "+ this.z5);
			System.out.println("MC count 0.5<Z bin <=0.6 " + this.z6);
			System.out.println("MC count 0.6<Z bin <=0.7 "+ this.z7);
			System.out.println("MC count 0.7<Z bin <=0.8 "+this.z8);
			System.out.println("MC count 0.8<Z bin <=0.9 "+ this.z9);
			System.out.println("MC count 0.9<Z bin <=1 " + this.z10);

		}
		else if ( s=="pt") {
			System.out.println(" --- PT ---- " );
			System.out.println("MC counts PT<=0.2 " + this.pt1);
			System.out.println("MC count 0.2<PT <=0.3 " + this.pt2);
			System.out.println("MC count 0.3<PT <=0.4 " + this.pt3);
			System.out.println("MC count 0.4<PT  <=0.6 " + this.pt4);
			System.out.println("MC count 0.6<PT <=0.8 "+ this.pt5);
			System.out.println("MC count 0.8<PT  <=1 " + this.pt6);
			System.out.println("MC count 1<PT  <=1.5 "+ this.pt7);

		}

	}

	public void setTarget(LorentzVector target) {
		// TODO Auto-generated method stub
		this.VT=target;

	}

	public void setBeam(LorentzVector beam) {
		// TODO Auto-generated method stub
		this.VB=beam;
	}

	/**
	 * 
	 * @return the MC electron Q2 for the event 
	 */
	public double getQ2()
	{
		return this.Q2_el_MC;
	}

	/**
	 * 
	 * @return the MC electron y for the event 
	 */
	public double getY()
	{
		return this.y_el_MC;
	}

	/**	 
	 * @return the MC electron Q2 for the event 
	 */
	public double getW()
	{
		return this.W_el_MC;
	}		

	public double getXB() {
		return this.XB_el_MC;
	}

	public void Histo(String time, String workout) {
		
		new File(workout+"/Plots").mkdir();
		String PathFolder = workout+"/Plots/"+time;
		File dir = new File(PathFolder);
		dir.mkdir();
		PathFolder=workout+"/Plots/"+time+"/LUND";
		new File (PathFolder).mkdir();


		System.out.println(" (Pi) - > I am Plotting LUND information");

		

		//I could save the invairiant mass somewhere.


	}

public int CountParticles()
{
	return CountsID;
}


/**
 * Save Histograms 
 */
public void SaveHistograms_MultiBins(String workdirout,String time) {
	 //CountsPion.draw(CountsParentID);
	// Producing the H1F for the Getting counts as function of Phi 
    System.out.println(" I am saving FinalHipo"); 
	writeroutputMC.close();
}

public  ArrayList<Double> Get_H_En() {
	return Hadron_en;
}

public  ArrayList<Double> Get_H_mom() {
	return Hadron_mom;
}

public  ArrayList<Double> Get_H_theta() {
	return Hadron_theta;
}
public  ArrayList<Double> Get_H_phi() {
	return Hadron_phi;
}

public ArrayList<Double> Get_H_phiTrento(){
	return Hadron_phiTrento;
}
public  ArrayList<Double> Get_H_z() {
	return Hadron_z;
}
public  ArrayList<Double> Get_H_PT() {
	return Hadron_pt;
}

public ArrayList<Integer> Get_H_PID(){
	return Hadron_pid;
}

public ArrayList<Integer> Get_H_Parent(){
	return Hadron_parent;
}

public ArrayList<Double> Get_H_XF(){
	return Hadron_xf;
}

public ArrayList<Double> Get_H_etaB(){
	return Hadron_eta;
}

public ArrayList<Double> Get_H_MX(){
	return Hadron_mx;
}
public ArrayList<Double> Get_H_vz(){
	return Hadron_vz;
}
 


public void ResetHadrons () {
	
	Hadron_pt.clear();Hadron_z.clear();Hadron_phi.clear();Hadron_theta.clear();Hadron_mom.clear();Hadron_en.clear();
	Hadron_phiTrento.clear();Hadron_parent.clear();Hadron_pid.clear();Hadron_xf.clear();Hadron_eta.clear();Hadron_mx.clear();
	Hadron_vz.clear();
	
}




}
