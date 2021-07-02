package PHD_Clean;

import org.jlab.groot.data.H1F;
import org.jlab.groot.data.TDirectory;
import org.jlab.groot.ui.TGCanvas;
import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.data.SchemaFactory;
import org.jlab.jnp.hipo4.io.HipoReader;

public class CountParticles {

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		String PathData = "/Users/gangelini/work/Comparisons/Data_Standard/";
		String PathMC = "/Users/gangelini/work/Comparisons/MC_Standard/";
		TGCanvas Q2C = new TGCanvas("Q2C","Q2C",1600,1400);
		Q2C.divide(2,2);
		H1F q2d = new H1F("q2d",300,1,3);
		H1F q2m = new H1F("q2m",300,1,3);
		H1F etheta  = new H1F("etheta",900,-180,180);
		H1F ep  = new H1F("epx",300,0,5);
		HipoReader readerData = new HipoReader(); readerData.open(PathData+ "Skimmed_Hipo.hipo");
		Event eventData=new Event(); 
		SchemaFactory schemaData = readerData.getSchemaFactory();
		Bank DSTData = new Bank(schemaData.getSchema("SIDIS::ePiX"));
		int CountData=0;
		while(readerData.hasNext()==true){
			if(CountData%100000 == 0) System.out.println(CountData + " Data events");

			readerData.nextEvent(eventData);
			eventData.read(DSTData); 
			double Q2 =DSTData.getFloat("Q2", 0);
			double W = DSTData.getFloat("W", 0);
			double x = DSTData.getFloat("x", 0);
			double y = DSTData.getFloat("y", 0);
			double mX=DSTData.getFloat("mX", 0);
			double pt=DSTData.getFloat("pt", 0);
			double z= DSTData.getFloat("z", 0);
			CountData++;
			q2d.fill(Q2);
		}		
		
		int CountMC=0;
		Event eventMC=new Event(); 

		HipoReader readerMC = new HipoReader(); readerMC.open(PathMC+ "Skimmed_Hipo.hipo");
		SchemaFactory schemaMC = readerMC.getSchemaFactory();
		Bank DSTMC = new Bank(schemaMC.getSchema("SIDIS::ePiX"));
		while(readerMC.hasNext()==true){
			//if(CountMC%100000 == 0) System.out.println(CountMC + " MC events");

			readerMC.nextEvent(eventMC);
			eventMC.read(DSTMC); 
			double Q2 =DSTMC.getFloat("Q2", 0);
			double W = DSTMC.getFloat("W", 0);
			double x = DSTMC.getFloat("x", 0);
			double y = DSTMC.getFloat("y", 0);
			double mX=DSTMC.getFloat("mX", 0);
			double pt=DSTMC.getFloat("pt", 0);
			double z= DSTMC.getFloat("z", 0);
			double e_p= DSTMC.getFloat("e_p", 0);
			double e_theta= DSTMC.getFloat("e_theta", 0);
			double e_phi= DSTMC.getFloat("e_phi", 0);
			double pi_p= DSTMC.getFloat("pi_p", 0);
			double pi_theta= DSTMC.getFloat("pi_theta", 0);
			double pi_phi= DSTMC.getFloat("pi_phi", 0);

			
			
			CountMC++;
			q2m.fill(Q2);
			if(Q2<2.05 && Q2>2.01) {
				if(Math.toDegrees(e_theta)>9 &&Math.toDegrees(e_theta)<10) { 
				//System.out.println(Math.toDegrees(e_theta));
				if(Math.toDegrees(e_theta)> 9.823 && Math.toDegrees(e_theta)<9.827  )
				{etheta.fill(Math.toDegrees(e_phi));
			//	System.out.println(" theta "+e_theta + " phi" +  e_phi+ " p "+e_p );
					if(Math.toDegrees(e_phi)> 127.6  && Math.toDegrees(e_phi)< 128 ) {
					System.out.println(" pion p " +pi_p + "theta "+pi_theta + " phi "+pi_phi  );
					}
				}
				}
		
			ep.fill(e_p);
			}
		}		
		
	//	System.out.println("Entries MC " +CountMC + "Entries Data" +CountData );
		Q2C.cd(0);Q2C.draw(q2m);
		Q2C.cd(1);Q2C.draw(q2d);
		Q2C.cd(2);Q2C.draw(etheta);
		Q2C.cd(3);Q2C.draw(ep);
	}

}
