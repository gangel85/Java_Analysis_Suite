package PHD_Clean;
import java.awt.geom.Path2D;
import java.awt.geom.Point2D;
import java.io.File;
import java.util.ArrayList;
import java.util.List;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.data.H3F;
import org.jlab.groot.graphics.EmbeddedCanvas;

public class MultiBins {

	/**
	 * The Constructor is used to create  multi dimensional bins up to 5 dimensions. 
	 * If all the parameters are set not to zero then the class creates a double array of H3F.
	 * If Parameter 5 is set to 0 the class creates a double array of H2F
	 * If Parameter 4 and 5 are set to 0 the class creates a double array of H1F
	 * -- To be Implemented : Different bin Size by using Custom Axis in the Histograms
	 * @param FirstBin Number of bins for the first variable
	 * @param SecondBin Number of bins for the second variable
	 * @param ThirdBin Number of bins for the third variable
	 * @param ForthBin Number of bins for the fourth variable
	 * @param FifhtBin Number of bins for the fifth variable
	 */
	public int Dimension=0;
	public int FstBin=0,SndBin=0,TrdBin=0,FrtBin=0, FthBin=0;
	public double  Fst_Min=0,Fst_Max=0 , Snd_Min=0, Snd_Max=0, Trd_Min=0,Trd_Max=0,Frt_Min=0, Frt_Max=0, Fth_Min=0, Fth_Max=0;
	
	private double[] bin1,bin2,bin3,bin4,bin5;
	
	private boolean even=true;
	
	
	
	public List<List<H3F>> MultiHisto3D = new ArrayList<List<H3F>>();
	public List<List<H2F>> MultiHisto2D = new ArrayList<List<H2F>>();
	public List<List<H1F>> MultiHisto1D = new ArrayList<List<H1F>>();
	private String name;
	public H2F Counts2D;
	
	private String Name1st=null;
	private String Name2nd=null;
	private String Name3rd=null;
	private String Name4th=null;
	private String Name5th=null;
	
	private int PolygonNr = 0;
	private List<Path2D> Polygons = new ArrayList<Path2D>(); 
	
	private String HistoName=null;

	private boolean PolygonalBinClas=false;
	
	public  MultiBins(int FirstBins, int SecondBins, int ThirdBins, int ForthBins , int FifhtBins) {
		
		this.FstBin=FirstBins; this.SndBin= SecondBins;
		if (FifhtBins !=0 &&ForthBins !=0  && ThirdBins!=0 && SecondBins!=0 &&FirstBins!=0) {
			this.FrtBin= ForthBins; this.FthBin=FifhtBins;	this.TrdBin=ThirdBins;
		
			this.Dimension=5;
		}
		else if (FifhtBins ==0 && ForthBins !=0  && ThirdBins!=0 && SecondBins!=0 &&FirstBins!=0 ) {
			this.FrtBin= ForthBins;	this.TrdBin=ThirdBins;
			this.Dimension=4;
		}
		else if (FifhtBins ==0 && ForthBins ==0  && ThirdBins!=0 && SecondBins!=0 &&FirstBins!=0 ) {
			this.TrdBin=ThirdBins;
			this.Dimension=3;
		}
		
		else if (FifhtBins ==0 && ForthBins ==0  && ThirdBins==0 && SecondBins!=0 &&FirstBins!=0  ) {
			this.Dimension=2;
			this.InizializeClasPolygons();
			
		}
		else System.out.println(" --> Attention only 5th, 4th and 3rd entry can be set to 0. ");
	}

	
	public static void main(String[] args) {
		

		MultiBins Test = new MultiBins(4,4,5,5,100);
		Test.SetBins(0, 10, 0, 1, 0, 1, 0, 1, 0.1, 0.5);
		Test.GenerateHistograms("tests");
		Test.InizializeClasPolygons();
		double testQ2 =1.7;
		double testXb = 0.12; 
		int BinX= Test.ComputePolygonalBins(testQ2,testXb).get(0);
		int BinY = Test.ComputePolygonalBins(testQ2, testXb).get(1);
		//int BinX= Test.ComputePolygonalBins(11, 0.8).get(0);
		//int BinY = Test.ComputePolygonalBins(11, 0.8).get(1);

	System.out.println(  " Bin Q2 " + BinX + " Bin x "+BinY);

		
		
		MultiBins Test2 = new MultiBins(5,6,5,5,5);
		 double[] ax1 = new double[]{1,3,5,7,10};
		 double[] ax2 = new double[]{1,3,5,7,10};
		 double[] ax3 = new double[]{1,3,5,7,10,20};
		 double[] ax4 = new double[]{1,3,5,7,10,20};
		 double[] ax5 = new double[]{1,3,5,7,10,20};
		Test2.SetUnevenBins(ax1,ax2,ax3,ax4,ax5);
		double Variable = 8;
		
		System.out.println(" My bin is " +Test2.ComputeBin(4, 8));
		/*
		// TODO Auto-generated method stub
		//Let's say we want to study invariant mass of photons in bin of Q2,Xb,Z,Pt . So in total we have 5 variables.
		// Assuming that we want 4 bins in Q2
		// Assuming that we want 4 bins in Xb
		// Assuming that we wants 5 bins in Z
		// Assuming that we wants 5 bins in Pt
		// Assuming we wants 100 bins in the invariant mass 
		MultiBins Test = new MultiBins(4,4,5,5,100);
		// Now we set the limits of each bin, we start with Q2 min and max etc.
		Test.SetBins(0, 10, 0, 1, 0, 1, 0, 1, 0.1, 0.5);
		// Now that we have bins number and limits I can generate the histograms
		Test.GenerateHistograms();
		// In my code I will fill the histogram. Since we have 5 variables, the first 2 are in a double array the second 3 are in an histo3F
		Test.GetH3F(0, 0).fill(3,2,1);
		Test.GetH3F(0, 0).getXAxis().setTitle("Prova");
		System.out.println(" THe histogram title is: " + Test.GetH3F(0,0).getXAxis().getTitle());
		// As you can see we can access to the histogram I need
		double Q2=2.32;
		double XB=0.41;
		Test.GetH3F(Test.ComputeBin(1, Q2), Test.ComputeBin(2, XB)).getXAxis().setTitle("Q2:2.32/Xb:0.41");
		System.out.println(" The histogram in bin Q2: " + Q2 +" and XB " +XB);
		System.out.println(" The title of that H3F is " + Test.GetH3F(Test.ComputeBin(1, Q2), Test.ComputeBin(2, XB)).getXAxis().getTitle());
// 
		MultiBins Test4D =  new MultiBins(4,4,5,5,0);
		Test4D.SetBins(0, 10, 0, 1, 0, 1, 0, 1, 0, 0);
		Test4D.GenerateHistograms();
        Test4D.GetH2F(2, 2).fill(0.1, 0.1);
        Test4D.GetH2F(2, 2).fill(0.1, 0.1);
        System.out.println(" Third bin center" + Test4D.GetH2F(2, 2).getXAxis().getBinCenter(2));
        System.out.println("The bin 1-1 of the H2F contains " + Test4D.GetH2F(2, 2).getBinContent(1, 1) + " counts, while the bin 1-2 contains" + Test4D.GetH2F(2, 2).getBinContent(1, 2) );


        MultiBins Test3D = new MultiBins(5,5,10,0,0);
        Test3D.SetBins(0, 1, 0, 1, 0, 0.5, 0, 0, 0, 0);
        Test3D.GenerateHistograms();
        System.out.println(" Test on the 3D dimension object max size" + Test3D.GetH1F(0, 0).getMaximumBin());
        */
	}
	
	public void SetPolygonal(boolean polygon) {
		this.PolygonalBinClas=polygon;
	}
	/**
	 * This functions set the limits of each variable 
	 * @param First_min
	 * @param First_max
	 * @param Second_min
	 * @param Second_max
	 * @param Third_min
	 * @param Third_max
	 * @param Forth_min
	 * @param Forth_max
	 */
	
	
	public void SetUnevenBins(double[]b1, double[]b2, double[]b3,double[]b4, double[]b5) {
		even=false;
		if(this.Dimension==5) {
		this.bin1=b1;
		this.bin2=b2;
		this.bin3=b3;
		this.bin4=b4;
		this.bin5=b5;
		}
		else if(this.Dimension==4) {
			this.bin1=b1;
			this.bin2=b2;
			this.bin3=b3;
			this.bin4=b4;
			}
		else if(this.Dimension==3) {
			this.bin1=b1;
			this.bin2=b2;
			this.bin3=b3;
			}
		else if(this.Dimension==2) {
			this.bin1=b1;
		//System.out.println("-->> Bins 1 "+bin1.length);	
		
			this.bin2=b2;
			}

	}
	
	public double[] GetUnevenBins(int i) {
		if(i==0)System.out.println(" Return bin uneven needs to start from 1");
		if(even==true) System.out.println(" Attention your are try to get uneven bins when you have even seaprated bins");
		if(i==2) return this.bin2;
		if(i==3) return this.bin3;
		if(i==4) return this.bin4;
		if(i==5) return this.bin5;
		else return this.bin1;
	}
	public void SetBins(double First_min, double First_max, double Second_min, double Second_max, double Third_min, double Third_max, double Forth_min, double Forth_max, double Fifth_min, double Fifth_max  )
	{
		even=true;
		this.Fst_Min=First_min;
		this.Fst_Max=First_max;
		this.Snd_Min=Second_min;
		this.Snd_Max=Second_max;
		this.Trd_Min=Third_min;
		this.Trd_Max=Third_max;
		this.Frt_Min=Forth_min;
		this.Frt_Max=Forth_max;
		this.Fth_Min=Fifth_min;
		this.Fth_Max=Fifth_max;
	}

	/**
	 * This class generates the Histograms needed.
	 * Remember to execute .SetBins() first so define the boundary of each histogram 
	 */
	public void GenerateHistograms(String Name) {
    this.HistoName=Name;
	if(even==true) {	
    if(this.Dimension==5) {
			for(int qq=0 ; qq<FstBin; qq++) {
				MultiHisto3D.add(new ArrayList<H3F>());
				for (int kk=0 ; kk<SndBin; kk ++) {
					MultiHisto3D.get(qq).add( new H3F(TrdBin, Trd_Min, Trd_Max, FrtBin, Frt_Min, Frt_Max, FthBin, Fth_Min, Fth_Max)); 
				
				}
			}
		}
		else if (this.Dimension==4){
			for(int qq=0 ; qq<FstBin; qq++) {
				MultiHisto2D.add(new ArrayList<H2F>());
				for (int kk=0 ; kk<SndBin; kk ++) {
					H2F histo = new H2F(Name,TrdBin, Trd_Min, Trd_Max,FrtBin, Frt_Min, Frt_Max);
					MultiHisto2D.get(qq).add(histo);
					//System.out.println(" This Bin" + this.TrdBin);
				// System.out.println("Bin 2 center in the loop " +histo.getXAxis().getBinCenter(2));
				}
			}
		}
		else if (this.Dimension==3) {
			for(int qq=0 ; qq<FstBin; qq++) {
				MultiHisto1D.add(new ArrayList<H1F>());
				for (int kk=0 ; kk<SndBin; kk ++) {
					MultiHisto1D.get(qq).add(new H1F(Name,TrdBin, Trd_Min, Trd_Max));
				}
			}

		}
		
		else if (this.Dimension==2) {
			this.Counts2D=new H2F(Name,FstBin,Fst_Min,Fst_Max,SndBin,Snd_Min,Snd_Max);
			System.out.println("creating an histo 2D");
			
		}
	}
	else if (even==false) {
		if(this.Dimension==5) {
			for(int qq=0 ; qq<FstBin; qq++) {
				MultiHisto3D.add(new ArrayList<H3F>());
				for (int kk=0 ; kk<SndBin; kk ++) {
					// Adding vector axis to the hitso
					MultiHisto3D.get(qq).add( new H3F(bin3,bin4,bin5)); 
				}
			}
		}
		else if (this.Dimension==4){
			for(int qq=0 ; qq<FstBin; qq++) {
				MultiHisto2D.add(new ArrayList<H2F>());
				for (int kk=0 ; kk<SndBin; kk ++) {
					// Adding vector axis to the hitso
					H2F histo = new H2F(Name,bin3,bin4);
					MultiHisto2D.get(qq).add(histo);
					//System.out.println(" This Bin" + this.TrdBin);
				// System.out.println("Bin 2 center in the loop " +histo.getXAxis().getBinCenter(2));
				}
			}
		}
		else if (this.Dimension==3) {
			for(int qq=0 ; qq<FstBin; qq++) {
				MultiHisto1D.add(new ArrayList<H1F>());
				for (int kk=0 ; kk<SndBin; kk ++) {
					// Adding vector axis to the hitso
					MultiHisto1D.get(qq).add(new H1F(Name,bin3));
				}
			}

		}
		
		else if (this.Dimension==2) {
			// new bin in X and Y : 
			//System.out.println("^^^^^ 1st bins "+ this.FstBin +" ^^^% 2nd bins "+ this.SndBin);
			double []bin1new= new double[this.FstBin];
			for(int i=0; i<this.FstBin; i++) {
				double binsize =(bin1[1]-bin1[0])/this.FstBin ;
				double minvalue=bin1[0];
				double conto=minvalue+binsize*i;
				//System.out.println(" conto "+conto +" minvalue "+minvalue + " binsize "+binsize);
				bin1new[i]=conto;
			}
			

			double []bin2new= new double[this.SndBin];
			for(int i=0; i<this.SndBin; i++) {
				double binsize =(bin2[1]-bin2[0])/this.SndBin ;
				double minvalue=bin2[0];
				double conto=minvalue+binsize*i;
				//System.out.println(" conto "+conto +" minvalue "+minvalue + " binsize "+binsize);
				bin2new[i]=conto;
			}
			
			
			this.Counts2D=new H2F(Name,bin1new,bin2new);
			//System.out.println("creating an histo 2D with uneven bins min. Nr of bins "+bin1new.length);
			//System.out.println("creating an histo 2D with uneven bins min. Nr of bins "+bin2new.length);

		}
		
	}
	}
	
	public void FillH2(double qval, double xval) {
		
		int BinX,BinY;
		if(this.PolygonalBinClas==true) {
			BinX= ComputePolygonalBins(qval, xval).get(0);
			BinY = ComputePolygonalBins(qval, xval).get(1);
			}
			else {
				BinX = ComputeBin(1, qval);
				BinY = ComputeBin(2,xval);
			}
		
		// Convert the shaped bins into standard bin equivalents:
		//
		
		
		double Min =0 , Max=0, Bins =0;
	
		if(even==true) {
			double XMax=this.Fst_Max;
			double XMin =this.Fst_Min;
			int XBins = this.FstBin;
			double XDimension= XMax - XMin;
			double XSize= XDimension/XBins;
			double XValue = XMin+(XSize/2)+BinX*XSize;
			
		
		
	
			double YMax=this.Snd_Max;
			double YMin =this.Snd_Min;
			int YBins = this.SndBin;
			double YDimension = YMax - YMin;
			double YSize= YDimension/YBins;
			double YValue = YMin+(YSize/2)+BinY*YSize;
			
			//System.out.println( "FillH2 in Bins , BinX from polygon "+ BinX + "Valore per histogramma" + XValue );
		this.Counts2D.fill(XValue,YValue);
		//Counts2D.getBinContent(0, 1);
		}
		else {

			double XMax=this.bin1[this.bin1.length-1];
			double XMin =this.bin1[0];
			int XBins = this.FstBin;
			double XDimension= XMax - XMin;
			double XSize= XDimension/XBins;
			double XValue = XMin+(XSize/2)+BinX*XSize;
			
			
			double YMax=this.bin2[this.bin2.length-1];
			double YMin =this.bin2[0];
			int YBins = this.SndBin;
			double YDimension = YMax - YMin;
			double YSize= YDimension/YBins;
			double YValue = YMin+(YSize/2)+BinY*YSize;
		//if(BinY>1) {
			//System.out.println("Valore minimo x "+XMin + " MAX "+ XMax + " BinX "+ BinX);
			//System.out.println("Valore minimo y "+YMin + " MAX "+ YMax+ "Bin Y "+ BinY);
			//System.out.println("qvale " + qval +  "Bin X from polygon "+ BinX + "Valore per histogramma" + XValue );
			//System.out.println("xval " +xval+ " Bin y from polygon "+ BinY + "Valore per histogramma" + YValue );
		   // System.out.println("XALUE X MIN " +this.Counts2D.getXAxis().min() + " MAX " +this.Counts2D.getXAxis().max()+" BINS X "+this.Counts2D.getXAxis().getNBins()+ " Y MIN " +this.Counts2D.getYAxis().min() +" Y MAX "+this.Counts2D.getYAxis().max()+" BINS Y "+this.Counts2D.getYAxis().getNBins());
		//}
		
		this.Counts2D.fill(XValue,YValue);
		}
		
	}
	
	public H2F GetH2() {
		return this.Counts2D;
	}
	/**
	 * 
	 * @param i the bin of the first variable
	 * @param j the bin for the second variable
	 * @return the Histo3D that occupies the ith-jth position in the array of histograms
	 */
	public H3F GetH3F(int i, int j) {
		if (this.Dimension==5) {
			return MultiHisto3D.get(i).get(j);
		}
		else {
			System.out.println(" ** Waring, You created Bins of dimensions different from 5" );
			return null; 
		}
	}

	/**
	 * 
	 * @param i the bin of the first variable
	 * @param j the bin for the second variable
	 * @return the Histo2D that occupies the ith-jth position in the array of histograms
	 */
	public H2F GetH2F(int i, int j) {
		if (this.Dimension==4) {
			return MultiHisto2D.get(i).get(j);
		}
		else {
			System.out.println(" ** Waring, You created Bins of dimensions different than 4 " );
			return null; 
		}
	}

	/**
	 * 
	 * @param i the bin of the first variable
	 * @param j the bin for the second variable
	 * @return the Histo3D that occupies the ith-jth position in the array of histograms
	 */
	public H1F GetH1F(int i, int j) {
		if (this.Dimension==3) {
			return MultiHisto1D.get(i).get(j);
		}
		else {
			System.out.println(" ** Waring, You created Bins of dimensions<5" );
			return null; 
		}
	}

	/**
	 * The function return the bin number for a value
	 * @param q the number representing the variable 
	 * @param value the value of the variable
	 * @return the bin 
	 */
	public int ComputeBin(int q, double value) {
		int bin=-1;
		if(even==true) {
		double Min =0 , Max=0, Bins =0;
		if (q==1) {
			Max=this.Fst_Max;
			Min =this.Fst_Min;
			Bins = this.FstBin;
		}
		else if (q==2) {
			Min= this.Snd_Min;
			Max =this.Snd_Max;
			Bins = this.SndBin;
		}
		else if (q==3) {
			Min = this.Trd_Min;
			Max = this.Trd_Max;
			Bins = this.TrdBin;
		}
		else if (q==4 && this.Dimension>3) {
			Min = this.Frt_Min;
			Max = this.Frt_Max;
			Bins = this.FrtBin;
		}

		else if (q==5 && this.Dimension>4) {
			Min = this.Fth_Min;
			Max = this.Fth_Max;
			Bins = this.FthBin;
		}
		double Dimension= Max - Min;
		double Size= Dimension/Bins;
		double positionmin=Min;
		double positionmax=Min+Size;
	
		for(int i=0;  i<Bins; i++ ) {
			positionmin=Min+i*Size;
			positionmax= Min+(i+1)*Size;
			if(value<=positionmax   && value>=positionmin) {
				bin=i;
				break;
			}
		}
		return bin;
		}
		
		else {
		
			double []arraybin = null;
			if (q==1) {
				arraybin=this.bin1;
			}
			else if (q==2) {
				arraybin=this.bin2;
			}
			else if (q==3) {
				arraybin=this.bin3;
			}
			else if (q==4 && this.Dimension>3) {
				arraybin=this.bin4;
			}

			else if (q==5 && this.Dimension>4) {
				arraybin=this.bin5;
			}
			for(int i=1; i<arraybin.length; i++) {
			//	System.out.println("-> MB:  my array length is "+ arraybin.length + " arraybin left "+arraybin[i-1] + " right "+ arraybin[i]  );
				//System.out.println("-> MB:  Value is "+ value  +" MB dimension" + this.Dimension + " Even? "+ even);
				if(value<=arraybin[i]   && value>=arraybin[i-1]) {
					bin=i-1;
					break;
				}
			}
			return bin;
		}
		
		
	}
	
	private List<Integer> MapPolygons(int i) {
		/*
		// First variable - Second Variable 
		List<Integer> Indici = new ArrayList<Integer>();
		//System.out.println("  MapPoly " + i);
		if (i==0) {Indici.add(0, 0); Indici.add(1, 0);}
		if (i==0) {Indici.add(0, 0); Indici.add(1, 0);}
		else if (i==1) {Indici.add(0, 1); Indici.add(1, 0);}
		else if (i==2) {Indici.add(0,2); Indici.add(1, 0);}
		
		else if (i==3) {Indici.add(0, 0); Indici.add(1, 1);}
		else if (i==4) {Indici.add(0, 1); Indici.add(1, 1);}
		else if (i==5) {Indici.add(0,2); Indici.add(1, 1);}
		
		
		else if (i==6) {Indici.add(0, 0); Indici.add(1, 2);}
		else if (i==7) {Indici.add(0, 1); Indici.add(1, 2);}
		else if (i==8) {Indici.add(0,2); Indici.add(1, 2);}
		
		return Indici;
		*/
	
		List<Integer> Indici = new ArrayList<Integer>();
		// 0 is Q2 and 1 is X 
		if (i==0) {Indici.add(0, 0); Indici.add(1, 0);}
		else if (i==1) {Indici.add(0, 1); Indici.add(1, 0);}
		else if (i==2) {Indici.add(0,1); Indici.add(1, 1);}
		else if (i==3) {Indici.add(0,1); Indici.add(1, 2);}
		else if (i==4) {Indici.add(0,2); Indici.add(1, 1);}
		else if (i==5) {Indici.add(0,2); Indici.add(1, 2);}
		else if (i==6) {Indici.add(0,2); Indici.add(1, 3);}
		else if (i==7) {Indici.add(0,2); Indici.add(1, 4);}
		else if (i==8) {Indici.add(0,3); Indici.add(1, 2);}
		else if (i==9) {Indici.add(0,3); Indici.add(1, 3);}
		else if (i==10) {Indici.add(0,3); Indici.add(1, 4);}
		else if (i==11) {Indici.add(0,4); Indici.add(1, 3);}
		else if (i==12) {Indici.add(0,4); Indici.add(1, 4);}
		else if (i==13) {Indici.add(0,5); Indici.add(1, 4);}
		
		/*
		// 0 is Q2 and 1 is X 
		if (i==0) {Indici.add(0, 0); Indici.add(1, 0);}
		else if (i==1) {Indici.add(0, 1); Indici.add(1, 0);}
		else if (i==2) {Indici.add(0,0); Indici.add(1, 1);}
		else if (i==3) {Indici.add(0,1); Indici.add(1, 1);}
		else if (i==4) {Indici.add(0,1); Indici.add(1, 2);}
		else if (i==5) {Indici.add(0,0); Indici.add(1, 3);}
		else if (i==6) {Indici.add(0,1); Indici.add(1, 3);}
		else if (i==7) {Indici.add(0,0); Indici.add(1, 4);}
		else if (i==8) {Indici.add(0,1); Indici.add(1, 4);}
		*/
		
		
		
		
		// OLD with 16 bins 
		/*
		if (i==0) {Indici.add(0, 0); Indici.add(1, 0);}
		else if (i==1) {Indici.add(0, 0); Indici.add(1, 1);}
		else if (i==2) {Indici.add(0,1); Indici.add(1, 1);}
		else if (i==3) {Indici.add(0,0); Indici.add(1, 2);}
		else if (i==4) {Indici.add(0,1); Indici.add(1, 2);}
		else if (i==5) {Indici.add(0,0); Indici.add(1, 3);}
		else if (i==6) {Indici.add(0,1); Indici.add(1, 3);}
		else if (i==7) {Indici.add(0,2); Indici.add(1, 3);}
		else if (i==8) {Indici.add(0,0); Indici.add(1, 4);}
		else if (i==9) {Indici.add(0,1); Indici.add(1, 4);}
		else if (i==10) {Indici.add(0,2); Indici.add(1, 4);}
		else if (i==11) {Indici.add(0,0); Indici.add(1, 5);}
		else if (i==12) {Indici.add(0,1); Indici.add(1, 5);}
		else if (i==13) {Indici.add(0,2); Indici.add(1, 5);}
		else if (i==14) {Indici.add(0,0); Indici.add(1, 6);}
		else if (i==15) {Indici.add(0,1 ); Indici.add(1, 6);}
*/
		return Indici;
		

	}
	public List<Integer> ComputePolygonalBins(double qval, double xval) {
		List<Integer> Indici = new ArrayList<Integer>();
		int indice=0;
     
	//	System.out.println(" ---- polygon size is  " + Polygons.size());
	for(int i=0; i<this.Polygons.size(); i++) {
		
			
		if ( Polygons.get(i).contains(qval,xval)==true) { 
		
			Indici =MapPolygons(i);
		//	System.out.println(" Got my values q2 "+Indici.get(0)+" xb "+Indici.get(1)  );
			return Indici;
		}
	}
	//System.out.println("->> Il BIN Non l ho trovato, il valore e' " +qval +" " + xval);
	// Overfloar or over Q2 in bin 5 , overfloat or under X in bin 6; 
	//Indici.add(0,53);
	Indici.add(0,5);
	Indici.add(1,6);
	return Indici;
		
	}
	/**
	 * The function print historams pictures
	 * @param time
	 */
	public void Print(String time) {
		String s1="/Users/gangelini/Desktop/";
		new File(s1+"/Plots").mkdir();
		String PathFolder = s1+"/Plots/"+time;
		File dir = new File(PathFolder);
		dir.mkdir();
		PathFolder=s1+"/Plots/"+time+"/MultiBins";
		new File (PathFolder).mkdir();
		System.out.println(" (Pi) - > I am Plotting Multi Bins Histograms");
        for (int j=0; j<this.FstBin ; j++)
        {
        	String number = Integer.toString(j);
    		PathFolder=s1+"/Plots/"+time+"/MultiBins/First_"+number;
    		File dirbin = new File(PathFolder);
    		dirbin.mkdir();
    		for(int k=0 ; k<this.SndBin; k++ ) {
    		String number2=Integer.toString(k);
        	PathFolder=s1+"/Plots/"+time+"/MultiBins/First_"+number+"/Second_"+number2;
        	File dirbin2 = new File(PathFolder);
        	dirbin2.mkdir();
        	
    		EmbeddedCanvas photons = new EmbeddedCanvas();
    		photons.setSize(1600,1000);
    		photons.divide(this.TrdBin,this.FrtBin);
    		photons.setAxisTitleSize(20);
    		photons.setAxisFontSize(20);
    		photons.setTitleSize(20);
    		
    		for(int i=0 ; i<(this.TrdBin*this.FrtBin); i++) {
    			photons.cd(i);
    			if(this.Dimension==3) {
				photons.draw(this.GetH1F(j, k));
    			}
    			}	   
    		    String strg0 = String.format("%s/Pi0s.png",PathFolder);
    			System.out.println("Saving plots in "+ PathFolder + " For 1st bin "+ j + " Second bin "+k);
    			photons.save(strg0);	
    		}
    		
       }
		
		// I plot z from  from 0.1 to 0.9 because first and last bin are useless
	}
	
	public void SetName_1stVariable( String Nome){this.Name1st = Nome;}
	public void SetName_2ndVariable(String Nome) {this.Name2nd = Nome;}
	public void SetName_3rdVariable(String Nome) {this.Name3rd = Nome;}
	public void SetName_4thVariable(String Nome) {this.Name4th = Nome;}
	public void SetName_5thVariable(String Nome) {this.Name5th = Nome;}
	
	public String GetName_1stVariable() { return this.Name1st;}
	public String GetName_2ndVariable() { return this.Name2nd;}
	public String GetName_3rdVariable() { return this.Name3rd;}
	public String GetName_4thVariable() { return this.Name4th;}
	public String GetName_5thVariable() { return this.Name5th;}
	
public String GetNameHistos() {
		
		// TODO Auto-generated method stub
		return this.HistoName;
	}


	
	

	public void InizializeClasPolygons() {
		// Fill the path from top right by going clockwise. 
// The initialization of standard binning for clas. 		
		// I am making different polygons for test. 
	//System.out.println(" I am inizialazing Polygons for the istograms " + this.HistoName)  ;
		
		

		// TO MODIFY IN HERE ! 
		   this.PolygonNr=14;
		     // I am trying to define squared bins to compare with what was happening before. 
		     Path2D poly1 = new Path2D.Double();
		     
		     Double valoresX[] = {2.00, 1.10, 1.10,2.0};
		     Double valoresY[] = {0.21,0.085,0.069,0.130};

		     poly1.moveTo(valoresX[0], valoresY[0]);
		     for(int i = 1; i < valoresX.length; ++i) {
		        poly1.lineTo(valoresX[i], valoresY[i]);
		     }
		     poly1.closePath();
		     
		     Path2D poly2 = new Path2D.Double();
		     
		     Double valores2X[] = {2.6, 2.0, 2.0};
		     Double valores2Y[] = {0.17,0.17,0.13};
		     
		     poly2.moveTo(valores2X[0], valores2Y[0]);
		     for(int i = 1; i < valores2X.length; ++i) {
		        poly2.lineTo(valores2X[i], valores2Y[i]);
		     }
		     poly2.closePath();
		
		     Path2D poly3 = new Path2D.Double();
		     
		     Double valores3X[] = {2.6,2.0,2.0,2.6};
		     Double valores3Y[] = {0.21,0.21,0.17,0.17};
		    
		     poly3.moveTo(valores3X[0], valores3Y[0]);
		     for(int i = 1; i < valores3X.length; ++i) {
		        poly3.lineTo(valores3X[i], valores3Y[i]);
		     }
		     poly3.closePath();
		     this.Polygons.add(poly1);
		     this.Polygons.add(poly2);
		     this.Polygons.add(poly3);
		     // -----
		     //here
		     Path2D poly4 = new Path2D.Double();
		     Double valores4X[] = {2.6,2.0,2.6};
		     Double valores4Y[] = {0.34,0.21,0.21};
		   
		     poly4.moveTo(valores4X[0], valores4Y[0]);
		     for(int i = 1; i < valores4X.length; ++i) {
		        poly4.lineTo(valores4X[i], valores4Y[i]);
		     }
		     poly4.closePath();
		     
		     Path2D poly5 = new Path2D.Double();
		     Double valores5X[] = {3.3,2.6,2.6};
		     Double valores5Y[] = {0.215,0.21,0.17};
		     poly5.moveTo(valores5X[0], valores5Y[0]);
		     for(int i = 1; i < valores5X.length; ++i) {
		        poly5.lineTo(valores5X[i], valores5Y[i]);
		     }
		     poly5.closePath();
		     
		     Path2D poly6 = new Path2D.Double();
		     Double valores6X[] = {3.3,2.6,2.6,3.3};
		     Double valores6Y[] = {0.26,0.26,0.21,0.215};
		     poly6.moveTo(valores6X[0], valores6Y[0]);
		     for(int i = 1; i < valores6X.length; ++i) {
		        poly6.lineTo(valores6X[i], valores6Y[i]);
		     }
		     poly6.closePath();
		   //  System.out.println( "Does Poly 6 contains 9, 0.5  ?" + poly6.contains(8.4,0.31) );
		     this.Polygons.add(poly4);
		     this.Polygons.add(poly5);
		     this.Polygons.add(poly6);
		  // -----
		     Path2D poly7 = new Path2D.Double();
		     Double valores7X[] =  {3.3,2.6,2.6,3.3};
		     Double valores7Y[] = {0.30,0.30,0.26,0.26};
		
		     poly7.moveTo(valores7X[0], valores7Y[0]);
		     for(int i = 1; i < valores7X.length; ++i) {
		        poly7.lineTo(valores7X[i], valores7Y[i]);
		     }
		     poly7.closePath();
		     
		     
		     Path2D poly8 = new Path2D.Double();
		     Double valores8X[] = {3.3,2.6,2.6,3.3};
		     Double valores8Y[] = {0.47,0.34,0.30,0.3};

		     poly8.moveTo(valores8X[0], valores8Y[0]);
		     for(int i = 1; i < valores8X.length; ++i) {
		        poly8.lineTo(valores8X[i], valores8Y[i]);
		     }
		     poly8.closePath();
		     
		     
		     Path2D poly9 = new Path2D.Double();
		     Double valores9X[] = {4.6,3.3,3.3};
		     Double valores9Y[] = {0.305,0.30,0.215};

		     poly9.moveTo(valores9X[0], valores9Y[0]);
		     for(int i = 1; i < valores9X.length; ++i) {
		        poly9.lineTo(valores9X[i], valores9Y[i]);
		     }
		     poly9.closePath();
		     
	   
		     this.Polygons.add(poly7);
		     this.Polygons.add(poly8);
		     this.Polygons.add(poly9);
			
		   
		     Path2D poly10 = new Path2D.Double();
		     Double valores10X[] = {4.6,3.3,3.3,4.6};
		     Double valores10Y[] = {0.38,0.38,0.30,0.30};

		     poly10.moveTo(valores10X[0], valores10Y[0]);
		     for(int i = 1; i < valores10X.length; ++i) {
		        poly10.lineTo(valores10X[i], valores10Y[i]);
		     }
		     poly10.closePath();
		     
		     Path2D poly11 = new Path2D.Double();
		     Double valores11X[] = {4.6,3.3,3.3,4.6};
		     Double valores11Y[] = {0.57,0.47,0.38,0.38};

		     poly11.moveTo(valores11X[0], valores11Y[0]);
		     for(int i = 1; i < valores11X.length; ++i) {
		        poly11.lineTo(valores11X[i], valores11Y[i]);
		     }
		     poly11.closePath();
		     
		     Path2D poly12 = new Path2D.Double();
		     Double valores12X[] = {6.4,4.6,4.6,6.4};
		     Double valores12Y[] = {0.45,0.45,0.305,0.421};

		     poly12.moveTo(valores12X[0], valores12Y[0]);
		     for(int i = 1; i < valores12X.length; ++i) {
		        poly12.lineTo(valores12X[i], valores12Y[i]);
		     }
		     poly12.closePath();
		     
		     
		     Path2D poly13 = new Path2D.Double();
		     Double valores13X[] = {6.4,4.6,4.6,6.4};
		     Double valores13Y[] = {0.66,0.57,0.45,0.45};

		     poly13.moveTo(valores13X[0], valores13Y[0]);
		     for(int i = 1; i < valores13X.length; ++i) {
		        poly13.lineTo(valores13X[i], valores13Y[i]);
		     }
		     poly13.closePath();
		     
		     Path2D poly14 = new Path2D.Double();
		     Double valores14X[] = {9.6,7.5,6.4,6.4,7.4};
		     Double valores14Y[] = {0.72,0.71,0.66,0.421,0.5};

		     poly14.moveTo(valores14X[0], valores14Y[0]);
		     for(int i = 1; i < valores14X.length; ++i) {
		        poly14.lineTo(valores14X[i], valores14Y[i]);
		     }
		     poly14.closePath();
		     
		     this.Polygons.add(poly10);
		     this.Polygons.add(poly11);
		     this.Polygons.add(poly12);
		     this.Polygons.add(poly13);
		     this.Polygons.add(poly14);
		
		
		
		// TO MODIFY IN HERE ! 
		/*
		   this.PolygonNr=9;
		     // I am trying to define squared bins to compare with what was happening before. 
		     Path2D poly1 = new Path2D.Double();
		     Double valoresX[] = {2.44, 1.38, 1.3,1.3};
		     Double valoresY[] = {0.15,0.15,0.12,0.078};
		     poly1.moveTo(valoresX[0], valoresY[0]);
		     for(int i = 1; i < valoresX.length; ++i) {
		        poly1.lineTo(valoresX[i], valoresY[i]);
		     }
		     poly1.closePath();
		     
		     Path2D poly2 = new Path2D.Double();
		     Double valores2X[] = {2.75, 1.5, 1.45,1.38,1.98};
		     Double valores2Y[] = {0.24,0.24,0.20,0.15,0.15};
		     
		     poly2.moveTo(valores2X[0], valores2Y[0]);
		     for(int i = 1; i < valores2X.length; ++i) {
		        poly2.lineTo(valores2X[i], valores2Y[i]);
		     }
		     poly2.closePath();
		
		     Path2D poly3 = new Path2D.Double();
		     Double valores3X[] = {3.88,2.75,1.98,2.44 };
		     Double valores3Y[] = {0.24,0.24,0.15,0.15};
		     poly3.moveTo(valores3X[0], valores3Y[0]);
		     for(int i = 1; i < valores3X.length; ++i) {
		        poly3.lineTo(valores3X[i], valores3Y[i]);
		     }
		     poly3.closePath();
		     this.Polygons.add(poly1);
		     this.Polygons.add(poly2);
		     this.Polygons.add(poly3);
		     // -----
		     //here
		     Path2D poly4 = new Path2D.Double();
		     Double valores4X[] = {3.63,1.6,1.5,2.75};
		     Double valores4Y[] = {0.34,0.34,0.24,0.24};
		   
		     poly4.moveTo(valores4X[0], valores4Y[0]);
		     for(int i = 1; i < valores4X.length; ++i) {
		        poly4.lineTo(valores4X[i], valores4Y[i]);
		     }
		     poly4.closePath();
		     
		     Path2D poly5 = new Path2D.Double();
		     Double valores5X[] = {5.49,3.63,2.75,3.88};
		     Double valores5Y[] = {0.34,0.34,0.24,0.24};
		     poly5.moveTo(valores5X[0], valores5Y[0]);
		     for(int i = 1; i < valores5X.length; ++i) {
		        poly5.lineTo(valores5X[i], valores5Y[i]);
		     }
		     poly5.closePath();
		     
		     Path2D poly6 = new Path2D.Double();
		     Double valores6X[] = {4.7,2.52,1.6,3.63};
		     Double valores6Y[] = {0.45,0.45,0.34,0.34};
		     poly6.moveTo(valores6X[0], valores6Y[0]);
		     for(int i = 1; i < valores6X.length; ++i) {
		        poly6.lineTo(valores6X[i], valores6Y[i]);
		     }
		     poly6.closePath();
		   //  System.out.println( "Does Poly 6 contains 9, 0.5  ?" + poly6.contains(8.4,0.31) );
		     this.Polygons.add(poly4);
		     this.Polygons.add(poly5);
		     this.Polygons.add(poly6);
		  // -----
		     Path2D poly7 = new Path2D.Double();
		     Double valores7X[] =  {7.22,4.7,3.63,5.49};
		     Double valores7Y[] = {0.45,0.45,0.34,0.34};
		
		     poly7.moveTo(valores7X[0], valores7Y[0]);
		     for(int i = 1; i < valores7X.length; ++i) {
		        poly7.lineTo(valores7X[i], valores7Y[i]);
		     }
		     poly7.closePath();
		     
		     
		     Path2D poly8 = new Path2D.Double();
		     Double valores8X[] = {7.42,5.4,4.05,3.05,2.52,4.7};
		     Double valores8Y[] = {0.708,0.64,0.57,0.5,0.45,0.45};

		     poly8.moveTo(valores8X[0], valores8Y[0]);
		     for(int i = 1; i < valores8X.length; ++i) {
		        poly8.lineTo(valores8X[i], valores8Y[i]);
		     }
		     poly8.closePath();
		     
		     
		     Path2D poly9 = new Path2D.Double();
		     Double valores9X[] = {11.5,9.25,7.42,4.7,7.22,9.3};
		     Double valores9Y[] = {0.79,0.75,0.708,0.45,0.45,0.58};

		     poly9.moveTo(valores9X[0], valores9Y[0]);
		     for(int i = 1; i < valores9X.length; ++i) {
		        poly9.lineTo(valores9X[i], valores9Y[i]);
		     }
		     poly9.closePath();
		     
	   
		     this.Polygons.add(poly7);
		     this.Polygons.add(poly8);
		     this.Polygons.add(poly9);
			
		   
		
		*/
		     double ok; 
		     // Questo era il vecchio 
		     
		     
		/*
	     this.PolygonNr=9;
	     // I am trying to define squared bins to compare with what was happening before. 
	     Path2D poly1 = new Path2D.Double();
	     Double valoresX[] = {4.66666666667, 1.0, 1.0, 4.66666666667};
	     Double valoresY[] = {0.3066666666666667, 0.3066666666666667, 0.01,0.01 };
	
	     poly1.moveTo(valoresX[0], valoresY[0]);
	     for(int i = 1; i < valoresX.length; ++i) {
	        poly1.lineTo(valoresX[i], valoresY[i]);
	     }
	     poly1.closePath();
	     Path2D poly2 = new Path2D.Double();
	     Double valores2X[] = { 8.33333336667,  4.66666666667, 4.66666666667,8.33333336667,};
	     Double valores2Y[] = {0.3066666666666667, 0.3066666666666667, 0.01, 0.01 };
	     poly2.moveTo(valores2X[0], valores2Y[0]);
	     for(int i = 1; i < valores2X.length; ++i) {
	        poly2.lineTo(valores2X[i], valores2Y[i]);
	     }
	     poly2.closePath();
	
	     Path2D poly3 = new Path2D.Double();
	     Double valores3X[] = {12.0,8.33333336667,  8.33333336667,12.0  };
	     Double valores3Y[] = {  0.3066666666666667,0.3066666666666667,0.01,0.01};
	     poly3.moveTo(valores3X[0], valores3Y[0]);
	     for(int i = 1; i < valores3X.length; ++i) {
	        poly3.lineTo(valores3X[i], valores3Y[i]);
	     }
	     poly3.closePath();
	     this.Polygons.add(poly1);
	     this.Polygons.add(poly2);
	     this.Polygons.add(poly3);
	     // -----
	     //here
	     Path2D poly4 = new Path2D.Double();
	     Double valores4X[] =  { 4.66666666667, 1.0,1.0, 4.66666666667};
	     Double valores4Y[] = {0.6033333333333334, 0.6033333333333334, 0.3066666666666667,  0.3066666666666667};
	   
	     poly4.moveTo(valores4X[0], valores4Y[0]);
	     for(int i = 1; i < valores4X.length; ++i) {
	        poly4.lineTo(valores4X[i], valores4Y[i]);
	     }
	     poly4.closePath();
	     
	     Path2D poly5 = new Path2D.Double();
	     Double valores5X[] = {  8.33333336667, 4.66666666667,4.66666666667, 8.33333336667};
	     Double valores5Y[] = {0.6033333333333334,0.6033333333333334,  0.3066666666666667, 0.3066666666666667};
	     poly5.moveTo(valores5X[0], valores5Y[0]);
	     for(int i = 1; i < valores5X.length; ++i) {
	        poly5.lineTo(valores5X[i], valores5Y[i]);
	     }
	     poly5.closePath();
	     
	     Path2D poly6 = new Path2D.Double();
	     Double valores6X[] = {12.0, 8.33333336667,  8.33333336667,12.0};
	     Double valores6Y[] = {0.6033333333333334, 0.6033333333333334,  0.3066666666666667,0.3066666666666667};
	     poly6.moveTo(valores6X[0], valores6Y[0]);
	     for(int i = 1; i < valores6X.length; ++i) {
	        poly6.lineTo(valores6X[i], valores6Y[i]);
	     }
	     poly6.closePath();
	   //  System.out.println( "Does Poly 6 contains 9, 0.5  ?" + poly6.contains(8.4,0.31) );
	     this.Polygons.add(poly4);
	     this.Polygons.add(poly5);
	     this.Polygons.add(poly6);
	  // -----
	     Path2D poly7 = new Path2D.Double();
	     Double valores7X[] =  { 4.66666666667,  1.0, 1.0,4.66666666667};
	     Double valores7Y[] = {0.9,0.9,0.6033333333333334,  0.6033333333333334};
	
	     poly7.moveTo(valores7X[0], valores7Y[0]);
	     for(int i = 1; i < valores7X.length; ++i) {
	        poly7.lineTo(valores7X[i], valores7Y[i]);
	     }
	     poly7.closePath();
	     
	     
	     Path2D poly8 = new Path2D.Double();
	     Double valores8X[] = {8.33333336667, 4.66666666667, 4.66666666667,8.33333336667};
	     Double valores8Y[] = {0.9, 0.9,0.6033333333333334,0.6033333333333334};

	     poly8.moveTo(valores8X[0], valores8Y[0]);
	     for(int i = 1; i < valores8X.length; ++i) {
	        poly8.lineTo(valores8X[i], valores8Y[i]);
	     }
	     poly8.closePath();
	     
	     
	     Path2D poly9 = new Path2D.Double();
	     Double valores9X[] = { 12.0, 8.33333336667, 8.33333336667,12.0};
	     Double valores9Y[] = {0.9,0.9,0.6033333333333334,0.6033333333333334};
	  
	     poly9.moveTo(valores9X[0], valores9Y[0]);
	     for(int i = 1; i < valores9X.length; ++i) {
	        poly9.lineTo(valores9X[i], valores9Y[i]);
	     }
	     poly9.closePath();
   
	     this.Polygons.add(poly7);
	     this.Polygons.add(poly8);
	     this.Polygons.add(poly9);
		// TODO Auto-generated method stub
		 * 
		 */
		
	}

	
	
}