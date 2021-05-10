/**
 * 
 */
package TIFSeq;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;

import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.util.Log;
import net.sf.picard.util.SolexaQualityConverter;

/**
 * @author wave
 * @date Oct 23, 2012
 */

class mergeRes{
	 static boolean pmerged=false;
	 static boolean uext=false;
	 static int ncsta=-1;
	 static boolean dext=false;
	 static int ncend=-1;
	 static int clusterLab=0;
		
	public mergeRes() {
		// TODO Auto-generated constructor stub
	}
	public mergeRes(boolean pmerged,boolean uext,int ncsta, boolean dext, int ncend,int clusterLab) {
		this.pmerged=pmerged;
		this.uext=uext;
		this.ncsta=ncsta;
		this.dext=dext;
		this.ncend=ncend;
		this.clusterLab=clusterLab;
		// TODO Auto-generated constructor stub
	}
}

public class TScluster extends CommandLineProgram {
    private static final Log LOG = Log.getInstance(TScluster.class);

    @Usage 
    public String USAGE = "\tTScluster V1.1.2.\n" +
    		"\tCluster the transcrption Start or End sites based on the window sizes...\n"+
        		"Usage: \n"
        +"\tjavasol TScluster F1=TSsFile F2=chrLengthFIle outF=clusters.txt \n" +
        		"\tjavasol TScluster F1=ypd_tss.txt F2=chrLength.txt outF=outClusters.txt"
        ;

    @Option(shortName="F1", doc="Input file for transcript start sites, or end sites data.")
    public String tsFile;

    @Option(shortName="F2", doc="Input the length of all chromosomes (contigs) that used....")
    public static String chrLenFile;

    @Option(shortName="outF", doc="output clluster file...")
    public static String outFile;

    @Option(shortName="FT", doc="output clluster file...")
    public static String format="bed";

    @Option(shortName = "WS", doc = "Half window size.. defalt = 2", optional = true)
    public static int hfwinsize=2 ;
    
    @Option(shortName = "DT", doc = "The max distance between two TSs.. default = 5", optional = true)
    public static int disthres=5 ;
    
    @Option(shortName = "OP", doc = "The minimal number of the peak.. defalt = 0", optional = true)
    public static int pthreshold=0 ;
    
    @Option(shortName = "OS", doc = "The minimal number of total number of reads in the cluster window..... default = 0", optional = true)
    public static int sumhold=0 ;
    /** Stock main method. */
    
    public static void main(final String[] argv) {
        System.exit(new TScluster().instanceMain(argv));
        
    }
	
	public static File relDirec=new File("/Volumes/ouzhou/projects/TSES/data/2012-09-09/cluster/");
//	public static String chrLenFile=relDirec +"/"+"chrLen.txt";
//	public static int hfwinsize=20; //half window size...
//	public static int disthres=40; //distance threshold...
	
	public static LinkedHashMap tsclus=new LinkedHashMap(); //cluster object of tss...

//	public static void main(String[] args){
//		TScluster.doTScluster(args);
////		TScluster.test();
//		
//	}
	public static void test(){
		String t5file=relDirec+"/"+"gal_t5_ord.txt";
//		String t5file=relDirec+"/"+"test.txt";
		
		TScluster.doCluster(t5file);
		String outfile=relDirec+"/chrstrand/";
//		TScluster.printCluster(outfile,3);
//		TScluster.printCluster(outfile,5);
		
		LinkedHashMap chrtssArray=TScluster.readInTSsArray(t5file);
		TScluster.pintClusterNumber(outfile,chrtssArray,3,50);
		
//		System.out.println(((int [])chrtssArray.get("chr2_+") )[278388-1]);
//		System.out.println(((int [])chrtssArray.get("chr2_-") )[278388-1]);
	}
	
	

	/* (non-Javadoc)
	 * @see net.sf.picard.cmdline.CommandLineProgram#doWork()
	 */
	@Override
	protected int doWork() {
		// TODO Auto-generated method stub
		TScluster.doCluster(tsFile);
		LinkedHashMap chrtssArray=TScluster.readInTSsArray(tsFile);
		TScluster.pintClusterNumber(outFile,chrtssArray,pthreshold,sumhold);
		
		return 0;
	}
	 
	public static void doTScluster(String[] args){

		
	}
	public static void doCluster(String tsfile){
		
		LinkedHashMap chrtss=TScluster.readInTSs(tsfile);
		
		Iterator it =chrtss.keySet().iterator();
        while(it.hasNext()){
        	
	        	String chrstrand=it.next().toString();
	       
	        	
	        	LinkedHashMap tss=(LinkedHashMap) chrtss.get(chrstrand);
	        	
//	        	if(chrstrand.equalsIgnoreCase("chr2_+")){
	        		LinkedHashMap tsclu=new LinkedHashMap(){};
	        		tsclus.put(chrstrand, tsclu);
	        		TScluster.clusterChr(chrstrand, tss);
	        	 	System.out.println(chrstrand);
	        	 	
	        	 
//	        }
	        
        }		
	}
	
	
	public static void printCluster(String outFileDir,int threshold){
		
		try{
			
			String outFile=outFileDir+"/"+"allCluster_p"+threshold+".bed";
//			String outFile=outFileDir+"/"+"allCluster_w65_p"+threshold+".bed";
			PrintWriter out=new PrintWriter(new BufferedWriter(new FileWriter(outFile)));
			
			Iterator its=tsclus.keySet().iterator();
			System.out.println(tsclus.keySet().size());
			while(its.hasNext()){
				
				String chrstrand=its.next().toString();
				System.out.println(chrstrand + "is printing out....");
				
				String bb[]=chrstrand.split("_");
	    			String chr=bb[0];
	    			String strand=bb[1];
	    			
	    			
//	    			String outFile=outFileDir+"/"+chrstrand+".bed";
//	    			PrintWriter out=new PrintWriter(new BufferedWriter(new FileWriter(outFile)));
//	    			System.out.println(outFile);
	//    			int chrn=Integer.parseInt(chr.substring(3));
		
				LinkedHashMap tsclu=(LinkedHashMap) tsclus.get(chrstrand);
				
				Iterator itt=tsclu.keySet().iterator();
				
				while(itt.hasNext()){
					
					int cluNum= Integer.parseInt( itt.next().toString());
				  	
				  	int[] cluInfo= (int[]) tsclu.get(cluNum);
					int csta=cluInfo[0];
					int cend=cluInfo[1];
					int ctss=cluInfo[2];
					int cnum=cluInfo[3];
					
	//				int chrn=cluInfo[4];
	//				int nstr=cluInfo[5];
//				  	System.out.println(chr+"\t"+(csta-1)+"\t"+cend+"\t"+cluNum+"_"+ctss+"\t"+cnum+"\t"+strand);
					if(cnum>=threshold){
						out.println(chr+"\t"+(csta-1)+"\t"+cend+"\t"+cluNum+"_"+ctss+"\t"+cnum+"\t"+strand);
//				  		out.println(chr+"\t"+(csta-1)+"\t"+cend+"\t"+cluNum+"\t"+ctss+"\t"+cnum+"\t"+strand);
					}
				  	
				}
			  	
//			  	out.close();
			}
			out.close();
			}catch(Exception ex){
				ex.printStackTrace();
				System.out.println("There is problem writing the new sequences.");}
				
	}
	
	public static void pintClusterNumber(String outFile, LinkedHashMap chrtssArray,int pthreshold, int sumhold){
		
		try{
			
//			String outFile=outFileDir+"/"+"allClusterNums_p"+threshold+".bed";
//			String outFile=outFileDir+"/"+"allClusterNums_w65_p"+threshold+".bed";
			PrintWriter out=new PrintWriter(new BufferedWriter(new FileWriter(outFile)));
			
			Iterator its=tsclus.keySet().iterator();
			System.out.println(tsclus.keySet().size());
			while(its.hasNext()){
				
				String chrstrand=its.next().toString();
				System.out.println(chrstrand + "is printing out....");
				
				String bb[]=chrstrand.split("_");
	    			String chr=bb[0];
	    			String strand=bb[1];
	    			
	    			
//	    			String outFile=outFileDir+"/"+chrstrand+".bed";
//	    			PrintWriter out=new PrintWriter(new BufferedWriter(new FileWriter(outFile)));
//	    			System.out.println(outFile);
	//    			int chrn=Integer.parseInt(chr.substring(3));
		
				LinkedHashMap tsclu=(LinkedHashMap) tsclus.get(chrstrand);
				int [] chrTSs= (int []) chrtssArray.get(chrstrand);
				
				Iterator itt=tsclu.keySet().iterator();
				
				while(itt.hasNext()){
					
					int cluNum= Integer.parseInt( itt.next().toString());
				  	
				  	int[] cluInfo= (int[]) tsclu.get(cluNum);
					int csta=cluInfo[0];
					int cend=cluInfo[1];
					int ctss=cluInfo[2];
					int cnum=cluInfo[3];
					
	//				int chrn=cluInfo[4];
	//				int nstr=cluInfo[5];
//				  	System.out.println(chr+"\t"+(csta-1)+"\t"+cend+"\t"+cluNum+"_"+ctss+"\t"+cnum+"\t"+strand);
//				  	out.println(chr+"\t"+(csta-1)+"\t"+cend+"\t"+cluNum+"_"+ctss+"\t"+cnum+"\t"+strand);
					
					int sum=0;
					for(int ii=csta-1;ii<cend;ii++){
						sum=sum+chrTSs[ii];
					}
					
					/*Debug the numbers...*/
					if(chrTSs[ctss-1]!=cnum){
						System.out.println("The corrdinates maybe wrong....");
						System.out.println(chrTSs[ctss-1]+"\t"+cnum);
						
						System.out.println(chr+"\t"+(csta-1)+"\t"+cend+"\t"+cluNum+"\t"+ctss+"\t"+cnum+"\t"+strand+"\t"+sum);
					  	
						}
					
					if(cnum>=pthreshold && sum>=sumhold){
//						out.println(chr+"\t"+(csta-1)+"\t"+cend+"\t"+cluNum+"\t"+ctss+"\t"+cnum+"\t"+strand+"\t"+sum);
						
						if(format.equalsIgnoreCase("bed")){
							out.println(chr+"\t"+(csta-1)+"\t"+cend+"\t"+cluNum+"_"+ctss+"_"+sum+"\t"+cnum+"\t"+strand);	
						}
						if(format.equalsIgnoreCase("tab")){
							out.println(chr+"\t"+csta+"\t"+cend+"\t"+cluNum+"\t"+ctss+"\t"+sum+"\t"+cnum+"\t"+strand);	
						}
					}
				  	
				  	
				  	
				}
			  	
//			  	out.close();
			}
			out.close();
			}catch(Exception ex){
				ex.printStackTrace();
				System.out.println("There is a problem writing the new sequences.");}
				
	}
	
	public static LinkedHashMap readInTSs(String t5file) {
		
		LinkedHashMap chrtss=new LinkedHashMap();
	
//		String prechrstrand="";
		
		try{
			
			BufferedReader br=new BufferedReader(new FileReader(t5file));
			
			while(br.ready()){
				String str=br.readLine();
				String a[]=str.split("\t");
				
				String chr=a[0];
				String strand=a[1];
				int pos=Integer.parseInt(a[2]);
				int num=Integer.parseInt(a[3]);
				
//				if(num>1) {System.out.println(str);}
				
				String chrstrand=chr+"_"+strand;
				if(chrtss.keySet().contains(chrstrand)){
//				if(chrtss.keySet().contains(chrstrand)){
							
					LinkedHashMap tss=(LinkedHashMap) chrtss.get(chrstrand);
					tss.put(chr+"_"+strand+"_"+pos, num);
					
				}else{
					LinkedHashMap tss=new LinkedHashMap();
					chrtss.put(chrstrand,tss);
					tss.put(chr+"_"+strand+"_"+pos, num);	
					
					
				}
				
			}
			
			System.out.println("All transcript sites are read into a linked hashmap...");
			System.out.println(chrtss.size());
			
			}catch(Exception ex){
			ex.printStackTrace();
			System.out.println("There is problem reading the file.");
			}
			
		return(chrtss);
	}
	
	public static LinkedHashMap readInTSsArray(String t5file) {
		
		LinkedHashMap chrtssArray=new LinkedHashMap();
		
		/*Create whole TSs arrays with 22 chromosomes...*/
		LinkedHashMap chrLens=TScluster.readChrLen();
		Iterator itt=chrLens.keySet().iterator();
			
		try{
			
			String [] strands=new String [2] ;
			strands[0]="+";
			strands[1]="-";
			
			while(itt.hasNext()){
				String chr=itt.next().toString();
				
				
				int chrLength=Integer.parseInt((String) chrLens.get(chr));
				for(int i=0;i<strands.length;i++){
					String strand=strands[i];			
//					System.out.println(chrLens);
					int [] chrTSs=new int [chrLength];
					for(int ii=0;ii<chrLength;ii++){
						chrTSs[ii]=0;
					}
					
					String chrstrand=chr+"_"+strand;
					chrtssArray.put(chrstrand, chrTSs);
					System.out.println(chrstrand+"\t"+chrLength);
				}
				
			}
			
			
			BufferedReader br=new BufferedReader(new FileReader(t5file));
//			System.out.println(t5file);
			while(br.ready()){
				String str=br.readLine();
				String a[]=str.split("\t");
				
				String chr=a[0];
				String strand=a[1];
				int pos=Integer.parseInt(a[2]);
				int num=Integer.parseInt(a[3]);
				String chrstrand=chr+"_"+strand;
//				System.out.println(chrstrand);	
				
//				if(chrtssArray.keySet().contains(chr)){
//				if(chrtss.keySet().contains(chrstrand)){
							
				int [] chrTSs=(int []) chrtssArray.get(chrstrand);
				chrTSs[pos-1]=num;
//				}else{}
				
			}
			
			System.out.println("All transcript sites are read into a linked hashmap with TSs arrays...");
//			System.out.println(chrtssArray.size());s
			
			}catch(Exception ex){
			ex.printStackTrace();
			System.out.println("There is problem reading the file.");
			}
			
		return(chrtssArray);
	}
	
	public static void clusterChr(String chrstrand, LinkedHashMap tss){
		
		
		String bb[]=chrstrand.split("_");
    		String chr=bb[0];
//    		String strand=bb[1];
//    		
//    		int chrn=Integer.parseInt(chr.substring(3));
//    		int nstr=0;if(strand=="+"){nstr=1;}
    		
		LinkedHashMap chrLens=TScluster.readChrLen();
		int chrLength=Integer.parseInt((String) chrLens.get(chr));
		
		int [] chrPosLabs=new int [chrLength];
		for(int ii=0;ii<chrLength;ii++){
			chrPosLabs[ii]=0;
		}
		int clusterLab=1;
		
		Iterator it =tss.keySet().iterator();
        int clsNum=1;
		while(it.hasNext()){
        	
			LinkedHashMap tsclu=(LinkedHashMap) tsclus.get(chrstrand);
			
        		String chrPos=it.next().toString();
        		
        		String[] cc=chrPos.split("_");
//        		System.out.println(chrPos+"\t"+tss.get(chrPos));
        		int pos=Integer.parseInt(cc[2]);
        		int num= Integer.parseInt(tss.get(chrPos).toString());
        		
        		
        		if(chrPosLabs[pos-1]>0){
        		/* The position is already in previous clusters...*/	
        		}else{
        			/* The position is not in any cluster yet, put it in a new cluster or merge into old ones...*/
                		
        				int j=1;
                		boolean pmerged=false;
                		boolean uext=true;
                		boolean dext=true;
                		
                		int ncsta=0;
                		int ncend=(int) 1e10;
                		
                		mergeRes ress=new mergeRes(pmerged,uext,ncsta,dext,ncend,0){};
                		
        				while( j< (hfwinsize+1) && ! pmerged && (uext || dext) ){
        					
        					int ulable=chrPosLabs[pos-1-j];
        					int dlable=chrPosLabs[pos-1+j];
//        					System.out.println(ulable+"\t"+dlable);
        					
//        					int maxlab=Math.max(ulable, dlable);
        					
        					/*Only extend one end*/
        					if(uext && !dext){
        						/*Only upstream...*/
        						ncsta=pos-j;
        						
        						/*for debugging...*/
//            					if(chr.endsWith("chr1") && pos==76417) {System.out.println(ncsta+"\t"+j+"\t"+ulable+"\t"+dlable);}
            					
        						if(ulable>0){
        							 ress=TScluster.sMerge(chrstrand,"up", ulable, pos, pmerged, uext, ncsta, dext, ncend,0);
        						}else{
        							ress.ncsta=ncsta;
        						}
        						
        					}
        					if(dext && !uext){
        						ncend=pos+j;
        						/*Only downstream...*/
        						if(dlable>0){
        							 ress=TScluster.sMerge(chrstrand,"down", dlable, pos, pmerged, uext, ncsta, dext, ncend,0);
        						}else{
        							/*else do nothing; to deal with the previous generated ress, make the ncsta and ncend are correct....*/
        							ress.ncend=ncend;
        						}
        						
        					}
        					
        					/*for debugging...*/
//        					if(chr.endsWith("chr1") && pos==76417) {System.out.println(ncsta+"\t"+ress.ncsta+"\t"+ncend+"\t"+ress.ncend+"\t"+ress.pmerged+"\t"+uext+"\t"+dext);}
        					
        					
        					/*Extend both ends...*/
        					if(dext && uext){
        						ncsta=pos-j;
        						ncend=pos+j;
        						
	        					if(ulable==0){
	        						if(dlable==0){
	        							/* both ends are fine for this extension; 
	        							 * do nothing and go for further extension... */
	        							
	        						}else{
	        							/* extension overlaps with downstream cluster, merge them or not*/
	        							 ress=TScluster.sMerge(chrstrand,"down", dlable, pos, pmerged, uext, ncsta, dext, ncend,0);  
	        						}
	        						
	        					}else{
	        						
	        						if(dlable==0){
	        							/* extension overlaps with upstream cluster, merge them or not*
	        							 * ulable>0 & dlable==0*/ 
	        							 ress=TScluster.sMerge(chrstrand,"up", ulable, pos, pmerged, uext, ncsta, dext, ncend,0); 
	        							
	        						}else{
	        							/* extension overlaps with both up/downstream clusters, decide which side to merge or not....*/
	        							
	        							/*for debug..*/
//	        							if(chr.endsWith("chr1") && pos==76406) {System.out.println(ncsta+"\t"+ncend+"\t");}
	                					
	        							/* extend to one side or not...*/
	        							 ress=TScluster.bMerge(chrstrand,ulable, dlable, pos, pmerged, uext, ncsta, dext, ncend,0); 
	        							
	        						}
	        					}
        					}
        					
        					pmerged=ress.pmerged;
        					uext=ress.uext;
        					ncsta=Math.max(ncsta,ress.ncsta);
        					dext=ress.dext;
        					ncend=Math.min(ncend,ress.ncend);
//        					clusterLab=ress.clusterLab;
        					j++;
        				}
        				
        				if(ress.pmerged){
        					/*Add the cluster number code into the chromosome array....*/
        					
        					for(int k=ress.ncsta-1;k<ress.ncend;k++){
            					chrPosLabs[k]=ress.clusterLab;
            				}
            				if(ress.clusterLab==clusterLab){clusterLab++;}
        				}else{
        					/*for debug*/
//        					if(chr.endsWith("chr1") && pos==76417) {System.out.println(ncsta+"\t"+ress.ncsta+"\t"+ncend+"\t"+ress.ncend+"\t"+ress.pmerged);}
        					
        					/*create one new cluster...*/
        					int[] cluInfo=new int [4];
        					cluInfo[0]=Math.max(ncsta,ress.ncsta);
        					cluInfo[1]=Math.min(ncend, ress.ncend);
        					cluInfo[2]=pos;
        					cluInfo[3]=num;
//        					cluInfo[4]=chrn;
//        					cluInfo[5]=nstr;
//        					cluInfo[4]=chr;
        					
//        					System.out.println(cluInfo[0]+"\t"+cluInfo[1]+"\t"+cluInfo[2]+"\t"+cluInfo[3]+"\t"+clusterLab);
        					
        					tsclu.put(clusterLab,cluInfo);
        					tsclus.put(chrstrand, tsclu);
        					/*Add the cluster number code into the chromosome array....*/
            				for(int k=cluInfo[0]-1;k<cluInfo[1];k++){
//            					System.out.println(k+1+"\t"+clusterLab);
            					chrPosLabs[k]=clusterLab;
            				}
        					
        					clusterLab++;
        				}
        				
        				
        		}
        		
        }
		
        
//		return(clu);
		
	}
	
	public static mergeRes sMerge(String chrstrand,String updown,int clulable,int pos,boolean pmerged,boolean uext,int ncsta, boolean dext, int ncend,int clusterLab){
		mergeRes ress=new mergeRes(pmerged,uext,ncsta,dext,ncsta,clusterLab){};
		LinkedHashMap tsclu=(LinkedHashMap) tsclus.get(chrstrand);
		int[] cluInfo= (int[]) tsclu.get(clulable);
		int csta=cluInfo[0];
		int cend=cluInfo[1];
		int ctss=cluInfo[2];
		int cnum=cluInfo[3];
		
//		System.out.println(Math.abs(ctss-pos)+disthres);
		
		if(Math.abs(ctss-pos)<disthres){
			/*merge, update the upstream cluster infor...*/
			if(updown=="up"){
				cluInfo[1]=pos;
				ncend=pos;
			}else{if(updown=="down"){
				cluInfo[0]=pos;
				ncsta=pos;
			}else{
				System.out.println("The parameter for the sMerge is wrong....");
			}}
			tsclu.put(clulable,cluInfo);
			tsclus.put(chrstrand, tsclu);
//			System.out.println(updown+" is merged...");
			pmerged=true;
			clusterLab=clulable;
		}else{
			/*Not merge, update the start boundary here
			 * Depends on up or down, write out up or down....*/
			if(updown=="up"){
				uext=false;
				ncsta=ncsta+1;
			}else{if(updown=="down"){
				dext=false;
				ncend=ncend-1;
			}else{
				System.out.println("The parameter for the sMerge is wrong....");
			}}
			
		}
		
		ress.pmerged=pmerged;
		ress.uext=uext;
		ress.ncsta=ncsta;
		ress.dext=dext;
		ress.ncend=ncend;
		ress.clusterLab=clusterLab;
		
		
		return(ress);
	}
	
	public static mergeRes bMerge(String chrstrand, int ulable,int dlable,int pos,boolean pmerged,boolean uext,int ncsta, boolean dext, int ncend,int clusterLab){
		LinkedHashMap tsclu=(LinkedHashMap) tsclus.get(chrstrand);
//		System.out.println("bMerge is needed...");
		mergeRes ress=new mergeRes(pmerged,uext,ncsta,dext,ncend,clusterLab){};
		
		int[] ucluInfo= (int[]) tsclu.get(ulable);
		int ucsta=ucluInfo[0];
		int ucend=ucluInfo[1];
		int uctss=ucluInfo[2];
		int ucnum=ucluInfo[3];
		
		int[] dcluInfo= (int[]) tsclu.get(dlable);
		int dcsta=dcluInfo[0];
		int dcend=dcluInfo[1];
		int dctss=dcluInfo[2];
		int dcnum=ucluInfo[3];
		
		int dUp=Math.abs(uctss-pos);
		int dDown=Math.abs(dctss-pos);
		
		boolean mergeUp=false;
		boolean mergeDown=false;
		
		/*Could merge in both sides; select one side to merge*/
		if(dUp<disthres && dDown<disthres){
			
			/*Closer to up...*/
			if(dUp<dDown){mergeUp=true;}
			
			/*Same distance to up and down; bigger number...*/
			if(dUp==dDown){
				if(ucnum>=dcnum){
					mergeUp=true;	
				}
				else{mergeDown=true;}}
			
			/*Closer to down...*/
			if(dUp>dDown){mergeDown=true;}
			
		}
		
		/*Merge to upstream cluster...*/
		if(dUp<disthres && dDown>=disthres){
			mergeUp=true;
		}
		
		/*Merge to downstream cluster...*/
		if(dUp>=disthres && dDown<disthres){
			mergeDown=true;
		}
		
		if(mergeUp){
			ucluInfo[1]=pos;
			ncend=pos;
			tsclu.put(ulable,ucluInfo);
//			System.out.println(updown+" is merged...");
			pmerged=true;
			clusterLab=ulable;
			tsclus.put(chrstrand, tsclu);
		}
		
		if(mergeDown){
//			ress=TScluster.sMerge("down", dlable, pos, pmerged, uext, ncsta, dext, ncend,0);
			dcluInfo[0]=pos;
			ncsta=pos;
			tsclu.put(dlable,dcluInfo);
//			System.out.println(updown+" is merged...");
			pmerged=true;
			clusterLab=dlable;
			tsclus.put(chrstrand, tsclu);
		}
		
		if(!pmerged){
			ncsta=ncsta+1;
			ncend=ncend-1;
			uext=false;
			dext=false;
		}
		
		ress.pmerged=pmerged;
		ress.uext=uext;
		ress.ncsta=ncsta;
		ress.dext=dext;
		ress.ncend=ncend;
		ress.clusterLab=clusterLab;
		
		
		return(ress);
	}
	
	
	public static HashMap tMerge(LinkedHashMap clu, int clulable,int pos){
		HashMap megres=new HashMap(){};
		boolean pmerged=false;
		
		
		
		
		return(megres);
		
	}
	public static int countArray( int[] chrPosLabs,int start,int end){
		/* start and end here are genome positions, with +1 ...
		 * */
		
		int sums=0;
		for(int i =start-1; i < end;  i++){
			sums=chrPosLabs[i]+sums;
		}
		return(sums);
	}
	
	 public static LinkedHashMap readChrLen(){
	    	
		 LinkedHashMap chrLens=new LinkedHashMap();
	    	
	    	try { 
	    	
	    		BufferedReader br = new BufferedReader(new FileReader(chrLenFile));
	    	
				while(br.ready()){
					
					String str=br.readLine();
					String a[]=str.split("\t");
					
					chrLens.put(a[0], a[1]);
					
				}
				
	    	} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
					System.out.print(e);
				}
				
	    	return(chrLens);
	    	
	 }


	
	
	

}
