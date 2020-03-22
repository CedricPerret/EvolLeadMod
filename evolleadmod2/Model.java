package evolleadmod2;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.io.FileUtils;
import org.apache.commons.math3.distribution.*;
import java.util.zip.*;
import org.apache.commons.math3.stat.descriptive.moment.Variance;
import org.apache.commons.math3.stat.descriptive.moment.Skewness;;

public class Model extends Thread {
    private int nGen;private int P;private int nIndIni;
    private double rI; private double rBMax; private double rBInc;
    private int K; private double bMax; private double bInc; private double bMid;
    private double C; private double lambda; private double costNego;
    private double fThr; private int nL;
    private double mu; private double sigma; private double d; 
    private double costM;; private double m;
    private boolean theta;
    private ZipOutputStream pw;
    //private PrintWriter pw;
    private long seed;
    private String nameFile;
    private int stepData;
    private int detail;
    private int firstPrint;
    
    public Model(
        int nGen,int P,int nIndIni,
        double rI, double rBMax, double rBInc,
        int K, double bMax, double bInc, double bMid,
        double C, double lambda, double costNego,
        double fThr, int nL,
        double mu, double sigma, double d, 
        double costM, double m,
        boolean theta,
        ZipOutputStream pw,
        //PrintWriter pw,
        long seed,
        String nameFile,
        int stepData,
        int detail,
        int firstPrint
    )
    {
        this.nGen = nGen;this.P = P;this.nIndIni = nIndIni;
        this.rI = rI; this.rBMax = rBMax; this.rBInc = rBInc;
        this.K = K; this.bMax = bMax; this.bInc = bInc; this.bMid = bMid;
        this.C = C; this.lambda = lambda; this.costNego = costNego;
        this.fThr = fThr; this.nL = nL;
        this.mu = mu; this.sigma = sigma; this.d = d; 
        this.costM = costM; this.m = m;
        this.theta = theta;
        this.pw = pw;
        this.seed = seed;
        this.nameFile = nameFile;
        this.stepData = stepData;
        this.detail = detail;
        this.firstPrint = firstPrint;
    }
    
    //To know when every models are finished
    public static int nbreModelOver = 0;
        public static int getNbreModelOver(){
        return nbreModelOver;
        }
    
    //@Override
    public void run(){
        double fVar; double thetaSum; double probTemp;int speaker; int nL2; double diffAlpha; int nEvent;double bTot; double r; double I; double w; double rB;
        int off; double alphaTemp; double thetaTemp; double zTemp; double fTemp; String resTemp =""; int printJump = 10;
        double[][] resFConsensus = new double[nGen+1][P];
        double[][] resTConsensus = new double[nGen+1][P];
        double[][] resBTot = new double[nGen+1][P];
        double[] bTotPre = new double[P];
        
        //Random generator
        Utility utility = new Utility(new Random(seed));
        
        //Create Entry inside the zip
        System.out.println(nameFile);
        
        ZipEntry e = new ZipEntry(nameFile + "-seed_ " + seed + ".txt");
        try {
            pw.putNextEntry(e);
        } catch (IOException ex) {
            Logger.getLogger(Model.class.getName()).log(Level.SEVERE, null, ex);
        }
        
            //Writing parameters
            
//            try {
//                pw.write(("nGen,P,nIndIni,rI,K,bMax,bInc,C,fIni,fThr,nL,mu,sigma,d,costM,m,lambda,\r\n"+
//                        nGen + P + nIndIni + rI + K + bMax + bInc + C + fIni +fThr +nL + mu +sigma + d + costM + m + lambda + "\r\n"+
//                        seed+"\r\n").getBytes()); 
//            } catch (IOException ex) {
//                Logger.getLogger(Model.class.getName()).log(Level.SEVERE, null, ex);
//            }
//        
//            if(seed == 0){pw.write(("nGen,P,nIndIni,rI,K,bMax,bInc,bMid,C,fIni,fThr,nL,mu,sigma,d,costM,m,lambda,\r\n"+
//                nGen + P + nIndIni + rI + K + bMax + bInc + bMid + C + fIni +fThr +nL + mu +sigma + d + costM + m + lambda + "\r\n"+
//                seed+"\r\n"));}
           
           //Initialization ====================================================
           //Table for f* and t*
            resFConsensus = new double[nGen+1][P];                             
            resTConsensus = new double[nGen+1][P];
            resBTot = new double[nGen+1][P];
            //Table for populations
            List<List<Individual>> popNow = new ArrayList<>();
            List<List<Individual>> popNext = new ArrayList<>();
            for(int j=0; j<P ;j++){
                    popNow.add(new ArrayList<>());
                    popNext.add(new ArrayList<>());
                    bTotPre[j] = 0;                                             //Initial bTotPre Array
            //Initial population t=0
                    for(int k=0; k<nIndIni; k++){     
                        popNow.get(j).add(new Individual(utility.randomDouble(),utility.randomDouble(),utility.randomDouble(),utility.randomDouble()));
                    }
                }
            // SIMULATIONS =====================================================
            for(int i=0; i<nGen; i++){
                if((i%Math.round(nGen/printJump)) == 0){System.out.println("Simul " + seed + " Gen " + (i*100)/nGen + "%");}  
                
                for(int j=0; j<P; j++){
                    //Empty patch-----------------------------------------------
                    //If we need to write but patch is empty, we write a list of NA
                    if(i>=firstPrint&(i==0 || (i+1)%stepData==0)) {
                        if(popNow.get(j).isEmpty())
                        {resTemp = resTemp + "NA,NA,NA,NA,NA"
                                //+ ",NA"
                                + ",NA,NA,NA," +j+","+ i+","+ seed+ "\r\n";
                        }
                    }
                    //Jump the rest
                    if(popNow.get(j).isEmpty()){continue;}        
                    
               
                    //NEGOTIATION PROCESS=======================================
                    int maxNegoEvent = 1;

                    double[] probSpeaker = new double[popNow.get(j).size()];    // Initial array of the probability of speaking of each individual
                    //Talkativeness Theta : probability of speaking-------------
                    thetaSum = 0;                                               //Sum of theta to then normalize
                    for(int k=0; k<popNow.get(j).size(); k++){                       
                    thetaSum += Math.pow(popNow.get(j).get(k).getTheta(),4);
                    }
                    probTemp = 0;                
                    for(int k=0; k<popNow.get(j).size(); k++){                  //We calculate the weigthed probabilities
                    probSpeaker[k]= probTemp + (Math.pow(popNow.get(j).get(k).getTheta(),4)/thetaSum);
                    probTemp = probSpeaker[k];
                    }
                    
                    //List of listeners----------------------------------------
                    int[] listenerList;
                    if(popNow.get(j).size() <= nL){nL2 = (popNow.get(j).size())-1;}     //To avoid looking for non existing listeners   
                    else {nL2 = nL;}
                    
                    for(int negoCounter = 0; negoCounter<maxNegoEvent; negoCounter++) {
	                    fVar = utility.variancePref(popNow.get(j));                 // Initial pref sd
	                    nEvent = 0;                                                 // Initial number of negotiation event
	                    //Simulations of negotiations-------------------------------
	                    while(fVar > fThr){
	                        speaker = utility.probSample(probSpeaker, utility.randomDouble());           //Speaker chosen in function of theta
	                        popNow.get(j).get(speaker).addCountNego();                              //Count number of negociation individual take part in
	                        listenerList = utility.randomSampleOtherList(popNow.get(j).size(), nL2,speaker);            //Listener chosen randomly in the rest of the population
	                        for(int l=0; l<nL2; l++){                               //Update the preference of each listener
	                            diffAlpha =  popNow.get(j).get(speaker).getAlpha() - popNow.get(j).get(listenerList[l]).getAlpha();
	                            if(diffAlpha <= 0.01){diffAlpha = 0.01;}             // If values of alpha are the same
	                            popNow.get(j).get(listenerList[l]).setFNego(popNow.get(j).get(listenerList[l]).getFNego() + diffAlpha * (popNow.get(j).get(speaker).getFNego()-popNow.get(j).get(listenerList[l]).getFNego()));
	                        }
	                        nEvent++;     
	                        fVar = utility.variancePref(popNow.get(j));		
	                        //if(nEvent == C){break;}                                 //To faster the code : If cost of consensus > 1, we stop
	                    }
	                    
	                    //Outcomes of the negotiation process----------------------- 
	                    resTConsensus[i][j] = resTConsensus[i][j]+(nEvent/maxNegoEvent);           //We write the tConsensus
	                    resFConsensus[i][j] = utility.meanPref(popNow.get(j));      //We write the fConsensus
	                    for(Individual ind : popNow.get(j)) {
	                    	ind.setBias(ind.getBias()+
	                    			(1-Math.abs(ind.getF()-resFConsensus[i][j])));
	                    	ind.setF(utility.randomDouble());
	                    	ind.setFNego(ind.getF());
	                    }
                    }

                    //COLLECTIVE ACTION=========================================
                    // Cost of institution proportional to t*
                    I = resTConsensus[i][j]*C;
                    //Total benefit of the collective action
                    	//With additive cost
                    		//Sigmoid
                    bTot = ((bMax/(1+Math.exp(-bInc*(popNow.get(j).size() - bMid))))-I)*1;if(bTot < 0){bTot = 0;}
                    		//Logistic
                    //bTot = (bMax*(1-Math.exp(-bInc*(popNow.get(j).size()))))-I;if(bTot < 0){bTot = 0;}
                    
                    if((bMax/(1+Math.exp(-bInc*(100 - bMid))))>(bMax/(1+Math.exp(-bInc*(200 - bMid))))) {
                    	System.out.println("Probleme with bTot function");
                    }
                    resBTot[i][j] = bTot; 
                    	
                    //DISTRIBUTION OF RESOURCES=================================
                    //Sorting individuals---------------------------------------
                    Collections.sort(popNow.get(j));                            //We sort individuals based on alpha
                    double distribSum = 0;
                    double[] valueDistribution = new double[popNow.get(j).size()]; 
                    
	                for(int k=0; k<popNow.get(j).size(); k++){                       
	                  distribSum += 1+d*popNow.get(j).get(k).getBias();
	                }               
	                for(int k=0; k<popNow.get(j).size(); k++){                  //We calculate the weigthed probabilities
	                  valueDistribution[k]= (1+d*popNow.get(j).get(k).getBias())/distribSum;
	                }
                    
                    //For distribution proportional to alpha
//                    for(int k=0; k<popNow.get(j).size(); k++){                       
//                    distribSum += (1+d*popNow.get(j).get(k).getAlpha());
//                    }               
//                    for(int k=0; k<popNow.get(j).size(); k++){                  //We calculate the weigthed probabilities
//                    valueDistribution[k]= ((1+d*popNow.get(j).get(k).getAlpha())/distribSum);
//                    }
//                    
                    //If distribution determined by group decision
                    //double[] valueDistribution = functionDistribution(popNow.get(j).size(),resFConsensus[i][j]);
	                
                    //Calcul of personal share
                    for(int k=0; k<popNow.get(j).size(); k++){
                        rB = rBMax*(1-Math.exp(-rBInc*(bTot+(bTotPre[j]*lambda)) * valueDistribution[k])) - (popNow.get(j).get(k).getCountNego()*costNego);      //Calcul increase of growth rate in function of collective action ressource at t and t-1  
                    	if(rB<0){rB=0;}
                        w = (rI/(1+(popNow.get(j).size()/K)))  + rB ;                     //Calcul fitness
                        popNow.get(j).get(k).setW(w);
                        PoissonDistribution offDistrib = new PoissonDistribution(popNow.get(j).get(k).getW());
                        off = offDistrib.sample();                                  // Number of offsprings sampled from a poisson distribution
                        popNow.get(j).get(k).setOff(off);

                        //WRITING-----------------------------------------------
                        if(detail==2) {
	                        if(i>=firstPrint&(i==0 || (i+1)%stepData==0)) {
	                            resTemp = resTemp + popNow.get(j).get(k).getAlpha()
	                            +","
	                            + popNow.get(j).get(k).getBias()
	                            +","
	                            + popNow.get(j).get(k).getF()
	                            +","
	                            + popNow.get(j).get(k).getW()
	                            //+","
	                            //+ popNow.get(j).get(k).getOff()
	                            +","
	                            + resFConsensus[i][j]
	                            + ","
	                            + resTConsensus[i][j]
	                            + ","
	                            + resBTot[i][j]
	                            +","
	                            + j
	                            +","
	                            + i
	                            +","
	                            + seed     
	                            + "\r\n";
	                        }
                        }
                    }
                    if(detail==1) {
                        if(i>=firstPrint&(i==0 || (i+1)%stepData==0)) {
                            resTemp = resTemp + utility.skewnessAlphaGlobal(popNow)
                            + ","
                            + utility.skewnessAlpha(popNow.get(j))
                            +","
                            + resTConsensus[i][j]
                            + ","
                            + resBTot[i][j]
                            + ","
                            + resFConsensus[i][j]
                            + ","
                            + popNow.get(j).size()
                            + ","
                            + j
                            +","
                            + i
                            +","
                            + seed     
                            + "\r\n";
                        }
                    }
                    
                    if(detail==2) {
	                    if(i>=firstPrint&(i==0 || (i+1)%stepData==0)) {try {
	                        pw.write(resTemp.getBytes());
	                                                pw.flush();
	                        } catch (IOException ex) {
	                            Logger.getLogger(Model.class.getName()).log(Level.SEVERE, null, ex);
	                        };resTemp = "";}
                    }
                    
                    bTotPre[j] = bTot;
                }
                //Print at each generation
                if(detail==0) {
                	double[] sizePatch = new double[P];
                	for(int j = 0; j<P;j++) {sizePatch[j]=(double)popNow.get(j).size();}
                	if(i>=firstPrint&(i==0 || (i+1)%stepData==0)) {
                		resTemp = resTemp + utility.skewnessAlphaGlobal(popNow)
                		+ ","
                		+ utility.meanSkewnessAlpha(popNow)
                		+ ","
                		+ utility.mean(resTConsensus[i])
                		+ ","
                		+ utility.mean(resBTot[i])
                		+ ","
                		+ utility.mean(resFConsensus[i])
                		+","
                		+ utility.mean(sizePatch)
                		+","
                		+ i
                		+","
                		+ seed     
                		+ "\r\n";
                	}
                }
                //Write at each generation if detail level is low enough
                if(detail!=2) {
                	if(i>=firstPrint&(i==0 || (i+1)%stepData==0)) {try {
                		pw.write(resTemp.getBytes());
                		pw.flush();
                	} catch (IOException ex) {
                		Logger.getLogger(Model.class.getName()).log(Level.SEVERE, null, ex);
                	};resTemp = "";}
                }
                //REPRODUCTION==================================================
                if(i == nGen){break;}                                           // No offspring for last generation
                for(int j=0; j<P; j++){
                    if(popNow.get(j).isEmpty()){continue;}                       //Empty patch (break different than continue)
                    for(int k=0; k<popNow.get(j).size(); k++){
                        alphaTemp = popNow.get(j).get(k).getAlpha();
                        thetaTemp = popNow.get(j).get(k).getTheta();
                        fTemp = utility.randomDouble();
                        zTemp = popNow.get(j).get(k).getZ(); 
                        for(int l=0; l<popNow.get(j).get(k).getOff(); l++){
                            //Mutations-----------------------------------------
                            if(utility.testProb(mu) == 1){alphaTemp = utility.mutation(alphaTemp,sigma,0.1,1);} 
                            if(utility.testProb(mu) == 1){zTemp = utility.mutation(zTemp,(sigma),0.0,1);}
                            if(theta == true){thetaTemp = alphaTemp;}
                            else{thetaTemp = 1;}
                            //Migration-----------------------------------------
                            if(P != 1 && utility.testProb(m) == 1){       
                                if(utility.testProb(costM)==0){     
                                    popNext.get(utility.randomSampleOther(P,j)).add(new Individual(alphaTemp,thetaTemp,fTemp,zTemp));}
                                }
                            else{popNext.get(j).add(new Individual(alphaTemp,thetaTemp,fTemp,zTemp));}
                        }
                    }
                }
                
                //We translate popNext to popNow
                popNow.clear();                             
                for(int j=0; j<P; j++){
                    popNow.add(new ArrayList<>(popNext.get(j)));
                    popNext.get(j).clear();
                }
                    
                }
                System.out.println("nGen,P,nIndIni,rI,K,bMax,bInc,bMid,C,fThr,nL,mu,sigma,d,costM,m,lambda,\r\n"+
                nGen + P + nIndIni + rI + K + bMax + bInc + bMid + C +fThr +nL + mu +sigma + d + costM + m + lambda + "\r\n"+
                seed+"\r\n");
                
               nbreModelOver++;
        try {
            pw.closeEntry();
        } catch (IOException ex) {
            Logger.getLogger(Model.class.getName()).log(Level.SEVERE, null, ex);
        }
              
      
        }
    }

