package evolleadmod2;

import static java.lang.Math.exp;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.stat.descriptive.moment.Skewness;

public class Utility {
    private Random randomGenerator;
    private Skewness skewnessCalculator;
    
    public Utility(Random pRandomGenerator){
        this.randomGenerator = pRandomGenerator;
        this.skewnessCalculator = new Skewness();
    }
    
    
    public double randomDouble(){
        return randomGenerator.nextDouble();
    }
    // Probability test
    public int testProb(double pMu){
        if(Math.random() >= pMu){
            return 0;                                                           //Nothing
        }
        else{
           return 1;                                                            //Event        
        }
    }
    //Writer
    public List writeFile(List pList){
        
        return(pList);
    }
    
    //Sample an individual except one
    public int randomSampleOther(int pArrayL, int pIndex){
        //if(pLArray == 1){System.out.println("BUG randomSampleOther List too short");}
        ArrayList<Integer> resList = new ArrayList();
        for(int i=0; i<pArrayL; i++){resList.add(i);}
        resList.remove(pIndex);
        int index = (int)(Math.floor(Math.random()*(pArrayL-1)));
        return resList.get(index);
        }
    
    //Sample a number of individuals except one
    public int[] randomSampleOtherList(int pArrayL, int pSampleSize, int pIndex){
        //if(pLArray == 1){System.out.println("BUG randomSampleOther List too short");}
        ArrayList<Integer> resList = new ArrayList();
        for(int i=0; i<pArrayL; i++){resList.add(i);}
        resList.remove(pIndex);
        int indexTemp;
        int[] indexL = new int[pSampleSize];
        for(int i=0; i<pSampleSize; i++){
            indexTemp = (int)(Math.floor(Math.random()*(resList.size())));
            indexL[i]= resList.get(indexTemp);
            resList.remove(indexTemp);
        }
        return indexL;
        }
    
    //Binary search algorithm for sampling probability
    public int probSample(double[] pArrayProb, double pKey){             
                            int lower = 0;
                            int upper = pArrayProb.length-1;
                            int mid;
                            while (lower < upper){ 
                                mid = (int)Math.floor((lower + upper )/2);      
                                if((pArrayProb[mid] - pKey) > 0){
                                    upper = mid;
                                }
                                else{
                                    lower = mid + 1;
                                }
                            }
                            return lower;
                        }
    
    //Calcul of variance
    public double sdPref(List<Individual> pList){
        double M= 0;
        double S = 0;
        double oldM;
        for(int i = 0; i < pList.size(); i++){
            oldM = M;
            M += (pList.get(i).getXNego()-M)/(i+1);
            S += (pList.get(i).getXNego()-M) * (pList.get(i).getXNego()-oldM);
        }
        return (Math.sqrt(S/(pList.size()-1)));
    }
    
    //Calcul of mean
    public double mean(double[] m) {
        double sum = 0;
        for (int i = 0; i < m.length; i++) {
            sum += m[i];
        }
        return sum / m.length;
    }
    
    //Calcul of mean of the preferences
    public double meanPref(List<Individual> pList){
        double sumX = 0;
        for(int i = 0; i < pList.size(); i++){
            sumX += pList.get(i).getXNego();
        }
        return(sumX/(pList.size()));
    }
    
    // Method to create a discrete distribution of negative exponential form
    public double[] functionDistribution(int pSize, double pFConsensus){
        double[] valueDistribution2 = new double[pSize];
        double sumDistrib = 0.0d;
        for(int i=0; i < pSize; i++){
            sumDistrib = sumDistrib + exp(pFConsensus * (-i));                  // Sum of the whole share to normalize the final values
        }
        for(int i=0; i < pSize; i++){
            valueDistribution2[i] = exp(pFConsensus * (-i))/sumDistrib;
        }
        return(valueDistribution2);
    }
    
        //Method to mutate a value (between 0 and 1) by normal distribution
    public double mutation(double pMean,double pSigma,double pMinValue,double pMaxValue){
        NormalDistribution mutDis = new NormalDistribution(pMean,pSigma);
        double res = mutDis.sample();
        if(res < pMinValue){res = pMinValue;}
        if(res > pMaxValue){res = pMaxValue;}
        return(res);
    }
    
    public double skewnessAlpha(List<Individual> pList) {
    	double[] alphaList = new double[pList.size()];
    	for(int i =0; i<pList.size(); i++) {alphaList[i]=pList.get(i).getAlpha();}
    	double res = this.skewnessCalculator.evaluate(alphaList);
    	return(res);
    }
    
    public double meanSkewnessAlpha(List<List<Individual>> pList) {
    	double res = 0;
    	for (int iPatch = 0; iPatch<pList.size();iPatch++) {
    		double[] alphaList = new double[pList.get(iPatch).size()];
    		for(int i =0; i<pList.get(iPatch).size(); i++) {alphaList[i]=pList.get(iPatch).get(i).getAlpha();}
    		res = res + this.skewnessCalculator.evaluate(alphaList)/pList.size();
    	}
    	return(res);
    }
    
    public double skewnessAlphaGlobal(List<List<Individual>> pList) {
    	int totalSize = 0;
    	
    	for (int i = 0; i<pList.size();i++) {totalSize+=pList.get(i).size();}
    	double[] alphaList = new double[totalSize];
    	
    	int k=0;
    	for(int i =0; i<pList.size(); i++) {
    		for(int j=0; j<pList.get(i).size();j++) {
    		alphaList[k]=pList.get(i).get(j).getAlpha();
    		k++;
    		}
    	}
    	double res = this.skewnessCalculator.evaluate(alphaList);
    	return(res);
    }
    
}
