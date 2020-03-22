package evolleadmod2;

import java.util.Comparator;

public class Individual implements Comparable<Individual> {
    // Communication skill
    private double alpha;
    //Talkativeness (Probability of speaking)
    private double theta;
    // Preference (fairness)
    private double f;
    // Preference Negotiated
    private double fNego;
    //Number of negotiations
    private double countNego;
    // Migration preference
    private double z;
    //fitness
    private double w;
    //Conditional migration ?
    //private int m;
    //Number of offsprings
    private int off;
    //Bias
    private double bias;
    
    public static int nbreIndividual = 0;
    public static int getNbreIndividual(){
        return nbreIndividual;
    }
    
    public Individual(){
        System.out.println("Creation of a default individual");
        alpha = 0;
        theta = 1;
        f = 0;
        fNego = 0;
        countNego = 0;
        z = 0;
        w = -1;
        //m = -1;
        off = -1;
        nbreIndividual++;
        bias = 0;
    }
    
    public Individual(double pAlpha, double pTheta, double pF, double pZ){
        alpha = pAlpha;
        theta = pAlpha;
        f = pF; 
        fNego = pF;
        countNego = 0;
        z = pZ;
        w = -1;
        //m = 0;
        off = -1;
        nbreIndividual++;
        bias = 0;
    }
    
    public Individual(Individual pIndividual){
        alpha = pIndividual.alpha;
        theta = pIndividual.theta;
        f = pIndividual.f;
        fNego = pIndividual.fNego;
        countNego = 0;
        z = pIndividual.z;
        w = -1;
        //m = 0;
        off = -1;
        nbreIndividual++;
        bias = 0;
    }
    //Setters
    public void setAlpha(double pAlpha){alpha = pAlpha;}
    public void setTheta(double pTheta){theta = pTheta;}
    public void setF(double pF){f = pF;}
    public void setFNego(double pFNego){fNego = pFNego;}
    public void setCountNego(double pCountNego){countNego = pCountNego;}
    public void setZ(double pZ){z = pZ;}
    public void setW(double pW){w = pW;}
    //public void setM(int pM){m = pM;}
    public void setOff(int pOff){off = pOff;}
    public void setBias(double pBias) {bias=pBias;}
    //Getters
    public double getAlpha(){return alpha;}
    public double getTheta(){return theta;}
    public double getF(){return f;}
    public double getFNego(){return fNego;}
    public double getCountNego(){return countNego;}
    public double getZ(){return z;}
    public double getW(){return w;}
    //public int getM(){return m;}
    public int getOff(){return off;}
    public double getBias() {return bias;}
    
    public String Description(){
        return "\nalpha = " + alpha + " \ntheta = " + theta + "\nw = " + w + "\nf = " + f + "\nfNego" + fNego + "\nz = " + z + "\noff = " + off ;
    }
    
    public void addCountNego(){
        this.countNego++;
    }
    
    @Override
    public int compareTo(Individual obj)
    {
        // compareTo returns a negative number if this is less than obj, 
        // a positive number if this is greater than obj, 
        // and 0 if they are equal.
        if(this.alpha < obj.alpha)
          return 1;
        else if(obj.alpha < this.alpha)
          return -1;
          return 0;
    }
    

	
    
   
}


