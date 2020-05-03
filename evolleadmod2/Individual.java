package evolleadmod2;

import java.util.Comparator;

public class Individual implements Comparable<Individual> {
    // Influence
    private double alpha;
    // Initial opinion/preference
    private double x;
    // Preference Negotiated
    private double xNego;
    //Number of negotiations
    private double countNego;
    // Migration preference (if conditional migration: not used here)
    private double z;
    //fitness
    private double w;
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
        x = 0;
        xNego = 0;
        countNego = 0;
        z = 0;
        w = -1;
        //m = -1;
        off = -1;
        nbreIndividual++;
        bias = 0;
    }
    
    public Individual(double pAlpha, double pX, double pZ){
        alpha = pAlpha;
        x = pX; 
        xNego = pX;
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
        x = pIndividual.x;
        xNego = pIndividual.xNego;
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
    public void setX(double pX){x = pX;}
    public void setXNego(double pXNego){xNego = pXNego;}
    public void setCountNego(double pCountNego){countNego = pCountNego;}
    public void setZ(double pZ){z = pZ;}
    public void setW(double pW){w = pW;}
    public void setOff(int pOff){off = pOff;}
    public void setBias(double pBias) {bias=pBias;}
    //Getters
    public double getAlpha(){return alpha;}
    public double getX(){return x;}
    public double getXNego(){return xNego;}
    public double getCountNego(){return countNego;}
    public double getZ(){return z;}
    public double getW(){return w;}
    public int getOff(){return off;}
    public double getBias() {return bias;}
    
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


