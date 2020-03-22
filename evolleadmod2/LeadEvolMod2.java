package evolleadmod2;


import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.lang.management.ManagementFactory;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.zip.*;
import org.apache.commons.io.FileUtils;
import org.apache.commons.math3.distribution.*;
import org.apache.commons.math3.stat.descriptive.moment.Variance;
import org.apache.commons.cli.*;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;

import sun.security.util.Length;


public class LeadEvolMod2 {
            
    public static void main(String[] args) throws FileNotFoundException, IOException, InterruptedException {
    	
    	//Initialize and parse parameters from options
    	Parameters parameters = new Parameters();
    	JCommander.newBuilder()
    	  .addObject(parameters)
    	  .build()
    	  .parse(args);
        
        //Working directory
        String wd;
    	// wd = "C:/Users/40011091/Phd-Thèse/A1-Projects/EvolLeadMod/C-Results/Res 4.0/";
        wd = System.getProperty("user.dir")+"/";

        //Name file   
        String nameFile= "";
        for(int i=4; i<args.length+1; i=i+2) {nameFile = nameFile + args[i-1] + "_" + args[i];}
        
        //nameFile = nameFile.substring(1) + "-";
        
        //Seed initialization
        long processID;
        long seed;
        
        System.out.println(nameFile);
        //Writer initialization
        //Might have to change to 0 and 1 when running on personal machine because java ignore first argument
        File f = new File(wd + args[1] + args[2] + nameFile + ".zip");
        ZipOutputStream out = new ZipOutputStream(new FileOutputStream(f));
        
        
        for(int iSimul =0; iSimul < parameters.nSimul; iSimul++){
        new Model(
                parameters.nGen, parameters.P, parameters.nIndIni, 
                parameters.rI, parameters.rBMax, parameters.rBInc,
                parameters.K, parameters.bMax, parameters.bInc, parameters.bMid,
                parameters.C, parameters.lambda, parameters.costNego,
                parameters.fThr, parameters.nL,
                parameters.mu, parameters.sigma, parameters.d,
                parameters.costM, parameters.m,
                parameters.theta,
                out, seed =  parameters.processID + iSimul, nameFile,
                parameters.stepData, parameters.detail, parameters.firstPrint
            ).run(); 
        };
         
        out.close();
     

    }
    
}
