package evolleadmod2;
import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.ZipOutputStream;

import com.beust.jcommander.Parameter;

public class Parameters {
	  @Parameter
	  public List<String> parameters = new ArrayList<>();
	  
	  //Simulation
	  @Parameter(names = {"-ID"}, description = "ID of process on cluster")
	  public Long processID = 0L;
	  @Parameter(names = {"-sd"}, description = "step for printing data in gen")
	  public Integer stepData;
	  @Parameter(names = {"-de"}, description = "print every individual or patch")
	  public Integer detail;
	  @Parameter(names = {"-fp"}, description = "first generation at which start to print")
	  public Integer firstPrint;
	  
	  //Evolution processes
	  @Parameter(names = {"-S"}, description = "Number of simulations")
	  public Integer nSimul;
	  @Parameter(names = {"--nGen","-G" }, description = "Number of generations")
	  public Integer nGen;
	  @Parameter(names = {"--sizePatchIni","-N" }, description = "Initial size of population in one patch")
	  public Integer nIndIni;
	  @Parameter(names = {"--mutationRate","-mu" }, description = "Mutation rate")
	  public Double mu;
	  @Parameter(names = {"--strengthMutation","-sigma" }, description = "strength of the mutations")
	  public Double sigma;
	  
	  //Ecological parameters
	  @Parameter(names = {"--nPatch","-P" }, description = "Number of patches")
	  public Integer P;
	  @Parameter(names = {"--carryingCapacity","-K" }, description = "Carrying capacity")
	  public Integer K;
	  
	  //Life history traits
	  @Parameter(names = {"--growthRateIntrinsic","-rI" }, description = "Intrinsic growth rate of one individual")
	  public Double rI;
	  @Parameter(names = {"--growthRateMax","-rMax" }, description = "Maximum additional growth rate")
	  public Double rBMax;
	  @Parameter(names = {"--growthRateIncrease","-rInc" }, description = "Steepness of the increase of growth rate")
	  public Double rBInc;
	  @Parameter(names = {"--migrationRate","-m" }, description = "Migration rate")
	  public Double m;
	  
	  //Collective action and additional resources
	  @Parameter(names = {"--collectiveRessourcesMax","-bMax" }, description = "Maximum ressources produced by collective action")
	  public Double bMax;
	  @Parameter(names = {"--collectiveRessourcesMid","-bMid" }, description = "Population size which produced half of the maximum possible collective ressources")
	  public Double bMid;
	  @Parameter(names = {"--collectiveRessourcesInc","-bInc" }, description = "Steepness of the increase of collective ressources")
	  public Double bInc;
	  @Parameter(names = {"--lambda","-L" }, description = "Amount of ressources remaining at the next generation")
	  public Double lambda;

	  //Collective decision-making process
	  @Parameter(names = {"--nListeners","-nL" }, description = "Number of listeners during a single negotiation event")
	  public Integer nL;
	  @Parameter(names = {"--consensusThreshold","-xThr" }, description = "Variance threshold under which consensus is considered reached")
	  public Double xThr;
	  @Parameter(names = {"-k" }, description = "Parameter controlling how much alpha affects talkativeness")
	  public Double kAlpha;
	  
	  //Collective institutions
	  @Parameter(names = {"--costCollectiveDecision","-C" }, description = "Inverse of Total cost of collective decision making")
	  public Double C;
	  @Parameter(names = {"--distribution","-d" }, description = "Skewness of distribution of ressources")
	  public Double d;
	  

	  
	  


      
}
