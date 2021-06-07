#@ MoleculeArchive archive

import de.mpg.biochem.mars.molecule.*
import de.mpg.biochem.mars.table.*
import de.mpg.biochem.mars.util.*
import org.scijava.table.*
import org.apache.commons.math3.stat.regression.SimpleRegression


def regression = new SimpleRegression(true)
def E_list = []
def S_list = []

//create a list with all data points to be used in the fit
archive.getMoleculeUIDs().stream().forEach({UID ->
  def molecule = archive.get(UID)
  if (molecule.hasTag("Active_single")){
  	double E = molecule.getParameter("iiiEapp")
    double S = 1 / molecule.getParameter("iiiSapp")
    S_list.add(S)
    E_list.add(E)}})

//add the datapoints to be included in the regression analysis and store found parameters
for (i in [E_list,S_list].transpose()){
	regression.addData(i[0],i[1])
}

double a = regression.getSlope()
double b = regression.getIntercept()
double beta = a + b - 1
double gamma = (b - 1) / (a + b - 1)

archive.metadata().forEach{metadata -> 
	metadata.setParameter("beta", beta)
	metadata.setParameter("gamma", gamma)
}

//calculate F_AA and F_DD values for each molecule and calculate the E and S values

archive.getMoleculeUIDs().stream().forEach({UID ->
  def molecule = archive.get(UID)
  if (molecule.hasTag("Active_single")){
  	double FAD = molecule.getParameter("F_AD")
  	double Idemdex = molecule.getParameter("iiIdemdex")
  	double Iaemaex = molecule.getParameter("iiIaemaex")
  	double FDD = gamma * Idemdex
  	double FAA = Iaemaex / beta
    molecule.setParameter("F_DD",FDD)
    molecule.setParameter("F_AA",FAA)
    double E = FAD / (FDD + FAD)
    double S = (FDD + FAD) / (FDD + FAD + FAA)
    molecule.setParameter("E",E)
    molecule.setParameter("S",S)
  }
    
 })








    






