#@ MoleculeArchive archive

import de.mpg.biochem.mars.molecule.*
import de.mpg.biochem.mars.table.*
import de.mpg.biochem.mars.util.*
import org.scijava.table.*
import org.apache.commons.math3.stat.regression.SimpleRegression


//Set correction factors manually
double beta = 1 //no correction value: 1, reference published value: 0.85
double gamma = 1  //no correction value: 1, reference published value: 1.14
def alpha_list = []
def delta_list = []

archive.metadata().forEach{metadata -> 
	alpha_list.add(metadata.getParameter("alpha"))
	metadata.setParameter("beta", beta)
	metadata.setParameter("gamma", gamma)
	delta_list.add(metadata.getParameter("delta"))
}

double alpha = alpha_list[0]
double delta = delta_list[0]

//Calculate corrected E and S values
archive.getMoleculeUIDs().stream().forEach({UID ->
  def molecule = archive.get(UID)
  if (molecule.hasTag("Active_single")){
  	double Iaemdex = molecule.getParameter("iiIaemdex")
    double Idemdex = molecule.getParameter("iiIdemdex")
    double Iaemaex = molecule.getParameter("iiIaemaex")
    double FAD = Iaemdex - alpha * Idemdex - delta * Iaemaex
    double iiiE = FAD / (FAD + Idemdex)
    double iiiS = (FAD + Idemdex) / (FAD + Idemdex + Iaemaex)
    molecule.setParameter("F_AD",FAD)
    molecule.setParameter("iiiEapp",iiiE)
    molecule.setParameter("iiiSapp",iiiS)
  }
})

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

